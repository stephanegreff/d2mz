#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""d2mz — batch converter Bruker .d (micrOTOF-Q / BAF) -> mzML / mzXML.

Features
  * pure-Python, standard library only (no libbaf2sql, no numpy, no pyteomics)
  * batch: converts every *.d folder found under the input path
  * lock-mass mass recalibration (quadratic ppm-drift model vs RT), with
    automatic per-file polarity detection and the Bruker sodium
    formate/acetate cluster series built in (nothing to type in)
  * Bruker-grade safety: a recalibration is trusted only with >= 7 lock
    observations (the number of calibration points Bruker DataAnalysis
    uses); if the calibrant is weak or absent in a run, the lock-mass
    model is borrowed from the adjacent calibrated file of the same
    series/polarity, and this is stated in the report
  * precursor re-phasing: each MS2 precursor m/z is re-aligned on the
    most intense recalibrated MS1 apex inside the isolation window
    (+-0.5 Da, automatic — no tolerance to tune; Garayev-type correction)
  * run start: the mzML <run> startTimeStamp carries the acquisition
    start (SampleInfo.xml CreationDateTime, converted to UTC;
    fallback: mtime of the .baf); mzXML 3.2 has no standard slot for
    it, so the same date is stamped as an XML comment there
  * output mzML 1.1.0 (indexedmzML) and/or mzXML 3.2, 32- or 64-bit
    binary arrays, optional zlib compression, Bruker centroid data only
  * UV/DAD data of the .d, when present: the FIRST acquired UV channel
    (the first acquired channel above 250 nm) of the file -- or the
    --uv-trace wavelength(s) -- is written
    as an mzML chromatogram (MS-only tools ignore them; mzXML stays free of
    them), and the FULL wavelength spectrum is kept only at the
    chromatographic apexes -- as human-readable mzML run-level
    userParams, so the mzML is self-contained (the .d never needs to
    be shared)
  * retention-time standards: when a 'parabens' run (.d whose name
    contains 'parabens') is part of the folder, the retention anchors of
    the standard series (GallAc, pHOHAc, Me->HepParaben,
    RI = 100 x carbons) are measured on it automatically and stamped
    into every mzML of the folder; manual --parabene name=rt entries
    add or override. The whole homologous series is integrated even
    when the current gradient cannot elute the last members. Stored as
    mzML run-level userParams (transparently
    ignored by MS software); report.json carries processing statistics
    only, no measurement data
  * per-file processing report (JSON) + precursor table (TSV)
  * two validated Impact-II maxima-block layouts, auto-detected by
    geometry -- MS1 blocks are 16n ([idx f64][intensity f32][SNR f32],
    SNR floor filter) while MS2/lock blocks are 12n ([idx f64]
    [intensity f32] with NO SNR array). MS2 conversion also
    declares its acquisition mode (d2mz:acquisition_mode userParam)
    and WARNS when the run is outside the validated scope (DDA CID,
    constant CE or CE ramps -- IsCID, targeted MS/MS lists and
    non-CID fragmentation modes are flagged, never silently
    converted)

Usage
  python3 d2mz.py INPUT_DIR [INPUT_DIR2 ...] \
      --format both --precision 32
  Nothing to set: for every .d file, d2mz automatically detects the
  polarity and locks on the built-in sodium formate/acetate cluster
  series of the standard Bruker calibrant (38 positive / 36 negative
  cluster ions, exact masses validated against the .baf calibrant
  tables). --polarity only forces a polarity when auto-detection is
  not wanted; there is no option to supply another calibrant.

Derived from public reverse engineering of the Bruker BAF format
(H3Jones/bruker-baf, itself derived from the R package 'rtms' by N. Twarog)
and from the baf2sql knowledge of the kbaf2mzml project (bjonnh).

Licence: GPL-3.0 (see LICENSE).
"""
import argparse
import base64
import bisect
import gzip
import hashlib
import json
import math
import os
import re
import sqlite3
import struct
import sys
import time
import zlib
from xml.sax.saxutils import escape, quoteattr

__version__ = "2026-11-25_cp"

LAST_RUN = {'n_warnings': 0, 'warnings': [],
           'n_errors': 0, 'errors': []}   # set by main()

# Bruker DataAnalysis fits its internal lock-mass model on 7 calibration
# points: d2mz keeps the same minimum before trusting a recalibration.
MIN_LOCK_POINTS = 7

# Retention-index standard series, injected in the same sequences as the
# samples ("parabens" runs). RI = 100 x number of carbons (gallic acid C6
# -> 600, ... heptylparaben C14 -> 1400). The WHOLE homologous series up
# to heptyl is integrated: anchors absent
# from the run (short gradient, low ACN) are simply not matched and not
# stamped -- whatever the chromatographic method, a wider gradient to
# 100% ACN will find them with no code change. Anchor ions: [M-H]- in
# negative, [M+H]+ in positive.
PARABEN_STANDARDS = [
    # name       [M-H]-     [M+H]+     RI
    ('GallAc',    169.0143, 171.0289, 600),
    ('pHOHAc',    137.0245, 139.0391, 700),
    ('MeParaben', 151.0401, 153.0547, 800),
    ('EtParaben', 165.0557, 167.0703, 900),
    ('PrParaben', 179.0714, 181.0860, 1000),
    ('BuParaben', 193.0870, 195.1016, 1100),
    ('PeParaben', 207.1027, 209.1173, 1200),
    ('HexParaben', 221.1183, 223.1330, 1300),
    ('HepParaben', 235.1340, 237.1486, 1400),
]
PARABEN_PPM = 15.0   # MS1 match window (recalibrated: plenty)

# Precursor re-phasing search window (Da, half-width). The quadrupole can
# only have isolated the true precursor within its isolation window
# around the setpoint; +-0.5 Da is the widest common window, so it is a
# hard physical cap — no user-tuned tolerance. Within the window the most
# intense recalibrated MS1 apex wins (deterministic; the quadrupole
# isolated the most abundant ion of that window).
PREC_WINDOW_DA = 0.5

# Bruker CollisionSweep acquisitions record, per MS2 spectrum, the TOP of
# the stepped collision-energy window (PerSpectrumVariables Variable 5).
# The windows below were validated empirically on files
# POS10-20_16977.d / POS20-40_17092.d / POS30-60_17100.d (Variable 5 =
# 20.0 / 40.0 / 60.0 respectively, matching the file names).
CE_WINDOWS = {10.0: (5.0, 10.0), 20.0: (10.0, 20.0), 40.0: (20.0, 40.0),
              60.0: (30.0, 60.0), 80.0: (40.0, 80.0),
              # the 2026-11-12 batch carries the SAME /2 sweep
              # pattern at 100-200 eV tops (CE_Act = ramp TOP,
              # validated: 20->10/20 ... 200->100/200)
              100.0: (50.0, 100.0), 120.0: (60.0, 120.0),
              140.0: (70.0, 140.0), 160.0: (80.0, 160.0),
              180.0: (90.0, 180.0), 200.0: (100.0, 200.0)}


def ce_window(ce):
    """Return 'lo-hi' (eV) for the recorded CollisionSweep setpoint, or None."""
    w = CE_WINDOWS.get(ce)
    return '%g-%g' % w if w else None


def method_ce_sweep(d_folder):
    """CE sweep mode DECLARED in the acquisition .method : is the
    per-scan collision energy a SWEEP (ramp lo->hi; the per-scan
    Collision_Energy_Act records the ramp TOP) or a CONSTANT? Also the
    sweep factor and the programmed CE lists (MSMS_ListCollisionEnergy
    / MSMSManual_ListCollisionEnergy). Returns a dict with keys among
    'sweep' (0/1), 'stepping' (0/1), 'sweep_factor', 'step_factor',
    'list_ce', 'manual_list_ce', 'sweep_offsets'."""
    mf = find_method_folder(d_folder)
    if not mf:
        return {}
    out = {}
    for fn in sorted(os.listdir(mf)):
        if not fn.endswith('.method') or 'acquisition' not in fn.lower():
            continue
        txt = None
        for enc in ('utf-8', 'latin-1'):
            try:
                with open(os.path.join(mf, fn), encoding=enc) as f:
                    txt = f.read()
                break
            except (UnicodeDecodeError, OSError):
                continue
        if txt is None:
            continue
        for perm, key in (
                ('CollisionSweep_Switch_Set', 'sweep'),
                ('MSMSAuto_SwitchCollisionEnergyStepping', 'stepping'),
                ('CollisionSweep_CollisionEnergy_Factor_Set',
                 'sweep_factor'),
                ('MSMSAuto_CollisionEnergySteppingFactor', 'step_factor')):
            m = re.search(r'<para[^>]*permname="%s"[^>]*>' % perm, txt)
            if not m:
                continue
            mv = re.search(r'value="([-0-9.eE+]+)"', m.group(0))
            if mv:
                out[key] = float(mv.group(1))
            else:      # value in the element text after the open tag
                mt = re.match(r'\s*([-0-9.eE+]+)\s*<',
                              txt[m.end():m.end() + 200])
                if mt:
                    out[key] = float(mt.group(1))
        for perm, key in (
                ('MSMS_ListCollisionEnergy', 'list_ce'),
                ('MSMSManual_ListCollisionEnergy', 'manual_list_ce'),
                ('Collision_ListCollisionSweep_Energy_Offset',
                 'sweep_offsets')):
            i = txt.find('permname="%s"' % perm)
            if i < 0:
                continue
            end = txt.find('</para_vec', i)
            chunk = txt[i:end] if end > 0 else txt[i:i + 2000]
            vals = [float(v) for v in re.findall(
                r'entry_double value="([0-9.eE+-]+)"', chunk)]
            if vals:
                out.setdefault(key, []).extend(vals)
    return out


def manual_list_targeted(prec_table, min_scans=100, min_repeats=50,
                         min_share=0.8):
    """v_co helper -- do the ACQUIRED MS2 precursors look like a REAL
    manual-target list? A manual list replays the exact same double
    every cycle: thousands of bit-identical declared values. Auto-DDA
    recomputes the precursor from each MS1 survey, so the value
    jitters cycle to cycle and never repeats exactly. True when at
    least min_scans MS2 were acquired and values repeated
    min_repeats+ times each make up at least min_share of all MS2."""
    if len(prec_table) < min_scans:
        return False
    pv = {}
    for r in prec_table:
        k = round(r['declared'], 9)   # exact double: a manual list
        pv[k] = pv.get(k, 0) + 1      # replays bit-identical values
    repeated = sum(c for c in pv.values() if c >= min_repeats)
    return repeated >= min_share * len(prec_table)


def ce_label(ce, ceinfo=None):
    """Unambiguous CE label for one MS2 scan. The mzML cvParam
    'collision energy' carries the Collision_Energy_Act value, which
    is the ramp TOP for a CollisionSweep run and the constant CE
    otherwise -- indistinguishable downstream. This label states the
    mode explicitly: 'ramp lo-hi eV (CollisionSweep, factor F)' or
    'constant X eV'; legacy factor-0.5 table guess when the .method
    declares no switch."""
    if ce is None:
        return None
    sw = (ceinfo or {}).get('sweep')
    if sw is not None:
        if sw:
            f = (ceinfo or {}).get('sweep_factor') or 0.5
            return 'ramp %g-%g eV (CollisionSweep, factor %g)' % (
                round(ce * f, 1), ce, f)
        return 'constant %g eV' % ce
    w = ce_window(ce)
    return 'ramp %s eV (factor 0.5 assumed)' % w if w else None

__date__ = "2026-11-09"

# ============================================================ BAF reader (pure)

# Default calibrant: sodium formate/acetate cluster series used by Bruker
# qTOFs (e.g. micrOTOF-Q, Impact II) with the formate/acetic-acid/NaOH
# calibrant solution. Exact masses validated against the calibrant tables
# shipped in the .baf acquisition files (<= 0.04 ppm). Used
# automatically: for every .d file, d2mz observes both series and
# detects the polarity itself (unless --polarity forces it).
BRUKER_FA_AC_NA_CLUSTERS = {
    'positive': [
        (90.97664, '[Na2(HCOO)1(AcO)0]+'),
        (158.96407, '[Na3(HCOO)2(AcO)0]+'),
        (226.95149, '[Na4(HCOO)3(AcO)0]+'),
        (294.93891, '[Na5(HCOO)4(AcO)0]+'),
        (362.92634, '[Na6(HCOO)5(AcO)0]+'),
        (430.91376, '[Na7(HCOO)6(AcO)0]+'),
        (104.99229, '[Na2(HCOO)0(AcO)1]+'),
        (186.99537, '[Na3(HCOO)0(AcO)2]+'),
        (268.99844, '[Na4(HCOO)0(AcO)3]+'),
        (351.00152, '[Na5(HCOO)0(AcO)4]+'),
        (172.97972, '[Na3(HCOO)1(AcO)1]+'),
        (240.96714, '[Na4(HCOO)2(AcO)1]+'),
        (308.95456, '[Na5(HCOO)3(AcO)1]+'),
        (376.94199, '[Na6(HCOO)4(AcO)1]+'),
        (444.92941, '[Na7(HCOO)5(AcO)1]+'),
        (512.91684, '[Na8(HCOO)6(AcO)1]+'),
        (580.90426, '[Na9(HCOO)7(AcO)1]+'),
        (648.89168, '[Na10(HCOO)8(AcO)1]+'),
        (716.87910, '[Na11(HCOO)9(AcO)1]+'),
        (784.86652, '[Na12(HCOO)10(AcO)1]+'),
        (852.85394, '[Na13(HCOO)11(AcO)1]+'),
        (254.98279, '[Na4(HCOO)1(AcO)2]+'),
        (322.97022, '[Na5(HCOO)2(AcO)2]+'),
        (390.95764, '[Na6(HCOO)3(AcO)2]+'),
        (458.94506, '[Na7(HCOO)4(AcO)2]+'),
        (526.93249, '[Na8(HCOO)5(AcO)2]+'),
        (594.91991, '[Na9(HCOO)6(AcO)2]+'),
        (662.90733, '[Na10(HCOO)7(AcO)2]+'),
        (730.89476, '[Na11(HCOO)8(AcO)2]+'),
        (798.88218, '[Na12(HCOO)9(AcO)2]+'),
        (866.86960, '[Na13(HCOO)10(AcO)2]+'),
        (934.85703, '[Na14(HCOO)11(AcO)2]+'),
        (336.98587, '[Na5(HCOO)1(AcO)3]+'),
        (404.97329, '[Na6(HCOO)2(AcO)3]+'),
        (472.96071, '[Na7(HCOO)3(AcO)3]+'),
        (540.94814, '[Na8(HCOO)4(AcO)3]+'),
        (608.93556, '[Na9(HCOO)5(AcO)3]+'),
        (676.92298, '[Na10(HCOO)6(AcO)3]+'),
        (744.91041, '[Na11(HCOO)7(AcO)3]+'),
        (812.89783, '[Na12(HCOO)8(AcO)3]+'),
        (880.88525, '[Na13(HCOO)9(AcO)3]+'),
        (948.87268, '[Na14(HCOO)10(AcO)3]+'),
        (418.98894, '[Na6(HCOO)1(AcO)4]+'),
    ],
    'negative': [
        (112.98563, '[Na1(HCOO)2(AcO)0]-'),
        (180.97305, '[Na2(HCOO)3(AcO)0]-'),
        (248.96047, '[Na3(HCOO)4(AcO)0]-'),
        (316.94790, '[Na4(HCOO)5(AcO)0]-'),
        (384.93532, '[Na5(HCOO)6(AcO)0]-'),
        (452.92274, '[Na6(HCOO)7(AcO)0]-'),
        (520.91016, '[Na7(HCOO)8(AcO)0]-'),
        (588.89758, '[Na8(HCOO)9(AcO)0]-'),
        (656.88500, '[Na9(HCOO)10(AcO)0]-'),
        (724.87242, '[Na10(HCOO)11(AcO)0]-'),
        (792.85984, '[Na11(HCOO)12(AcO)0]-'),
        (860.84726, '[Na12(HCOO)13(AcO)0]-'),
        (928.83468, '[Na13(HCOO)14(AcO)0]-'),
        (996.82210, '[Na14(HCOO)15(AcO)0]-'),
        (141.01693, '[Na1(HCOO)0(AcO)2]-'),
        (223.02000, '[Na2(HCOO)0(AcO)3]-'),
        (305.02307, '[Na3(HCOO)0(AcO)4]-'),
        (127.00128, '[Na1(HCOO)1(AcO)1]-'),
        (194.98870, '[Na2(HCOO)2(AcO)1]-'),
        (262.97612, '[Na3(HCOO)3(AcO)1]-'),
        (330.96355, '[Na4(HCOO)4(AcO)1]-'),
        (398.95097, '[Na5(HCOO)5(AcO)1]-'),
        (466.93839, '[Na6(HCOO)6(AcO)1]-'),
        (534.92582, '[Na7(HCOO)7(AcO)1]-'),
        (209.00435, '[Na2(HCOO)1(AcO)2]-'),
        (276.99177, '[Na3(HCOO)2(AcO)2]-'),
        (344.97920, '[Na4(HCOO)3(AcO)2]-'),
        (412.96662, '[Na5(HCOO)4(AcO)2]-'),
        (480.95404, '[Na6(HCOO)5(AcO)2]-'),
        (548.94147, '[Na7(HCOO)6(AcO)2]-'),
        (616.92889, '[Na8(HCOO)7(AcO)2]-'),
        (684.91631, '[Na9(HCOO)8(AcO)2]-'),
        (752.90374, '[Na10(HCOO)9(AcO)2]-'),
        (820.89116, '[Na11(HCOO)10(AcO)2]-'),
        (888.87859, '[Na12(HCOO)11(AcO)2]-'),
        (956.86601, '[Na13(HCOO)12(AcO)2]-'),
        (291.00742, '[Na3(HCOO)1(AcO)3]-'),
        (358.99485, '[Na4(HCOO)2(AcO)3]-'),
        (426.98227, '[Na5(HCOO)3(AcO)3]-'),
        (494.96969, '[Na6(HCOO)4(AcO)3]-'),
        (562.95712, '[Na7(HCOO)5(AcO)3]-'),
        (630.94454, '[Na8(HCOO)6(AcO)3]-'),
        (698.93197, '[Na9(HCOO)7(AcO)3]-'),
        (373.01050, '[Na4(HCOO)1(AcO)4]-'),
        (440.99792, '[Na5(HCOO)2(AcO)4]-'),
        (508.98534, '[Na6(HCOO)3(AcO)4]-'),
    ],
}

# --- Impact II (qTOF) quadratic index -> m/z calibration -------------------
# Bruker impact-generation .baf files carry NO mode 4/5 frequency
# calibration record (their .baf_xtr only holds mode-10 transformator
# records). Their centroid/maxima TOF index i maps to m/z through a
# per-file QUADRATIC  m = a*i^2 + b*i + c, which Bruker itself fits on
# the injected calibrant (formate/acetic-acid/NaOH clusters) and which we
# refit here per file from the MS1 maxima (see fit_impact_calibration).
# Validated on the Feb-2025 pollen runs: 0.25 ppm median / 0.38 ppm sigma /
# 0.92 ppm max residual against 30 calibrant targets (Saule NEG run
# 22456, calibrant segment). The bootstraps below are coarse quadratics
# (same hardware family) used ONLY to assign the first calibrant targets
# to raw indexes; the final model is always refitted per file.
IMPACT_QUAD_BOOTSTRAP = {
    # m/z = A*i^2 + B*i + C on the RAW centroid index (i up to ~350000)
    'positive': (6.15294136e-09, 0.00110772993, 49.8417111),
    'negative': (6.174825e-09, 0.001109141, 49.83071),
}
IMPACT_FIT_WINDOWS = (600.0, 200.0, 70.0)   # assignment windows (index units)
IMPACT_TARGET_MIN_FRAC = 0.20  # target must appear in >=20% of scans
IMPACT_IDX_MIN, IMPACT_IDX_MAX = 1500.0, 353000.0


def get_raw_vectors(path):
    out = []
    with open(path, 'rb') as f:
        while True:
            if len(f.read(4)) < 4:
                break
            _ = f.read(4)
            bs = struct.unpack('<i', f.read(4))[0]
            if bs == 0:
                break
            if bs < 24:
                raise ValueError('bad vector size %d' % bs)
            payload = f.read(bs - 24)
            _ = f.read(12)
            out.append(payload)
    return out


def read_index(idx_path):
    sections, header_off = {}, None
    for p in (x for x in get_raw_vectors(idx_path) if len(x) >= 24):
        id_a, id_b = p[0], p[1]
        flag = struct.unpack_from('<i', p, 8)[0]
        off = struct.unpack_from('<i', p, 16)[0]
        page = struct.unpack_from('<i', p, 20)[0]
        abs_off = (page << 32) + (off & 0xFFFFFFFF)
        if (id_a, id_b) == (1, 48):
            header_off = abs_off
        elif (id_a, id_b) == (2, 32):
            sections.setdefault(flag, {})['metadata'] = abs_off
        elif (id_a, id_b) == (1, 16):
            sections.setdefault(flag, {})['spectrum'] = abs_off
        elif (id_a, id_b) == (2, 16):
            sections.setdefault(flag, {})['maxima'] = abs_off
    return header_off, sections


def calibration_candidates(xtr_path):
    """All header-shaped records of the baf_xtr, with their offsets."""
    out = []
    for p in get_raw_vectors(xtr_path):
        if len(p) < 8:
            continue
        if struct.unpack_from('<i', p, 0)[0] != len(p):
            continue
        hsize = struct.unpack_from('<i', p, 4)[0]
        if hsize < 0:
            continue
        bso = 4 + hsize
        if bso + 4 > len(p):
            continue
        bsize = struct.unpack_from('<i', p, bso)[0]
        if bsize < hsize:
            continue
        fo = bso + bsize - hsize
        if fo + 64 > len(p):
            continue
        out.append((p, fo))
    return out


def probe_calibration(xtr_path, log=print, full=True):
    """Diagnostic dump: describe every calibration-record candidate so the
    real Bruker layout can be identified from a --probe run."""
    cands = calibration_candidates(xtr_path)
    log('probe (d2mz %s): %d calibration-record candidate(s) in %s'
        % (__version__, len(cands), os.path.basename(xtr_path)))
    for k, (p, fo) in enumerate(cands[:5]):
        hsize = struct.unpack_from('<i', p, 4)[0]
        bsize = struct.unpack_from('<i', p, 4 + hsize)[0]
        ints = struct.unpack_from('<8i', p, fo)
        dbls = struct.unpack_from('<8d', p, fo)
        log('  cand %d: len=%d hsize=%d bsize=%d fo=%d' %
            (k, len(p), hsize, bsize, fo))
        log('    ints@fo : %s' % list(ints))
        log('    dbls@fo : %s' % ['%.6g' % d for d in dbls])
        for off in (fo + 8, fo + 12, fo + 16):
            if off + 8 <= len(p):
                log('    dbls@%d: %s' % (off, ['%.6g' % d for d in
                    struct.unpack_from('<2d', p, off)]))
    if not full:
        return cands
    # full dump of the first 2 records of each distinct shape
    by_shape = {}
    for k, (p, fo) in enumerate(cands):
        hsize = struct.unpack_from('<i', p, 4)[0]
        bsize = struct.unpack_from('<i', p, 4 + hsize)[0]
        by_shape.setdefault((len(p), hsize, bsize), []).append(k)
    log('probe: %d distinct shape(s): %s'
        % (len(by_shape), sorted(by_shape.keys())))
    for key, ks in sorted(by_shape.items()):
        for k in ks[:2]:
            p, fo = cands[k]
            log('== FULL DUMP cand %d (len=%d, fo=%d) ==' % (k, len(p), fo))
            log('  header ints @0..%d: %s' % (fo - 1,
                list(struct.unpack_from('<%di' % (fo // 4), p, 0))
                if fo % 4 == 0 else
                ' '.join('%02x' % b for b in p[:fo])))
            n_tail = len(p) - fo
            log('  tail doubles @fo (%d B): %s' % (n_tail,
                ['%.6g' % d for d in struct.unpack_from(
                    '<%dd' % (n_tail // 8), p, fo)]
                if n_tail % 8 == 0 else
                ' '.join('%02x' % b for b in p[fo:])))
            log('  tail ints @fo (%d B): %s' % (n_tail,
                list(struct.unpack_from('<%di' % (n_tail // 4), p, fo))))
    return cands


TAG_ACQUISITION = 0xBFA03001
TAG_STRINGS = 0xBFA03003
TAG_CALIBRATION = 0xBFA00006


def _probe_xtr_overview(xtr_path, log=print):
    """Whole .baf_xtr census: raw-vector count and shapes, plus the
    per-record header/lock table of the calibration candidates."""
    vecs = get_raw_vectors(xtr_path)
    counts = {}
    for p in vecs:
        declared = struct.unpack_from('<i', p, 0)[0] if len(p) >= 4 else -1
        counts[(len(p), declared)] = counts.get((len(p), declared), 0) + 1
    log('probe: xtr census: %d raw vector(s), (len, declared) -> count: %s'
        % (len(vecs), sorted(counts.items())))
    cands = calibration_candidates(xtr_path)
    log('probe: per-record table (first 40 candidates):')
    for k, (p, fo) in enumerate(cands[:40]):
        hid = struct.unpack_from('<i', p, 20)[0]
        h7 = struct.unpack_from('<i', p, 28)[0]
        h8 = struct.unpack_from('<i', p, 32)[0]
        if len(p) >= fo + 8 and struct.unpack_from('<i', p, fo)[0] == 10:
            lock = '%.6g' % struct.unpack_from('<f', p, fo + 4)[0]
        else:
            lock = '-'
        log('  rec %d: len=%d ID=%d h7=%d h8=%d lock=%s'
            % (k, len(p), hid, h7, h8, lock))


def _probe_blocks(baf, sections, log=print, n_scans=3, cap=40,
                  extra_scans=()):
    """Dump the maxima (peak) blocks and spectrum-block headers of the
    first scans (plus the lock-mass scans), so the index->m/z
    calibration can be identified."""
    allp = [(f, o) for f, o in sorted(sections.items())
            if 'maxima' in o]
    picks = allp[:n_scans]
    for pick in allp:
        if pick[0] in extra_scans and pick not in picks:
            picks.append(pick)
    if not picks:
        log('probe: no maxima section in .baf_idx')
    for flag, offs in picks:
        try:
            off = offs['maxima']
            with open(baf, 'rb') as f:
                f.seek(off)
                bs, _d, hs = struct.unpack('<i4si', f.read(12))
                f.seek(off + hs)
                n, po = struct.unpack('<ii', f.read(8))
                log('probe: scan %s maxima: bs=%d hs=%d n=%d po=%d'
                    % (flag, bs, hs, n, po))
                f.seek(off)
                log('  maxima head hex: %s'
                    % f.read(min(hs + po + 160, bs)).hex())
                tail_n = min(256, bs)
                f.seek(off + bs - tail_n)
                log('  maxima tail hex (last %d B): %s'
                    % (tail_n, f.read(tail_n).hex()))
                f.seek(off + bs - 256)
                tv = struct.unpack('<32d', f.read(256))
                log('  maxima tail as 16 pairs:')
                for j in range(16):
                    log('    %.10g  %.10g' % (tv[2 * j], tv[2 * j + 1]))
                if 0 < n <= 200000 and po >= 8:
                    f.seek(off + hs + po)
                    idx = struct.unpack('<%dd' % n, f.read(n * 8))
                    itn = struct.unpack('<%dd' % n, f.read(n * 8))
                    top = sorted(range(n), key=lambda k: itn[k],
                                 reverse=True)[:cap]
                    log('  top-%d by intensity (idx, int):' % len(top))
                    for k in top:
                        log('    %.10g  %.10g' % (idx[k], itn[k]))
        except (OSError, struct.error) as e:
            log('probe: maxima dump failed for scan %s: %s' % (flag, e))
        so = offs.get('spectrum')
        if so is not None:
            try:
                with open(baf, 'rb') as f:
                    f.seek(so)
                    bs, _d, hs = struct.unpack('<i4si', f.read(12))
                    f.seek(so)
                    log('  spectrum head hex: %s' % f.read(64).hex())
            except (OSError, struct.error) as e:
                log('probe: spectrum dump failed: %s' % e)


def _probe_idx_flags(idx_path, baf, log=print):
    """Whole .baf_idx inventory: every (id_a, id_b) record type with
    entry count, flag range and a hex dump of the first target."""
    groups = {}
    for p in (x for x in get_raw_vectors(idx_path) if len(x) >= 24):
        key = (p[0], p[1])
        flag = struct.unpack_from('<i', p, 8)[0]
        off = struct.unpack_from('<i', p, 16)[0]
        page = struct.unpack_from('<i', p, 20)[0]
        abs_off = (page << 32) + (off & 0xFFFFFFFF)
        g = groups.setdefault(key, [0, abs_off, None, None])
        g[0] += 1
        g[2] = flag if g[2] is None else min(g[2], flag)
        g[3] = flag if g[3] is None else max(g[3], flag)
    log('probe: .baf_idx (id_a,id_b) inventory: %d type(s)'
        % len(groups))
    for key in sorted(groups):
        cnt, off, fmin, fmax = groups[key]
        line = ('  idx %s: %d entries, flags %s..%s, first @%d'
                % (key, cnt, fmin, fmax, off))
        try:
            with open(baf, 'rb') as f:
                f.seek(off)
                line += '\n    hex: %s' % f.read(96).hex()
        except OSError as e:
            line += ' (read failed: %s)' % e
        log(line)


def _probe_acquisition_header(baf, header_off, log=print):
    """Dump the acquisition-header record of analysis.baf: ASCII strings
    (instrument, operator, method...) and every tagged sub-record, with a
    full hex dump of the reference-calibration payloads (0xBFA00006)."""
    if not header_off:
        log('probe: no header section offset in .baf_idx')
        return
    with open(baf, 'rb') as f:
        f.seek(header_off)
        size, tag = struct.unpack('<II', f.read(8))
        if tag != TAG_ACQUISITION:
            log('probe: acquisition-header tag=0x%08x (expected 0x%08x)'
                % (tag, TAG_ACQUISITION))
            return
        f.seek(header_off)
        rec = f.read(size)
    log('probe: acquisition header: %d bytes at %d' % (size, header_off))
    for k, m in enumerate(re.finditer(rb'[\x20-\x7e]{6,}', rec)):
        if k >= 12:
            log('  ... (ascii strings truncated)')
            break
        log('  ascii @%d: %r' % (m.start(), m.group().decode('ascii')))
    seen = 0
    for o in range(8, len(rec) - 8, 4):
        sz, tg = struct.unpack_from('<II', rec, o)
        if (8 <= sz <= len(rec) - o
                and 0xBFA00000 <= tg <= 0xBFA0FFFF):
            log('  subrec @%d: size=%d tag=0x%08x' % (o, sz, tg))
            if tg == TAG_CALIBRATION:
                log('    hex: %s' % rec[o:o + sz].hex())
            seen += 1
            if seen >= 60:
                log('  ... (sub-record list truncated)')
                break


def _probe_sqlite(d_folder, log=print):
    """List the tables of analysis.sqlite and sample the interesting
    ones (properties, calibration, spectra...)."""
    sq = os.path.join(d_folder, 'analysis.sqlite')
    if not os.path.exists(sq):
        log('probe: no analysis.sqlite in %s' % d_folder)
        return
    c = sqlite3.connect('file:%s?mode=ro' % sq.replace('\\', '/'), uri=True)
    try:
        names = [r[0] for r in c.execute(
            "select name from sqlite_master where type in ('table','view') "
            "order by name")]
        log('probe: sqlite tables: %s' % names)
        for name in names:
            if name.startswith('sqlite_'):
                continue
            try:
                n = c.execute('select count(*) from "%s"' % name
                              ).fetchone()[0]
            except sqlite3.Error:
                continue
            log('  table %s: %d rows' % (name, n))
            if re.search(r'transform|variab|support|referenc|perseg',
                         name, re.I):
                cols = [r[1] for r in
                        c.execute('pragma table_info("%s")' % name)]
                log('    cols: %s' % cols)
                try:
                    cap = 100 if re.search(
                        r'transform|variab|support|referenc|perseg',
                        name, re.I) else 3
                    for r in c.execute(
                            'select * from "%s" limit %d' % (name, cap)):
                        log('    row: %s'
                            % [(v if not isinstance(v, bytes)
                                else v[:40].hex()) for v in r])
                except sqlite3.Error as e:
                    log('    (row dump failed: %s)' % e)
        # spectra line-record ids, transformator decode, misc tables
        try:
            cols = [r[1] for r in
                    c.execute('pragma table_info(Spectra)')]
            log('  Spectra cols: %s' % cols)
            nn = {}
            for col in cols:
                if col != 'Id' and col.endswith('Id'):
                    nn[col] = c.execute(
                        'select count(*) from "Spectra" '
                        'where "%s" is not null' % col).fetchone()[0]
            log('  Spectra non-null per Id column: %s' % nn)
            row = c.execute(
                'select * from Spectra where LineMzId is not null '
                'limit 1').fetchone()
            if row is not None:
                log('  Spectra first row with LineMzId: %s' % list(row))
            idcol = 'Id' if 'Id' in cols else cols[0]
            for sid in (25, 26, 35, 36, 37, 38):
                r = c.execute('select * from Spectra where "%s"=?'
                              % idcol, (sid,)).fetchone()
                if r is not None:
                    log('  Spectra scan %s: %s' % (sid, list(r)))
        except sqlite3.Error as e:
            log('  (Spectra detail failed: %s)' % e)
        try:
            rows = {}
            for r in c.execute('select Id, Blob from Transformators'):
                rows[r[0]] = r[1]
            log('  Transformators: %d rows, %d distinct blob(s)'
                % (len(rows), len(set(rows.values()))))

            def _dec(b):
                b = bytes(b)
                i0 = struct.unpack_from('<i', b, 0)[0]
                dd = [struct.unpack_from('<d', b, o)[0]
                      for o in (4, 12, 20, 28)]
                tl = struct.unpack_from('<i', b, 36)[0]
                return (i0, ['%.12g' % x for x in dd], tl)
            for sid in (1, 25, 26, 34, 35, 36, 37, 38, 39, 40):
                if sid in rows:
                    log('    T[%d] = %s' % (sid, _dec(rows[sid])))
            for rid, ak, tid in c.execute(
                    'select Id, AcquisitionKey, TransformatorId '
                    'from ReferenceTransformators'):
                if tid in rows:
                    log('    T[ref id=%d ak=%s] = %s'
                        % (tid, ak, _dec(rows[tid])))
        except (sqlite3.Error, struct.error) as e:
            log('  (Transformators decode failed: %s)' % e)
        for nm in ('AcquisitionKeys', 'Info', 'Properties', 'Steps'):
            try:
                for r in c.execute('select * from "%s" limit 12' % nm):
                    log('  %s row: %s' % (nm,
                        [v[:48].hex() if isinstance(v, bytes) else v
                         for v in r]))
            except sqlite3.Error:
                pass
    finally:
        c.close()


def _read_baf_params(baf, header_off):
    """Acquisition-header parameter table: {row_index: (type, fields)}.
    Ported from bruker-baf (H3Jones) read_baf_file_header."""
    params = {}
    with open(baf, 'rb') as f:
        f.seek(header_off)
        bs = struct.unpack('<i', f.read(4))[0]
        _ = f.read(8)                       # delimiter + subdelimiter
        f.read(4)                            # header size
        nxt = struct.unpack('<i', f.read(4))[0]
        if bs < 20 or nxt < 20:
            raise ValueError('bad acquisition header (%d/%d)' % (bs, nxt))
        f.seek(header_off + nxt)
        mbs = struct.unpack('<i', f.read(4))[0]
        f.seek(header_off + nxt + mbs)
        tbs = struct.unpack('<i', f.read(4))[0]
        _ = f.read(6)
        nrow = struct.unpack('<H', f.read(2))[0]
        for _r in range(nrow):
            rs = struct.unpack('<H', f.read(2))[0]
            if rs < 20:
                raise ValueError('bad header row size %d' % rs)
            row = b'\x00\x00' + f.read(rs - 2)
            ri, rt = struct.unpack_from('<Hi', row, 2)
            offs = struct.unpack_from('<6H', row, 8)
            fields = [None] * 6
            for i, o in enumerate(offs):
                if o == 0:
                    continue
                nxt_f = next((x for x in offs[i + 1:] if x), rs)
                fields[i] = row[o:nxt_f].split(b'\x00')[0]
            params[ri] = (rt, fields)
    return params


def _read_metadata_values(baf, off, params, log=print, known=None):
    """Per-scan metadata block: the VALUES of the acquisition-header
    parameters for one scan (numeric + string rows). Ported from
    bruker-baf read_baf_metadata_block. Returns {name: raw_value}.
    known is None -> log every value (and fill the dict in place);
    known is a dict -> log only the values that changed vs that dict."""
    values = {}

    def note(ri, raw):
        entry = params.get(ri)
        name = ('row%d' % ri) if entry is None else (
            entry[1][0].decode('utf-8', 'replace') if entry[1][0]
            else 'row%d' % ri)
        values[name] = raw
        if known is not None:
            if known.get(name) != raw:
                log('    %s: %r -> %r' % (name, known.get(name), raw))
            known[name] = raw
            return
        label_val = raw
        if entry is not None:
            vals = entry[1][3]
            if vals and isinstance(raw, float):
                try:
                    for item in vals.decode('utf-8', 'replace').split(';'):
                        code, sep, lab = item.partition(':')
                        if sep and float(code.strip()) == round(raw):
                            label_val = '%s (%r)' % (lab, raw)
                            break
                except ValueError:
                    pass
        log('    %s = %r' % (name, label_val))

    with open(baf, 'rb') as f:
        f.seek(off)
        mbs, _md, mhs = struct.unpack('<i4si', f.read(12))
        if mhs < 12 or mhs > mbs:
            raise ValueError('bad metadata block sizes %d/%d' % (mbs, mhs))
        f.seek(off + mhs)
        ns = f.tell()
        nbs, _nd, nhs, ncnt = struct.unpack('<i4sii', f.read(16))
        if nhs < 16 or nhs > nbs:
            raise ValueError('bad numeric metadata header')
        f.seek(ns + nhs)
        for _i in range(ncnt):
            # the row TYPE of the acquisition-header table says
            # how to read the 8 value bytes -- codes ending in 0x05
            # (261/517/65285/130821...) are INTEGER rows, 0x0A are
            # doubles. The denormals (SegmentId=5e-324,
            # Capillary=1.16e-320...) were those integers read as
            # f64; reinterpreting them fixes SegmentId=1 etc.
            ri, raw8 = struct.unpack('<H8s', f.read(10))
            num = struct.unpack('<d', raw8)[0]
            entry = params.get(ri)
            if entry is not None and (entry[0] & 0xFF) == 5:
                num = struct.unpack('<q', raw8)[0]
            note(ri, num)
        f.seek(ns + nbs)
        sbs, _sd, shs, scnt = struct.unpack('<i4sii', f.read(16))
        st = f.tell()
        f.seek(st + shs)
        for _i in range(scnt):
            rs = struct.unpack('<H', f.read(2))[0]
            if rs < 8:
                log('    [v8 metadata: short string row rs=%d, stop]'
                    % rs)
                break
            row = b'\x00\x00' + f.read(rs - 2)
            ri = struct.unpack_from('<H', row, 2)[0]
            dlen = struct.unpack_from('<i', row, 4)[0]
            slen = dlen - 1
            if slen < 0 or slen > rs - 8:
                log('    [v8 metadata: odd string row %d rs=%d dlen=%d '
                    'hex=%s]' % (ri, rs, dlen, row[:40].hex()))
                slen = max(0, min(slen, rs - 8))
            txt = row[8:8 + slen].decode('utf-8', 'replace')
            note(ri, txt)
    return values


def _maxima_layout(bs, hs, n, po):
    """classify a maxima block by its geometry. Two families were
    VALIDATED byte-exact on a real 2026-11-06 batch
    (A1/B1 Impact II, probe dumps 2026-11-12):
      - 3-array '16n' (MS1 scans, metadata len 412):
        body = bs-hs-po = 16n = [idx f64 8n][INTENSITY f32 4n][SNR f32 4n]
        (scan 1 of 37320: 52+40+16*8925 = 142892 = bs, exact)
      - 2-array '12n' (MS2 & lock-mass scans, metadata len 560):
        body = 12n = [idx f64 8n][INTENSITY f32 4n] and NO SNR ARRAY
        (scan 23: 52+40+12*972 = 11756 = bs; scan 24: 12176; scan 25:
        11120 -- all exact). read_maxima used to read a 'SNR'
        array at +12n i.e. PAST THE END of every MS2 block -- random
        bytes of the next section -- then filtered peaks on ita<2,
        silently dropping ~25-50% of the peaks of EVERY MS2 scan,
        intense fragments included (the 'missing 317.066' symptom).
    Returns 'A3' (16n), 'A2' (12n) or None (block too small)."""
    body = bs - hs - po
    if n <= 0 or body <= 0:
        return None
    if body == 16 * n:
        return 'A3'
    if body == 12 * n:
        return 'A2'
    # best effort for irregular blocks: largest layout that still fits
    if body > 16 * n:
        return 'A3'
    if body > 12 * n:
        return 'A2'
    return None


def _read_maxima_raw(baf, off, f32_ints=False):
    """Validated maxima-block layout (rtms / H3Jones bruker-baf):
    [i32 bs][4s tag][i32 hs] header, then at off+hs: [i32 n][i32 po],
    n index doubles at off+hs+po, then the intensity arrays.

    f32_ints=False: micrOTOF layout — n intensity doubles (rtms).

    f32_ints=True: Impact II layout (RE-PROBED BYTE-BY-BYTE
    2026-11-07): [n idx f64][n INTENSITY f32][n SNR f32]
    as THREE CONSECUTIVE arrays. Validated on real calibrant blocks
    (GA/RA 2Hz runs 2026 + Feb-2025 1uL RA runs + Taraxacum 4161):
    intensity 20..4.6e6 (calibrant clusters are the STRONGEST peaks),
    SNR 1..16. The former interleaved [field][SNR]-per-peak reading
    mixed the two arrays and fed the calibration fit garbage weights
    (garbage fits at 45-295 ppm). The intensity
    is returned in 'itn' (used everywhere as the peak weight); the SNR
    array is kept in 'ita' for inspection."""
    with open(baf, 'rb') as f:
        f.seek(off)
        bs, _tag, hs = struct.unpack('<i4si', f.read(12))
        if bs <= 0 or hs < 12 or hs > bs:
            return None
        f.seek(off + hs)
        n, po = struct.unpack('<ii', f.read(8))
        if n < 0 or po < 8 or hs + po > bs:
            return None
        lay = _maxima_layout(bs, hs, n, po)
        if lay is None:
            return None
        f.seek(off + hs + po)
        idx = struct.unpack('<%dd' % n, f.read(n * 8))
        if f32_ints:
            # Impact: intensity is f32; SNR array only in the 16n layout
            itn = struct.unpack('<%df' % n, f.read(n * 4))
            ita = (struct.unpack('<%df' % n, f.read(n * 4))
                   if lay == 'A3' else None)
        elif lay == 'A2':
            itn = struct.unpack('<%df' % n, f.read(n * 4))
            ita = None
        else:
            itn = struct.unpack('<%dd' % n, f.read(n * 8))
            ita = None
    return dict(bs=bs, hs=hs, n=n, po=po, idx=idx, itn=itn, ita=ita,
                layout=lay, off=off)


def _dec_transformator(blob):
    b = bytes(blob)
    if len(b) < 40:
        return None
    i0 = struct.unpack_from('<i', b, 0)[0]
    dd = [struct.unpack_from('<d', b, o)[0] for o in (4, 12, 20, 28)]
    tl = struct.unpack_from('<i', b, 36)[0]
    return (i0, dd, tl)


def _probe_v7(d, baf, sections, log=print):
    """v7 model-oriented probe.

    a) per-scan transformator blobs from analysis.sqlite (i0, doubles);
    b) Steps.Mass locks per scan;
    c) maxima blocks decoded with the validated layout for the first MS
       scans + lock scans: first/top (idx, intensity) pairs, candidate
       index->m/z model values, intensity sanity stats;
    d) full hex of the small baf_xtr vectors (20/28/76 B) and of the
       first len=560 record;
    e) hex walk of the (4,64) per-spectrum containers (likely home of
       the Line*/Profile* arrays);
    f) Spectra binary-array Ids in hex for a few scans.
    """
    log('probe: ===== v7 calibration-model probe =====')
    xtr = os.path.join(d, 'analysis.baf_xtr')

    # --- sqlite: locks, spectra rows, transformators -------------------
    spectra_rows, steps, locks, transf, spec_ids = [], {}, {}, {}, {}
    try:
        spectra_rows, steps, _psv = load_spectra_table(d, log)
        for sid, sts in steps.items():
            for st in sts:
                m = st.get('Mass')
                if m and m > 0:
                    locks.setdefault(sid, []).append(m)
        locks = {k: sorted(set(v)) for k, v in locks.items()}
        log('probe v7: Steps.Mass per lock scan: %s'
            % {k: ['%.5f' % m for m in v] for k, v in sorted(locks.items())})
    except Exception as e:
        log('probe v7: sqlite Steps load failed: %s' % e)
    try:
        c = sqlite3.connect('file:%s?mode=ro'
                             % os.path.join(d, 'analysis.sqlite')
                             .replace('\\', '/'), uri=True)
        try:
            for r in c.execute('select Id, Blob from Transformators'):
                transf[r[0]] = bytes(r[1])
            log('probe v7: Transformators: %d rows, %d distinct blob(s)'
                % (len(transf), len(set(transf.values()))))
            for sid in (1, 25, 26, 34, 35, 36, 37, 38):
                r = c.execute('select Id, TransformatorId, ProfileMzId, '
                              'ProfileIntensityId, LineMzId, '
                              'LineIntensityId, LineIndexId, LinePeakAreaId, '
                              'LineSnrId from Spectra where Id=?',
                              (sid,)).fetchone()
                if r is not None:
                    spec_ids[r[0]] = r[1]
                    log('  Spectra id=%d Tid=%s ProfMz=%s ProfInt=%s '
                        'LineMz=%s LineInt=%s LineIdx=%s Area=%s Snr=%s'
                        % tuple([r[0]] +
                                ['%s' % (('0x%016x' % v) if v is not None
                                         else '-') for v in r[1:]]))
                    tid = r[1]
                    if tid is not None and tid in transf:
                        log('    T[%d] = %s (hex: %s)'
                            % (tid, _dec_transformator(transf[tid]),
                               transf[tid].hex()))
        finally:
            c.close()
    except Exception as e:
        log('probe v7: sqlite Transformators failed: %s' % e)

    # global fallback transformator values (ref blob, if any)
    tref = None
    for blob in transf.values():
        dec = _dec_transformator(blob)
        if dec and dec[0] == 12:
            tref = dec
            break
    if tref is not None:
        log('probe v7: reference transformator (i0=12): %s' % (tref,))

    # --- maxima blocks with the validated layout ----------------------
    allp = [(f, o) for f, o in sorted(sections.items()) if 'maxima' in o]
    lock_scans = sorted(locks)[:6]
    picks = [p for p in allp if p[0] in lock_scans]
    picks += [p for p in allp if p[0] not in lock_scans][:3]
    for flag, offs in picks:
        try:
            r = _read_maxima_raw(baf, offs['maxima'])
        except (OSError, struct.error) as e:
            log('probe v7: maxima read failed for scan %s: %s' % (flag, e))
            continue
        if r is None:
            log('probe v7: scan %s maxima: block rejected (layout)' % flag)
            continue
        n = r['n']
        itn = r['itn']
        idx = r['idx']
        bad = sum(1 for v in itn if not (v == v) or v in (float('inf'),
                                                          float('-inf')))
        huge = sum(1 for v in itn if v == v and v > 1e7)
        log('probe v7: scan %s maxima (bs=%d n=%d po=%d): '
            'intensities: %d huge(>1e7), %d non-finite'
            % (flag, r['bs'], n, r['po'], huge, bad))
        log('  first 8 (idx, int): %s'
            % '  '.join('%.6g/%.6g' % (idx[k], itn[k])
                        for k in range(min(8, n))))
        order = sorted(range(n), key=lambda k: itn[k], reverse=True)
        top = [k for k in order if itn[k] == itn[k] and itn[k] < 1e25][:10]
        tvals = (tref[1] if tref else None)
        log('  top-10 by intensity (sane):')
        for k in top:
            i = idx[k]
            alts = ['int=%.6g idx=%.6g' % (itn[k], i)]
            if tvals:
                t0, big = tvals[2], tvals[3]
                if i > 0:
                    alts.append('C/i=%.5f' % (big / i))
                if i > t0:
                    alts.append('C/(i-t0)=%.5f' % (big / (i - t0)))
            log('    %s' % '  '.join(alts))
        if flag in locks:
            log('  known locks for scan %s: %s'
                % (flag, ['%.5f' % m for m in locks[flag]]))

    # --- small xtr vectors + first 560-record hex ----------------------
    try:
        for p in get_raw_vectors(xtr):
            if len(p) in (20, 28, 76):
                log('probe v7: xtr vector len=%d: %s' % (len(p), p.hex()))
        first560 = None
        for p in get_raw_vectors(xtr):
            if len(p) == 560:
                first560 = p
                break
        if first560 is not None:
            log('probe v7: first len=560 xtr vector hex (%d B): %s'
                % (len(first560), first560.hex()))
    except Exception as e:
        log('probe v7: xtr small-vector dump failed: %s' % e)

    # --- (4,64) containers: hex walk of the first scans ----------------
    try:
        idx_path = os.path.join(d, 'analysis.baf_idx')
        cont = {}
        for p in (x for x in get_raw_vectors(idx_path) if len(x) >= 24):
            if (p[0], p[1]) == (4, 64):
                flag = struct.unpack_from('<i', p, 8)[0]
                off = struct.unpack_from('<i', p, 16)[0]
                page = struct.unpack_from('<i', p, 20)[0]
                cont[flag] = (page << 32) + (off & 0xFFFFFFFF)
        log('probe v7: (4,64) containers: %d entr%s'
            % (len(cont), 'y' if len(cont) == 1 else 'ies'))
        shown = 0
        for flag in sorted(cont):
            if flag in lock_scans or flag in (25, 26):
                if shown >= 2:
                    break
                shown += 1
                off = cont[flag]
                with open(baf, 'rb') as f:
                    f.seek(off)
                    sz, tg = struct.unpack('<II', f.read(8))
                    log('probe v7: (4,64) scan %s @%d: size=%d tag=0x%08x'
                        % (flag, off, sz, tg))
                    f.seek(off)
                    rec = f.read(min(sz, 4096))
                log('  head hex (first 512 B): %s' % rec[:512].hex())
                o = 0
                walked = 0
                while o + 8 <= len(rec) and walked < 12:
                    s2, t2 = struct.unpack_from('<II', rec, o)
                    if 8 <= s2 <= len(rec) - o and \
                            0xBFA00000 <= t2 <= 0xBFA0FFFF:
                        log('  subrec @%d: size=%d tag=0x%08x hex=%s'
                            % (o, s2, t2, rec[o:o + min(s2, 64)].hex()))
                        o += s2
                        walked += 1
                    else:
                        o += 4
    except Exception as e:
        log('probe v7: (4,64) container walk failed: %s' % e)


def _probe_v8(d, baf, sections, log=print):
    """v8 decisive dump for the index->m/z fit.

    a) full maxima peak lists ((idx, intensity), every peak) of the first
       lock scans; for the first MS scans the strongest peaks with
       intensity < 1e7 (the >1e7 values are descriptor garbage, not
       intensities);
    b) per-scan transformator blob next to each dumped scan;
    c) the (1,16) spectrum-block payload of scan 25 as raw float32
       (head hex + values) to see whether profile m/z is stored there;
    Rows are prefixed 'v8p' for easy offline extraction.
    """
    log('probe: ===== v8 full-peak dump =====')
    locks = {}
    try:
        _rows, steps, _psv = load_spectra_table(d, log)
        for sid, sts in steps.items():
            for st in sts:
                m = st.get('Mass')
                if m and m > 0:
                    locks.setdefault(sid, []).append(m)
        locks = {k: sorted(set(v)) for k, v in locks.items()}
    except Exception as e:
        log('probe v8: Steps load failed: %s' % e)
    transf, spec_tid = {}, {}
    try:
        c = sqlite3.connect('file:%s?mode=ro'
                            % os.path.join(d, 'analysis.sqlite')
                            .replace('\\', '/'), uri=True)
        try:
            for r in c.execute('select Id, Blob from Transformators'):
                transf[r[0]] = bytes(r[1])
            for r in c.execute('select Id, TransformatorId from Spectra'):
                spec_tid[r[0]] = r[1]
        finally:
            c.close()
    except Exception as e:
        log('probe v8: transformator load failed: %s' % e)
    picks = [f for f in sorted(locks)[:8]
             if f in sections and 'maxima' in sections[f]]
    picks += [f for f in (25, 26, 27)
              if f not in locks and f in sections
              and 'maxima' in sections[f]][:2]
    for flag in picks:
        try:
            r = _read_maxima_raw(baf, sections[flag]['maxima'])
        except (OSError, struct.error) as e:
            log('probe v8: scan %s maxima failed: %s' % (flag, e))
            continue
        if r is None:
            log('probe v8: scan %s maxima: block rejected (layout)'
                % flag)
            continue
        n, idx, itn = r['n'], r['idx'], r['itn']
        tid = spec_tid.get(flag)
        tb = _dec_transformator(transf[tid]) if tid in transf else None
        log('probe v8: scan %s maxima n=%d lock=%s T=%s'
            % (flag, n, ['%.5f' % m for m in locks.get(flag, [])], tb))
        if flag in locks:
            rows = list(range(n))
        else:
            rows = sorted((k for k in range(n) if itn[k] < 1e7),
                          key=lambda k: itn[k], reverse=True)[:1500]
        for k in rows:
            log('  v8p %d %.4f %.8g' % (flag, idx[k], itn[k]))
    # spectrum-block payload of scan 25
    o = sections.get(25, {}).get('spectrum')
    if o:
        try:
            with open(baf, 'rb') as f:
                f.seek(o)
                bs, tg, hs = struct.unpack('<i4si', f.read(12))
                f.seek(o)
                data = f.read(min(bs, 2048))
            log('probe v8: scan 25 spectrum block bs=%d hs=%d '
                'payload_f32=%d' % (bs, hs, (bs - hs) // 4))
            log('  head hex after header (256 B): %s'
                % data[hs:hs + 256].hex())
            nfv = (len(data) - hs) // 4
            if nfv > 0:
                vals = struct.unpack_from('<%df' % nfv, data, hs)
                log('  first 64 f32: %s'
                    % ' '.join('%.6g' % v for v in vals[:64]))
        except (OSError, struct.error) as e:
            log('probe v8: spectrum block dump failed: %s' % e)


def _probe_sweep_ms2(d, baf, header_off, sections, log=print, slog=None):
    """targeted, COMPACT dump of the MS2 maxima blocks, written
    to a SEPARATE small file (d2mz_probe_sweep_<name>.txt) because the
    full probe report exceeds the 4 MiB upload limit of a
    chat. Two goals, in ONE round-trip:
    (a) confirm the maxima-block layout -- MS2/lock maxima blocks are
        12n two-array [idx f64][INTENSITY f32] with NO SNR array, so
        an earlier layout reader read the 'SNR' past the block end and filtered
        intense fragments away at random;
    (b) settle the remaining question: the INTENSITY SCALE of the
        high-CE sweep scans (DataAnalysis shows the base fragment
        317.066 at ~3.5e5 on EVERY ramp; the converted mzML only
        reaches that scale on the low ramp). Every MS2 scan whose
        precursor falls in the watch window is dumped with the top
        peaks and enough raw hex to re-interpret the arrays offline."""
    slog = slog or log
    def quiet(*_a, **_k):
        pass
    slog('probe (d2mz %s): sweep-MS2 dump -- folder %s'
         % (__version__, d))
    cei = method_ce_sweep(d)
    slog('probe: method CE info: %r' % cei)
    try:
        spectra_rows, steps, psv = load_spectra_table(d, quiet)
    except Exception as e:
        slog('probe: spectra table failed (%s) -- no sweep dump' % e)
        return
    # watch window: flavonoid glycosides (A: 479.118 /
    # B: 641.172); fall back to the first real MS2 scans when the
    # window is empty (e.g. another compound)
    watch = [s for s in spectra_rows
             if any(st.get('Mass') and 478.0 < st['Mass'] < 480.5
                    for st in steps.get(s.get('Id'), []))]
    ms2_all = [s for s in spectra_rows
               if any(st.get('Mass') and st['Mass'] > 0
                      for st in steps.get(s.get('Id'), []))]
    ms1 = next((s for s in spectra_rows
                if not any(st.get('Mass') and st['Mass'] > 0
                           for st in steps.get(s.get('Id'), []))), None)
    picks = watch[:24] if watch else ms2_all[:10]
    if ms1 is not None:
        picks = [ms1] + picks
    slog('probe: %d watch-window MS2 (prec 478.0-480.5), %d MS2 total'
         % (len(watch), len(ms2_all)))
    for s in picks:
        sid = s.get('Id')
        if sid is None:
            continue
        sec = sections.get(sid, {})
        prec = [st['Mass'] for st in steps.get(sid, [])
                if st.get('Mass') and st['Mass'] > 0]
        slog('')
        slog('probe: [scan %s] CE_Act=%r precursor=%r rt=%r'
             % (sid, psv.get(sid, {}).get(5), prec[:1], s.get('Rt')))
        moff = sec.get('maxima')
        if moff is None:
            slog('  no maxima section')
            continue
        try:
            with open(baf, 'rb') as f:
                f.seek(moff)
                bs, _t, hs = struct.unpack('<i4si', f.read(12))
                f.seek(moff + hs)
                n, po = struct.unpack('<ii', f.read(8))
                body = bs - hs - po
                slog('  maxima: bs=%d hs=%d n=%d po=%d body=%d '
                     '16n=%d 12n=%d -> %s'
                     % (bs, hs, n, po, body, 16 * n, 12 * n,
                        _maxima_layout(bs, hs, n, po)))
                if not (0 < n <= 200000 and po >= 8
                        and hs + po <= bs):
                    continue
                f.seek(moff + hs + po)
                head_hex = f.read(96).hex()
                f.seek(moff + hs + po)
                idx = struct.unpack('<%dd' % n, f.read(n * 8))
                itn = struct.unpack('<%df' % n, f.read(n * 4))
                ita = None
                if body >= 16 * n:
                    ita = struct.unpack('<%df' % n, f.read(n * 4))
                f.seek(moff + max(hs + po + 12 * n, 0))
                tail_hex = f.read(96).hex()
                f.seek(moff + bs)
                after_hex = f.read(32).hex()
                slog('  arrays head hex (96 B): %s' % head_hex)
                slog('  block tail hex (96 B): %s' % tail_hex)
                slog('  32 B after block end (what an old reader used as '
                     '"SNR"): %s' % after_hex)
                if ita is None:
                    # what an old reader read as 'SNR': 4n bytes PAST the
                    # block end -- shown here to document the bug
                    f.seek(moff + hs + po + 12 * n)
                    raw = f.read(4 * n)
                    if len(raw) == 4 * n:
                        ita_old = struct.unpack('<%df' % n, raw)
                        nlow = sum(1 for v in ita_old
                                   if v < IMPACT_SNR_FLOOR)
                        slog('  PAST-END snr: %d/%d '
                             'below floor %g'
                             % (nlow, n, IMPACT_SNR_FLOOR))
                elif len(ita):
                    nlow = sum(1 for v in ita if v < IMPACT_SNR_FLOOR)
                    slog('  real snr array: %d/%d below floor %g'
                         % (nlow, n, IMPACT_SNR_FLOOR))
                top = sorted(range(n), key=lambda j: itn[j],
                             reverse=True)[:25]
                slog('  top-25 (idx, intensity_f32, snr_or_none):')
                for j in top:
                    slog('    %.10g  %.10g  %s'
                         % (idx[j], itn[j],
                            '%.10g' % ita[j] if ita is not None
                            else 'none'))
                fin = [v for v in itn if v == v and abs(v) < 1e30]
                if fin:
                    fin.sort()
                    slog('  intensity quantiles: p05=%.6g p25=%.6g '
                         'median=%.6g p75=%.6g p95=%.6g max=%.6g'
                         % (fin[max(0, len(fin) // 20)],
                            fin[len(fin) // 4],
                            fin[len(fin) // 2],
                            fin[(3 * len(fin)) // 4],
                            fin[min(len(fin) - 1,
                                    (19 * len(fin)) // 20)],
                            fin[-1]))
        except (OSError, struct.error) as e:
            slog('  maxima dump failed: %s' % e)
        soff = sec.get('spectrum')
        if soff is not None:
            try:
                with open(baf, 'rb') as f:
                    f.seek(soff)
                    bs2, _t2, hs2 = struct.unpack('<i4si', f.read(12))
                    f.seek(soff)
                    slog('  spectrum section: bs=%d hs=%d head hex: %s'
                         % (bs2, hs2, f.read(min(160, bs2)).hex()))
            except (OSError, struct.error) as e:
                slog('  spectrum dump failed: %s' % e)


def probe_folder(d, log=print, sweep_log=None):
    """Full diagnostic dump of a .d folder: calibration-record candidates
    (baf_xtr), scan blocks (maxima + spectrum), acquisition header and
    analysis.sqlite inventory. Read-only, nothing is written."""
    log('probe (d2mz %s): folder %s' % (__version__, d))
    xtr = os.path.join(d, 'analysis.baf_xtr')
    idx = os.path.join(d, 'analysis.baf_idx')
    baf = open_baf(d)
    header_off, sections = read_index(idx)

    def _header_params_and_metadata():
        params = _read_baf_params(baf, header_off)
        log('probe: acquisition parameters: %d row(s)' % len(params))
        for ri in sorted(params):
            nm = params[ri][1][0]
            if nm:
                log('  %d %s' % (ri, nm.decode('utf-8', 'replace')))
        # sample scans SPREAD ACROSS THE WHOLE RUN (the first
        # scans are always MS1 -- the MS/MS rows, when present, only
        # appear on later scans). First sampled scan: full dump; the
        # others: only the rows that CHANGE (the constant rows are
        # instrument status shared by every scan).
        all_md = sorted(f for f, o in sections.items() if 'metadata' in o)
        n_md = len(all_md)
        if n_md:
            idxs = sorted(set([0, n_md // 8, n_md // 4, n_md // 2,
                               (3 * n_md) // 4, n_md - 2, n_md - 1]))
            md_flags = [all_md[i] for i in idxs]
        else:
            md_flags = []
        known = {}
        for k, flag in enumerate(md_flags):
            log('probe: metadata of scan %s [%d/%d] (values of the '
                'parameters above):%s'
                % (flag, k + 1, len(md_flags),
                   '' if k == 0 else ' (only changes shown)'))
            got = _read_metadata_values(baf, sections[flag]['metadata'],
                                        params, log,
                                        known if k else None)
            if k == 0:
                known.update(got)

    lock_scans = []
    try:
        for p, fo in calibration_candidates(xtr):
            if (len(p) >= fo + 8
                    and struct.unpack_from('<i', p, fo)[0] == 10):
                v = struct.unpack_from('<f', p, fo + 4)[0]
                if v != -1.0:
                    lock_scans.append(struct.unpack_from('<i', p, 32)[0])
    except Exception:
        pass
    lock_scans = sorted(set(lock_scans))[:3]
    if lock_scans:
        log('probe: lock-mass scans detected: %s' % lock_scans)

    for name, fn in (
            ('xtr candidates', lambda: probe_calibration(xtr, log,
                                                          full=False)),
            ('xtr overview', lambda: _probe_xtr_overview(xtr, log)),
            ('scan blocks', lambda: _probe_blocks(
                baf, sections, log, extra_scans=tuple(lock_scans))),
            ('idx flags', lambda: _probe_idx_flags(idx, baf, log)),
            ('acquisition header', lambda: _probe_acquisition_header(
                baf, header_off, log)),
            ('header params + metadata', _header_params_and_metadata),
            ('sqlite', lambda: _probe_sqlite(d, log)),
            ('v7 model probe', lambda: _probe_v7(d, baf, sections, log)),
            ('v8 full-peak dump', lambda: _probe_v8(d, baf, sections,
                                                    log)),
            ('sweep-MS2 dump', lambda: _probe_sweep_ms2(
                d, baf, header_off, sections,
                log=sweep_log if sweep_log is not None else log,
                slog=sweep_log))):
        try:
            fn()
        except Exception as e:  # probe must survive any single failure
            log('probe: %s failed: %s' % (name, e))


def read_calibration(xtr_path):
    matches = calibration_candidates(xtr_path)
    if not matches:
        raise ValueError('calibration record not found in %s' % xtr_path)
    good = [(p, fo) for p, fo in matches
            if struct.unpack_from('<i', p, fo + 8)[0] in (4, 5)]
    if not good:
        seen = sorted({struct.unpack_from('<i', p, fo + 8)[0]
                       for p, fo in matches})
        # impact-generation qTOF files carry only mode-10 transformator
        # records: no frequency calibration. The quadratic model is
        # fitted per file from the calibrant (fit_impact_calibration).
        return dict(mode=seen[0] if seen else 0, modes=seen)
    p, fo = good[0]
    mode = struct.unpack_from('<i', p, fo + 8)[0]
    alpha = struct.unpack_from('<d', p, fo + 12)[0] if mode == 5 else 0.0
    beta = struct.unpack_from('<d', p, fo + 20)[0]
    if mode == 4:
        beta = -beta
    return dict(
        mode=mode, alpha=alpha, beta=beta,
        freq_high=struct.unpack_from('<d', p, fo + 28)[0],
        size=struct.unpack_from('<i', p, fo + 40)[0],
        freq_low=struct.unpack_from('<d', p, fo + 48)[0],
        freq_width=struct.unpack_from('<d', p, fo + 56)[0])


def index_to_mz(cal, index):
    if cal.get('quad') is not None:
        a, b, c = cal['quad']
        return a * index * index + b * index + c
    den = cal['freq_width'] * (cal['size'] - index) / cal['size'] - cal['beta']
    return cal['freq_high'] / den + cal['alpha'] / cal['freq_high']


def quad_mz_to_index(quad, mz):
    """Inverse of the impact quadratic  m = a*i^2 + b*i + c  (i >= 0)."""
    a, b, c = quad
    disc = b * b - 4.0 * a * (c - mz)
    if disc < 0.0:
        return None
    return 0.5 * (-b + math.sqrt(disc)) / a


def _solve3(rows, rhs):
    """Solve a 3x3 linear system (Gaussian elimination, pure stdlib)."""
    m = [list(r) + [v] for r, v in zip(rows, rhs)]
    for col in range(3):
        piv = max(range(col, 3), key=lambda r: abs(m[r][col]))
        if abs(m[piv][col]) < 1e-300:
            return None
        m[col], m[piv] = m[piv], m[col]
        pv = m[col][col]
        m[col] = [v / pv for v in m[col]]
        for r in range(3):
            if r != col and m[r][col]:
                f = m[r][col]
                m[r] = [v - f * w for v, w in zip(m[r], m[col])]
    return [m[r][3] for r in range(3)]


def _quad_lstsq(points):
    """Weighted quadratic least squares on [(x, y, w), ...] -> (a, b, c).
    y = a*x^2 + b*x + c, fitted on the CENTERED variable u = x - xbar.

    CENTERING IS MANDATORY: raw TOF
    indices reach ~350000, so the raw normal-equation matrix holds
    S4 ~ 1e23 next to S0 ~ 1e2 -> condition ~1e21 and the 3x3 solve
    loses ALL double precision (the catastrophic 45-295 ppm junk fits
    and even 1e10-ppm divergences reproduced on real calibrant
    anchors). On the centered variable the matrix is well-conditioned;
    the result is converted back to raw (a, b, c) exactly."""
    sw = sum(w for _x, _y, w in points)
    if sw <= 0 or not points:
        return None
    xbar = sum(w * x for x, _y, w in points) / sw
    su0 = su1 = su2 = su3 = su4 = sy = suy = su2y = 0.0
    for x, y, w in points:
        u = x - xbar
        u2 = u * u
        su0 += w
        su1 += w * u
        su2 += w * u2
        su3 += w * u * u2
        su4 += w * u2 * u2
        sy += w * y
        suy += w * u * y
        su2y += w * u2 * y
    sol = _solve3([(su4, su3, su2), (su3, su2, su1), (su2, su1, su0)],
                  [su2y, suy, sy])
    if not sol:
        return None
    A, B, C = sol
    a = A
    b = B - 2.0 * A * xbar
    c = C + A * xbar * xbar - B * xbar
    return (a, b, c)


def _robust_quad_fit(points, rounds=6):
    """Robust (MAD-rejection) weighted quadratic fit.
    points: [(idx, mz, weight)]. Returns (a, b, c, med_ppm, max_ppm, n)."""
    pts = list(points)
    fit = _quad_lstsq(pts)
    if fit is None:
        return None
    for _ in range(rounds):
        a, b, c = fit
        res = sorted((a * x * x + b * x + c - y) / y * 1e6
                     for x, y, _w in pts)
        if not res:
            return None
        mad = res[len(res) // 2]
        mad = abs(mad) if mad else 1e-9
        keep = [(x, y, w) for x, y, w in pts
                if abs((a * x * x + b * x + c - y) / y * 1e6)
                < max(3.0 * mad, 2.0)]
        if len(keep) < 6 or len(keep) == len(pts):
            break
        pts = keep
        fit = _quad_lstsq(pts)
        if fit is None:
            return None
    a, b, c = fit
    res = [abs((a * x * x + b * x + c - y) / y * 1e6) for x, y, _w in pts]
    res.sort()
    return (a, b, c, res[len(res) // 2], res[-1], len(pts))


class CalibrantAbsentError(ValueError):
    """Impact file with no usable calibrant segment (no Segment==2 scan,
    or no calibrant ladder matched there). The file is NOT
    ground through the whole run looking for calibrant ions -- the file is
    DEFERRED and its calibration is borrowed at the end from a
    same-mode (same polarity) file of the same folder."""


def _impact_scan_list(baf, sections, spectra_rows, steps, min_int=3.0,
                      cal_only=False):
    """MS1 scans as raw [(rt, sorted-idx-array, (idx, SNR) pairs, is_cal)]
    records for the calibration fit (single read). is_cal marks the
    CALIBRANT SEGMENT scans (Spectra.Segment == 2 -- the
    calibrant sits in segment 2, not segment 1; validated on the lab's own runs: segment-2 MS1 scans form a short
    burst at the start of the run, RT ~6-18 s).
    cal_only=True reads ONLY the calibrant-segment maxima blocks
    (~10-60 blocks instead of ~1000-3000: the whole point of knowing
    where the calibrant lives).
    min_int is a signal-to-noise threshold (noise peaks sit at SNR=1)."""
    scans = []
    n_sp = len(spectra_rows)
    for k, s in enumerate(spectra_rows, 1):
        _progress(k, n_sp, 'calibrating scans', step=100)
        sid, rt = s.get('Id'), s.get('Rt', s.get('RT'))
        if sid is None or rt is None:
            continue
        if any(st.get('Mass') and st['Mass'] > 0 for st in steps.get(sid, [])):
            continue  # MS2
        is_cal = s.get('Segment') == 2
        if cal_only and not is_cal:
            continue
        moff = sections.get(sid, {}).get('maxima')
        if moff is None:
            continue
        m = _read_maxima_raw(baf, moff, f32_ints=True)
        if m and m['n']:
            pairs = [(i, v) for i, v in zip(m['idx'], m['itn'])
                     if v >= min_int]
            if pairs:
                pairs.sort()
                # (rt, sorted idx array, pairs, is_cal) -- the sorted
                # array enables bisect in _assign_targets, and keeping
                # the raw records lets fit_lock_model observe the lock
                # masses on the SAME data without re-reading the blocks
                scans.append((float(rt),
                              [i for i, _v in pairs], pairs, is_cal))
    _progress_end()
    return scans


def _assign_targets(scans, targets, quad, win):
    """For each calibrant target: per scan keep the best peak scored
    SNR * exp(-(d/(win/2.5))^2) within +/-win of the predicted index
    (intensity-weighted proximity, as validated at 0.25 ppm on the
    Saule reference extraction), then SNR-weighted mean index across
    scans. Targets seen in fewer than IMPACT_TARGET_MIN_FRAC of the
    scans are dropped. Returns anchor points [(idx, mz, w)]."""
    nscan = len(scans)
    min_seen = max(2, int(IMPACT_TARGET_MIN_FRAC * nscan))
    sigma = win / 2.5
    anchors = []
    for mz in targets:
        it = quad_mz_to_index(quad, mz)
        if it is None or not (IMPACT_IDX_MIN < it < IMPACT_IDX_MAX):
            continue
        wsum = 0.0
        xsum = 0.0
        seen = 0
        for xs, pairs in scans:
            best = None
            # bisect: peaks are sorted by index, so only the slice
            # [it-win, it+win] is scored (O(log n) instead of O(n))
            lo = bisect.bisect_left(xs, it - win)
            hi = bisect.bisect_right(xs, it + win)
            for i, v in pairs[lo:hi]:
                d = abs(i - it)
                if d <= win:
                    sc = v * math.exp(-(d / sigma) ** 2)
                    if best is None or sc > best[0]:
                        best = (sc, i, v)
            if best is not None:
                seen += 1
                w = best[2]      # raw intensity: cluster-to-junk contrast
                wsum += w        # is what wins the assignment (no cap)
                xsum += best[1] * w
        if seen >= min_seen and wsum > 0:
            anchors.append((xsum / wsum, mz, wsum / seen))
    return anchors


def _fit_quad_series(scans, targets, boot, log=None):
    """Bootstrap -> shrinking assignment windows -> robust quadratic.
    Returns (quad, stats dict) or None."""
    quad = boot
    stats = None
    for win in IMPACT_FIT_WINDOWS:
        if log is not None:
            log('      window %.0f: assigning %d targets on %d scans...'
                % (win, len(targets), len(scans)))
        anchors = _assign_targets(scans, targets, quad, win)
        if len(anchors) < 4:
            return stats
        fit = _robust_quad_fit(anchors)
        if fit is None:
            return stats
        quad = fit[:3]
        stats = dict(quad=quad, n_targets=fit[5], med_ppm=fit[3],
                     max_ppm=fit[4], anchors=anchors)
    return stats


def _cluster_families(ions):
    """Group calibrant cluster targets into (AcO)-ladder families.
    A Na-formate run only transmits the AcO=0 ladder; matching absent
    families produces junk anchors, so each family is fitted on its
    own and the best one wins."""
    fams = {}
    for mz, formula in ions:
        m = re.search(r'\(AcO\)(\d+)', formula)
        key = int(m.group(1)) if m else 0
        fams.setdefault(key, []).append(mz)
    return [(k, fams[k]) for k in sorted(fams)]


IMPACT_FIT_MAX_PPM = 2.0     # sanity gate on the per-family fit: a REAL
                             # calibrant segment fit lands at 0.08-0.73 ppm
                             # median (validated on real GA4
                             # POS, RA NEG and Parabens runs); 2.0 keeps a
                             # ~3x margin. The former 5.0 ppm gate let a
                             # partially-coincidental family print as
                             # 'fits' at 4.22 ppm median / 10.28 max on a
                             # real NEG file (the
                             # best-family selection still won, but the
                             # line lied)
IMPACT_FIT_MIN_W = 1.0       # base-peak anchor weight floor (INTENSITY
                             # semantics: real clusters
                             # sit at 1e4..4e6; the old SNR-scaled floors
                             # 7/12 were tuned for the wrong interleaved
                             # parse and would reject the Taraxacum-scale
                             # intensities 20..272)
IMPACT_SNR_CAP = None        # NO winsorization any more: the peak weight
                             # is the RAW INTENSITY and the cluster-to-junk
                             # contrast (~200x+) is exactly what makes the
                             # assignment robust; capping destroyed it
                             # (Scoring = intensity * proximity, the
                             # validated reference-script recipe)
IMPACT_SNR_FLOOR = 2.0       # per-peak SNR floor (read_maxima): the
                             # separate SNR array marks the digitizer-noise
                             # bucket at ~1.0 (thousands per scan)

# EXACT calibrant lists, ANCHORED in the script as a fallback. They
# are only used when the .m acquisition
# method carries NO Calibration_LastCalibrationReferenceMass vector.
# SOURCE: the 38-mass vector read in the lab's OWN .m methods (both the
# Parabens POS and NEG acquisitions carry the same 38 cation clusters:
# Na(HCOONa)n 1-6, Na(CH3COONa)n 1-4, Na(NaAc)x(NaF)y up to 948.872683
# + Na(NaAc)4(NaF)1). The .m vector is POLARITY-INDEPENDENT (cation
# formulas even in negative methods): the negative run transmits the
# ANION counterparts, anion = cation - 2 Na + 2 e (validated on a
# real DataAnalysis capture: cation 158.964069 -> anion
# 112.98563, spacing of the transmitted [Na_n(HCOO)_{n+1}]- ladder).
NA_MONO = 22.98976928      # monoisotopic sodium
E_MONO = 0.00054857990     # electron mass
HCOO_MONO = 44.99765427    # formate anion, monoisotopic
ACO_MONO = 59.01330434     # acetate anion, monoisotopic


def _anion_counterparts(masses, min_mz=50.0):
    """Cation reference masses -> the anion counterparts transmitted in
    negative polarity (same cluster series minus 2 Na + 2 e)."""
    d = 2.0 * NA_MONO - 2.0 * E_MONO
    return [m - d for m in masses if m - d >= min_mz]


def _cation_counterparts(masses, min_mz=50.0):
    """Anion reference masses -> the cation counterparts (same series
plus 2 Na - 2 e). Used when a calibration table of the WRONG polarity
was referenced in the acquisition sequence (a
post-processing file mis-referenced in the sequence -- the observed
negative .m carries the cation table, and the same mistake can hit a
positive run carrying the anion table)."""
    d = 2.0 * NA_MONO - 2.0 * E_MONO
    return [m + d for m in masses if m + d >= min_mz]


def _matches_series(masses, series, tol=0.01):
    """True when >= 75% of the reference vector masses match a member
of the given built-in cluster series (mis-referenced-table detector)."""
    if len(masses) < 4:
        return False
    hits = sum(1 for m in masses
               if any(abs(m - s) <= tol for s in series))
    return hits >= max(4, int(round(0.75 * len(masses))))


def _theoretical_cluster_targets(min_mz=50.0, max_mz=1100.0):
    """ALL theoretical calibrant clusters, detected or not, listed by
    Bruker or not (integrating the undetected/unlisted
    clusters consolidates the calibration at high masses). Enumeration of the sodium formate/acetate cluster
    compositions, verified to reproduce (a) the 38 masses of the .m
    vector EXACTLY, (b) the observed continuation of the
    calibrant-segment captures (444.9294 -> 906.8258):
    positive [Na_{b+c+1}(HCOO)_b(AcO)_c]+ = (b+c+1) Na + b HCOO + c AcO
    - e; negative counterpart [Na_{b+c-1}(HCOO)_b(AcO)_c]- = the cation
    - 2 Na + 2 e (the validated transform). Returns
    {'positive': [...], 'negative': [...]} sorted ascending."""
    out = {'positive': [], 'negative': []}
    seen = {'positive': set(), 'negative': set()}
    n = int(max_mz / min(HCOO_MONO, NA_MONO)) + 2
    for b in range(0, n):
        for c in range(0, n):
            if b + c < 1:
                continue
            core = b * HCOO_MONO + c * ACO_MONO
            cat = core + (b + c + 1) * NA_MONO - E_MONO
            ani = core + (b + c - 1) * NA_MONO + E_MONO
            for pol, m in (('positive', cat), ('negative', ani)):
                if min_mz <= m <= max_mz:
                    k = round(m, 6)
                    if k not in seen[pol]:
                        seen[pol].add(k)
                        out[pol].append(m)
    out['positive'].sort()
    out['negative'].sort()
    return out


IMPACT_EXACT_TARGETS = {
    # positive: the TargetScreener_NaFormateAcetate_POS.ref v2 list
    # (validated against observed clusters): the 38-mass
    # .m vector + the 15 observed continuation clusters
    # (498.901189 -> 906.825727)
    'positive': [
        90.976645, 104.992295, 158.964069, 172.979719, 186.995369, 226.951493,
        240.967143, 254.982793, 268.998443, 294.938917, 308.954567, 322.970217,
        336.985867, 351.001517, 362.926341, 376.941991, 390.957641, 404.973291,
        418.988941, 430.913765, 444.929415, 458.945065, 472.960715, 486.976365,
        498.901189, 512.916839, 526.932489, 540.948139, 554.963789, 566.888613,
        568.979439, 580.904263, 594.919913, 608.935563, 634.876033, 648.891683,
        662.907337, 676.922987, 702.863456, 716.879106, 730.894761, 744.910411,
        770.850880, 784.866530, 798.882185, 812.897835, 838.838303, 852.853953,
        866.869609, 880.885259, 906.825727, 934.857033, 948.872683
    ],
    # negative: the TargetScreener_NaFormateAcetate_NEG.ref v2 list
    # (validated against observed clusters): the 2019
    # actives + HCOO(NaCOOH)5-12 (384.935 -> 860.847) and
    # Ac(NaF)8-11 (602.913 -> 806.876)
    'negative': [
        112.985627, 127.001277, 141.016927, 180.973051, 194.988701, 209.004351,
        223.020001, 248.960475, 262.976125, 276.991775, 291.007425, 305.023075,
        316.947899, 330.963549, 344.979199, 358.994849, 373.010499, 384.935323,
        398.950973, 412.966623, 426.982273, 440.997923, 452.922747, 466.938397,
        480.954047, 494.969697, 508.985347, 520.910170, 534.925821, 548.941471,
        562.957121, 576.972771, 588.897594, 602.913244, 616.928895, 630.944545,
        644.960195, 656.885018, 670.900668, 684.916319, 698.931969, 712.947619,
        724.872442, 738.888092, 752.903742, 792.859866, 806.875516, 820.891166,
        860.847290, 888.878590, 956.866014
    ],
}


def fit_impact_calibration(baf, sections, spectra_rows, steps, opts,
                           method_pol=None, log=print, method_targets=None):
    """Impact-generation qTOF: per-file quadratic index -> m/z calibration
    fitted on the CALIBRANT SEGMENT (Segment==2) ONLY, on the EXACT
    calibrant list -- the reference-mass vector declared by the .m
    acquisition method (Calibration_LastCalibrationReferenceMass) when
    it has one, otherwise the declared polarity's Na-formate/acetic
    cluster list. There is NO other fit path:
    1. segment 2, 2. another file -- when the segment or
    the list does not produce a usable fit the file raises
    CalibrantAbsentError and is DEFERRED, its calibration being borrowed
    at the end from a same-mode file of the same folder (all files of
    the project share the same structure; a quad from a same-polarity
    neighbour is the Bruker practice). The adduct-family scan, the
    base-cluster anchor and the legacy all-MS1 fit were REMOVED.
    Returns a cal dict with quad=(a, b, c) and segment_fit=True."""
    n_ms1 = sum(1 for r in spectra_rows
                if not any(st.get('Mass') and st['Mass'] > 0
                           for st in steps.get(r.get('Id'), [])))
    recs = _impact_scan_list(baf, sections, spectra_rows, steps,
                             cal_only=True)
    if not recs:
        raise CalibrantAbsentError(
            'no calibrant segment (Segment==2) scan in this run')
    segment_fit = True
    cal_scans = [(xs, pairs) for _rt, xs, pairs, _c in recs]
    log('    calibrant segment (Segment==2): %d scans (of %d MS1)'
        % (len(cal_scans), n_ms1))
    # TARGET LIST: the EXACT .m reference vector when the method
    # declares one; otherwise the declared polarity's cluster list.
    # ONE single fit -- no adduct-family scan.
    if method_targets and len(method_targets) >= 4:
        log('    calibrant reference list read from .m method: %d '
            'masses -> direct fit' % len(method_targets))
        # The sequence can reference a calibration table of the WRONG
        # polarity in the .m method (a post-
        # processing file mis-referenced in the acquisition sequence).
        # Detect the FORM of the vector against the built-in cluster
        # series and apply the counterpart transform when the declared
        # polarity does not match the vector's form.
        vec_is_anionic = _matches_series(
            method_targets, IMPACT_EXACT_TARGETS['negative'])
        vec_is_cationic = _matches_series(
            method_targets, IMPACT_EXACT_TARGETS['positive'])
        if method_pol == 'negative':
            if vec_is_anionic:
                targets = list(method_targets)
                log('    negative run: reference vector is already '
                    'anionic -> fitting it directly')
            else:
                # cation table referenced in a negative method: the run
                # transmits the ANION counterparts (same series
                # - 2 Na + 2 e, validated on a real capture:
                # 158.964069 -> 112.98563)
                targets = _anion_counterparts(method_targets)
                log('    negative run: fitting the anion counterparts '
                    '(%d masses >= 50 Da)' % len(targets))
        elif method_pol == 'positive' and vec_is_anionic \
                and not vec_is_cationic:
            # anion table mis-referenced in a positive method: fit the
            # CATION counterparts (same series + 2 Na - 2 e)
            targets = _cation_counterparts(method_targets)
            log('    positive run: reference vector is anionic '
                '(mis-referenced calibration table) -> fitting the '
                'cation counterparts (%d masses >= 50 Da)'
                % len(targets))
        else:
            targets = list(method_targets)
        boot_pol = method_pol or 'positive'
    elif method_pol in ('positive', 'negative'):
        targets = list(IMPACT_EXACT_TARGETS[method_pol])
        log('    no reference vector in .m method: using the built-in '
            'exact %s list (%d masses)'
            % (method_pol, len(targets)))
        boot_pol = method_pol
    else:
        raise CalibrantAbsentError(
            'no calibrant reference vector in .m method and no '
            'declared polarity')
    # Integrate the clusters NOT detected
    # (or not even listed by Bruker) into the target list, to
    # consolidate the calibration at high masses. The extra targets
    # simply stay unmatched when the peak is absent; when the high-mass
    # clusters ARE transmitted (they are, on the lab's calibrant
    # segments), they become additional quadratic anchors.
    _th = _theoretical_cluster_targets()
    _add = [m for m in _th[boot_pol]
            if all(abs(m - t) > 0.002 for t in targets)]
    if _add:
        targets = sorted(targets + _add)
        log('    extended with %d theoretical clusters not detected/'
            'unlisted by Bruker (up to 1100 Da)' % len(_add))
    st = _fit_quad_series(cal_scans, targets,
                          IMPACT_QUAD_BOOTSTRAP[boot_pol], None)
    reason = None
    if not st:
        reason = 'no match'
    elif st['n_targets'] < 4:
        reason = 'too few targets (%d)' % st['n_targets']
    elif st['med_ppm'] > IMPACT_FIT_MAX_PPM:
        reason = 'gate failed (med %.1f ppm)' % st['med_ppm']
    elif st['anchors'][0][2] < IMPACT_FIT_MIN_W:
        reason = 'base peak too weak (matched noise)'
    if reason is not None:
        # Diagnostic: dump the strongest calibrant-segment maxima with
        # their coarse bootstrap m/z so the user can SEE what
        # segment 2 actually contains instead of a silent junk fit
        a0, b0, c0 = IMPACT_QUAD_BOOTSTRAP['positive']
        seen = {}
        # ALL peaks of the first scans (a pairs[:15] slice would only
        # ever show the LOWEST-index peaks: pairs are index-sorted)
        for _xs, pairs in cal_scans[:3]:
            for i, v in pairs:
                seen[i] = max(seen.get(i, 0.0), v)
        tops = sorted(seen.items(), key=lambda t: -t[1])[:12]
        if tops:
            log('    calibrant fit rejected (%s); strongest segment-2 '
                'maxima (coarse bootstrap m/z):' % reason)
            for i, v in tops:
                log('      m/z ~ %.4f  (idx %.0f, SNR %.1f)'
                    % (a0 * i * i + b0 * i + c0, i, v))
        raise CalibrantAbsentError('calibrant fit rejected (%s)' % reason)
    best = st
    best['polarity'] = method_pol or boot_pol
    log('    impact quadratic calibration (%s series): %d targets, '
        'residual %.2f ppm median / %.2f max'
        % (best['polarity'], best['n_targets'], best['med_ppm'],
           best['max_ppm']))
    cal = dict(mode=10, quad=best['quad'], polarity=best['polarity'],
                fit=dict(n_targets=best['n_targets'],
                         med_ppm=best['med_ppm'], max_ppm=best['max_ppm']),
                segment_fit=segment_fit,
                _impact_scans=recs)
    # _impact_scans: raw MS1 records, freed right after the pre-pass
    # (segment fits no longer observe lock masses over the run)
    return cal



def read_maxima(baf, off, cal):
    """Bruker peak-picked (maxima) block -> [(mz, intensity), ...]

    Impact generation (layout re-probed 2026-11-07):
    [n idx f64][n INTENSITY f32][n SNR f32] consecutive arrays (see
    _read_maxima_raw). The true per-peak INTENSITY is returned; the
    separate SNR array filters the digitizer-noise bucket (SNR ~1)."""
    impact = cal.get('quad') is not None or cal.get('mode') not in (4, 5)
    with open(baf, 'rb') as f:
        f.seek(off)
        bs, _, hs = struct.unpack('<i4si', f.read(12))
        if bs <= 0 or hs < 12 or hs > bs:
            return []
        f.seek(off + hs)
        head = f.read(8)
        n, po = struct.unpack('<ii', head)
        if n < 0 or n > 200000 or po < 8 or hs + po > bs:
            return []
        lay = _maxima_layout(bs, hs, n, po)
        if lay is None:
            return []      # block smaller than the 12n layout
        f.seek(off + hs + po)
        idx = struct.unpack('<%dd' % n, f.read(n * 8)) if n else ()
        if impact:
            itn = struct.unpack('<%df' % n, f.read(n * 4)) if n else ()
            ita = (struct.unpack('<%df' % n, f.read(n * 4))
                   if (n and lay == 'A3') else ())
        elif lay == 'A2':
            itn = struct.unpack('<%df' % n, f.read(n * 4)) if n else ()
            ita = ()
        else:
            itn = struct.unpack('<%dd' % n, f.read(n * 8)) if n else ()
            ita = ()
    out = []
    for k in range(len(itn)):
        if impact and ita and ita[k] < IMPACT_SNR_FLOOR:
            continue   # digitizer-noise bucket (SNR ~ 1.0); A2 blocks
                       # have ita == () -> NO filtering (the
                       # code read garbage past the MS2 block end here
                       # and dropped intense fragments at random)
        try:
            out.append((index_to_mz(cal, idx[k]), itn[k]))
        except (ZeroDivisionError, ValueError):
            continue
    return out


def _dad_open(db):
    """Open a chromatography-data(-pre).sqlite and locate the DAD cube.
    Returns (conn, wl_axis, time_col, intens_col) or None when the
    file does not carry the expected Spectra/SpectrumSources tables."""
    try:
        c = sqlite3.connect(db)
    except sqlite3.Error:
        return None
    try:
        tables = [r[0] for r in c.execute(
            "select name from sqlite_master where type='table'")]
        if 'Spectra' not in tables or 'SpectrumSources' not in tables:
            c.close()
            return None
        wl = None
        for row in c.execute('select * from SpectrumSources'):
            for val in row:
                if isinstance(val, bytes) and len(val) >= 16:
                    cand = list(struct.unpack(
                        '<%dd' % (len(val) // 8), val[:8 * (len(val) // 8)]))
                    if cand and 150.0 < cand[0] < 250.0 and \
                            cand[-1] > cand[0] and len(cand) >= 100:
                        wl = cand       # wavelength axis (nm, ascending)
                        break
            if wl:
                break
        if not wl:
            c.close()
            return None
        cols = [r[1] for r in c.execute('pragma table_info("Spectra")')]
        tcol = next((k for k in cols if 'time' in k.lower()), None)
        icol = next((k for k in cols if 'intens' in k.lower()), None)
        if not tcol or not icol:
            c.close()
            return None
        return c, wl, tcol, icol
    except (sqlite3.Error, struct.error):
        try:
            c.close()
        except Exception:
            pass
        return None


def _read_dad_cube(db, wanted):
    """DAD traces at the requested wavelengths ONLY. The cube holds
    ~48000 time rows x 611 wavelengths (~29M float values) per file;
    materialising it whole is what made the conversion crawl. Here
    each row unpacks just len(wanted) floats (one struct.unpack_from
    per wavelength), and the full spectrum is re-read on demand for
    the few apex scans only (see _read_dad_full_at).
    Returns [(rt_sec, {wl_nm: intensity})]."""
    op = _dad_open(db)
    if op is None:
        return []
    c, wl, tcol, icol = op
    n = len(wl)
    picks = {}
    for w in wanted:
        try:
            w = float(w)
        except (TypeError, ValueError):
            continue
        i = min(bisect.bisect_left(wl, w), n - 1)
        if i > 0 and abs(wl[i - 1] - w) <= abs(wl[i] - w):
            i -= 1
        picks[w] = i
    if not picks:
        c.close()
        return []
    out = []
    try:
        q = 'select "%s", "%s" from "Spectra"' % (tcol, icol)
        for t, blob in c.execute(q):
            if t is None or not isinstance(blob, bytes) or len(blob) < 4 * n:
                continue
            out.append((float(t),
                        {w: struct.unpack_from('<f', blob, 4 * i)[0]
                         for w, i in picks.items()}))
    except sqlite3.Error:
        c.close()
        return []
    c.close()
    out.sort(key=lambda p: p[0])
    return out


def _read_dad_full_at(db, rts):
    """Full wavelength spectra of the cube rows nearest each requested
    rt. Reads the time column once (cheap), then unpacks the intensity
    blob only for the few nearest rows."""
    if not rts:
        return []
    op = _dad_open(db)
    if op is None:
        return []
    c, wl, tcol, icol = op
    n = len(wl)
    out = []
    try:
        rows = c.execute('select rowid, "%s" from "Spectra"' % tcol).fetchall()
        rows = [(rid, float(t)) for rid, t in rows if t is not None]
        if not rows:
            c.close()
            return []
        for rt_a in rts:
            rid = min(rows, key=lambda p: abs(p[1] - rt_a))[0]
            r = c.execute('select "%s", "%s" from "Spectra" '
                          'where rowid = ?' % (tcol, icol),
                          (rid,)).fetchone()
            if r is None:
                continue
            t, blob = r
            if t is None or not isinstance(blob, bytes) or len(blob) < 4 * n:
                continue
            vals = struct.unpack('<%df' % n, blob[:4 * n])
            out.append((float(t), dict(zip(wl, vals))))
    except (sqlite3.Error, struct.error):
        pass
    finally:
        try:
            c.close()
        except Exception:
            pass
    out.sort(key=lambda p: p[0])
    return out


def _dbl_blob(blob):
    """Little-endian float64 array from a sqlite blob ([] when not a
    blob or empty)."""
    if not isinstance(blob, bytes) or len(blob) < 8:
        return []
    return list(struct.unpack('<%dd' % (len(blob) // 8), blob))


def _flt_blob(blob):
    """Little-endian float32 array from a sqlite blob ([] when not a
    blob or empty). The real TraceChunks store Times as f64 and
    Intensities as f32."""
    if not isinstance(blob, bytes) or len(blob) < 4:
        return []
    return list(struct.unpack('<%df' % (len(blob) // 4), blob))


def read_all_trace_chunks(d_folder):
    """ALL named traces of chromatography-data(-pre).sqlite
    (TraceSources/TraceChunks: the format of the lab's recent
    otoControl runs): UV_VIS_1, UV_VIS_2, Temp_Lamphouse,
    ColumnOven_Temp, Pump_Pressure, BPC, ... Every existing trace goes
    into the mzML (cheap in space). Returns [(description, instrument, [(rt, intensity)])]
    sorted by trace Id ([] when the file carries no named trace)."""
    for cand in ('chromatography-data.sqlite',
                 'chromatography-data-pre.sqlite'):
        p = os.path.join(d_folder, cand)
        if not os.path.exists(p):
            continue
        try:
            c = sqlite3.connect(p)
            tables = [r[0] for r in c.execute(
                "select name from sqlite_master where type='table'")]
            if 'TraceSources' not in tables or 'TraceChunks' not in tables:
                c.close()
                continue
            srcs = c.execute(
                'select Id, Description, Instrument from TraceSources '
                'order by Id').fetchall()
            out = []
            for tid, desc, inst in srcs:
                pts = []
                for tb, ib in c.execute(
                        'select Times, Intensities from TraceChunks '
                        'where Trace = ? order by rowid', (tid,)):
                    ts = _dbl_blob(tb)
                    its = _flt_blob(ib)
                    if len(ts) != len(its):
                        its = (its + [0.0] * len(ts))[:len(ts)]
                    pts.extend(zip(ts, its))
                pts = [(t, v) for t, v in pts if v is not None]
                pts.sort(key=lambda q: q[0])
                if pts:
                    out.append((str(desc), str(inst), pts))
            c.close()
            if out:
                return out
        except sqlite3.Error:
            continue
    return []


def hystar_uv_wavelengths(d_folder):
    """Wavelength of each acquired UV channel from the Hystar LC method
    (the .m folder of the .d: hystar.method, e.g. UV_VIS_1.Wavelength =
    280 nm). Returns {channel_number: wavelength_nm} ({} when absent --
    the method file is double-XML-escaped)."""
    mdir = find_method_folder(d_folder)
    if not mdir:
        return {}
    out = {}
    for fn in sorted(os.listdir(mdir)):
        if 'hystar' not in fn.lower():
            continue
        try:
            txt = open(os.path.join(mdir, fn), encoding='utf-8',
                       errors='replace').read()
        except OSError:
            continue
        txt = txt.replace('&amp;', '&').replace('&lt;', '<') \
                 .replace('&gt;', '>')   # double-escaped XML
        for m in re.finditer(
                r'<Name>UV_VIS_(\d+)\.Wavelength</Name>.*?'
                r'<Value>([\d.]+)</Value>', txt, re.S):
            try:
                out[int(m.group(1))] = float(m.group(2))
            except ValueError:
                pass
        if out:
            break
    return out


def hystar_lc_method(d_folder):
    """LC program (gradient, flow, solvents, oven setpoint) from the
    HyStar method (hystar.method, in the .m folder of the .d).
    the gradient is NOT recorded as a data channel by
    HyStar (the chromatography database only holds UV_VIS_1/2,
    Temp_Lamphouse, ColumnOven_Temp, Pump_Pressure, BPC) -- but the
    PROGRAM is in the method: 'PumpModule Time Table' rows give
    time -> Pump.%X.Value / Pump.Flow.Nominal / Pump.Curve.
    Returns None when no pump gradient table is found. The method
    file is double-XML-escaped (two unescape rounds needed)."""
    mdir = find_method_folder(d_folder)
    if not mdir:
        return None
    for fn in sorted(os.listdir(mdir)):
        if 'hystar' not in fn.lower():
            continue
        try:
            txt = open(os.path.join(mdir, fn), encoding='utf-8',
                       errors='replace').read()
        except OSError:
            continue
        for _ in (1, 2):        # double-escaped XML
            txt = txt.replace('&amp;', '&').replace('&lt;', '<') \
                     .replace('&gt;', '>')
        flow, grad, curve = [], {}, []
        for name, body in re.findall(
                r'<ListValues ListName="([^"]*)">(.*?)</ListValues>',
                txt, re.S):
            if name != 'PumpModule Time Table':
                continue
            v1 = re.search(r'<Value1>([^<]*)</Value1>', body)
            v2 = re.search(r'<Value2>([^<]*)</Value2>', body)
            v3 = re.search(r'<Value3>([^<]*)</Value3>', body)
            if not (v1 and v2 and v3):
                continue
            try:
                t = float(v1.group(1))
            except ValueError:
                continue
            func, val = v2.group(1).strip(), v3.group(1).strip()
            num = re.match(r'([-.\d]+)', val)
            m = re.match(r'Pump\.([%][A-D])\.Value$', func)
            if m and num:
                grad.setdefault(m.group(1), []).append(
                    (t, float(num.group(1))))
            elif func.endswith('Flow.Nominal') and num:
                flow.append((t, float(num.group(1))))
            elif func.endswith('Curve') and num:
                curve.append(float(num.group(1)))
        if not grad:
            continue
        solvents = {}
        for comp in ('%A', '%B', '%C', '%D'):
            m = re.search(r'<Struct Name="%s"[^>]*>(.*?)</Struct>'
                          % re.escape(comp), txt, re.S)
            if m:
                eq = re.search(r'Name="Equate"[^>]*Value="([^"]*)"',
                               m.group(1))
                if eq and eq.group(1) != comp:
                    solvents[comp] = eq.group(1)
        oven = None
        m = re.search(r'SymbolPath value="ColumnOven.Temperature.Nominal"'
                      r'\s*/><Value value="([-\d.]+)', txt)
        if m:
            oven = float(m.group(1))
        return dict(solvents=solvents, flow=flow, gradient=grad,
                    curve=(curve[0] if curve else None), oven_c=oven,
                    t_end=max(t for pts in grad.values()
                              for t, _ in pts))
    return None


def lc_gradient_points(lc, step_sec=2.0):
    """Programmed %B as a chromatographic trace [(rt_sec, pct)] --
    linear interpolation between the gradient table points (HyStar
    Curve 5 = linear). %B missing -> 100 - %A."""
    if not lc:
        return []
    table = lc['gradient'].get('%B')
    if not table:
        table = lc['gradient'].get('%A')
        if not table:
            return []
        table = [(t, 100.0 - v) for t, v in table]
    table = sorted(table)
    t_end = lc['t_end'] * 60.0
    out, i = [], 0
    t = 0.0
    while t <= t_end:
        while i + 1 < len(table) and table[i + 1][0] * 60.0 < t:
            i += 1
        if i + 1 < len(table):
            (t0, v0), (t1, v1) = table[i], table[i + 1]
            f = (t / 60.0 - t0) / (t1 - t0) if t1 > t0 else 0.0
            f = max(0.0, min(1.0, f))
        else:
            v0 = v1 = table[-1][1]
            f = 0.0
        out.append((t, v0 + f * (v1 - v0)))
        t += step_sec
    return out


def lc_gradient_str(lc):
    """Compact one-line gradient for userParams / logs."""
    if not lc or not lc['gradient']:
        return 'n/a'
    # only components that actually vary (%C/%D are usually all-zero)
    comps = sorted(c for c, pts in lc['gradient'].items()
                   if any(abs(v) > 1e-9 for _t, v in pts))
    if not comps:
        comps = sorted(lc['gradient'])[:1]
    per_t = {}
    for c in comps:
        for t, v in lc['gradient'][c]:
            per_t.setdefault(t, {})[c] = v
    return ' ; '.join('%g min: %s' % (
        t, ', '.join('%s %g%%' % (c, per_t[t][c]) for c in comps
                     if c in per_t[t])) for t in sorted(per_t))


def uv_first_wavelength(d_folder, spectra_rows):
    """First acquired UV channel of the .d ABOVE 250 nm: UV1 when it runs
    at 254/280 nm, but UV2 when UV1 sits at 210 nm. Cube source: the first stored wavelength
    above 250 nm. Legacy analysis.sqlite source: the smallest
    wavelength label above 250 nm found in the UV/DAD tables. Falls
    back to the plain first channel when none is above 250 nm.
    Returns [wl_nm] or []."""
    for cand in ('chromatography-data.sqlite',
                 'chromatography-data-pre.sqlite'):
        p = os.path.join(d_folder, cand)
        if os.path.exists(p):
            op = _dad_open(p)
            if op is not None:
                c, wl, _t, _i = op
                c.close()
                above = [round(w, 1) for w in wl if w > 250.0]
                return above[:1] or [round(wl[0], 1)]
    db = os.path.join(d_folder, 'analysis.sqlite')
    if not os.path.exists(db):
        return []
    rt_of = {r.get('Id') for r in spectra_rows}
    wls = []
    try:
        c = sqlite3.connect(db)
        tables = [r[0] for r in c.execute(
            "select name from sqlite_master where type='table'")]
        for t in tables:
            cols = [r[1] for r in c.execute('pragma table_info("%s")' % t)]
            if len(cols) != 3 or 'spectr' not in cols[0].lower():
                continue
            t_is_uv = 'uv' in t.lower() or 'dad' in t.lower() \
                or 'pda' in t.lower()
            try:
                rows = c.execute('select "%s", "%s", "%s" from "%s"'
                                 % (cols[0], cols[1], cols[2], t)).fetchall()
            except sqlite3.Error:
                continue
            for sid, ch, val in rows:
                if sid not in rt_of or val is None:
                    continue
                label = str(ch).strip()
                m = re.search(r'\d+(?:\.\d+)?', label)
                if m is None:
                    continue
                if re.search(r'(?i)(uv|nm)', label) is not None or \
                        (t_is_uv and re.fullmatch(r'\d+(?:\.\d+)?', label)):
                    wls.append(float(m.group(0)))
        c.close()
    except sqlite3.Error:
        return []
    above = sorted(w for w in wls if w > 250.0)
    return [above[0]] if above else ([min(wls)] if wls else [])


def read_uv_scans(d_folder, spectra_rows, wanted):
    """Per-scan UV/DAD wavelength data, if any. Primary source: the
    cube DAD inside chromatography-data(-pre).sqlite (recent otoControl
    format) -- read lambda-targeted (only the `wanted` wavelengths are
    unpacked per row, NOT the whole 611-lambda cube); falls back to the
    legacy 3-column analysis.sqlite heuristic (older formats).
    Returns (scans, cube_db): scans is a time-sorted list
    [(rt_sec, {wavelength_nm: intensity})] ([] when none, then nothing
    changes) and cube_db is the chromatography-data sqlite path when
    the cube source was used (None for the legacy fallback), so the
    caller can pull full apex spectra on demand via _read_dad_full_at.
    """
    for cand in ('chromatography-data.sqlite',
                 'chromatography-data-pre.sqlite'):
        p = os.path.join(d_folder, cand)
        if os.path.exists(p):
            scans = _read_dad_cube(p, wanted)
            if scans:
                return scans, p
    db = os.path.join(d_folder, 'analysis.sqlite')
    if not os.path.exists(db):
        return [], None
    rt_of = {r.get('Id'): r.get('Rt', r.get('RT')) for r in spectra_rows}
    scans = {}
    try:
        c = sqlite3.connect(db)
        tables = [r[0] for r in c.execute(
            "select name from sqlite_master where type='table'")]
        for t in tables:
            cols = [r[1] for r in c.execute('pragma table_info(\"%s\")' % t)]
            if len(cols) != 3:
                continue
            sid_col, ch_col, val_col = cols
            if 'spectr' not in sid_col.lower():
                continue
            t_is_uv = 'uv' in t.lower() or 'dad' in t.lower() \
                or 'pda' in t.lower()
            try:
                rows = c.execute('select "%s", "%s", "%s" from "%s"'
                                 % (sid_col, ch_col, val_col, t)).fetchall()
            except sqlite3.Error:
                continue
            for sid, ch, val in rows:
                label = str(ch).strip()
                m = re.search(r'\d+(?:\.\d+)?', label)
                if m is None or val is None:
                    continue
                wl = float(m.group(0))
                is_uv = re.search(r'(?i)(uv|nm)', label) is not None or \
                    (t_is_uv and re.fullmatch(r'\d+(?:\.\d+)?', label))
                if not is_uv:
                    continue
                rt = rt_of.get(sid)
                if rt is None:
                    continue
                scans.setdefault(float(rt), {})[wl] = float(val)
        c.close()
    except sqlite3.Error:
        return [], None
    return [(rt, wl_int) for rt, wl_int in sorted(scans.items())], None


def uv_traces_from_scans(scans, wavelengths):
    """Chromatographic trace at each requested wavelength (nearest
    available wavelength per scan). Returns {label: [(rt, inten)]}."""
    traces = {}
    for wl in wavelengths:
        pts = []
        for rt, wl_int in scans:
            if not wl_int:
                continue
            nearest = min(wl_int, key=lambda w: abs(w - wl))
            pts.append((rt, wl_int[nearest]))
        traces['%g nm' % wl] = pts
    return traces


def uv_apices(trace, rel_threshold=0.05):
    """Interior local maxima of a UV trace, above rel_threshold of the
    trace maximum. Returns [(rt, intensity)] (chromatographic apexes)."""
    n = len(trace)
    if n < 3:
        return []
    mx = max(v for _, v in trace)
    out = []
    for i in range(1, n - 1):
        rt, v = trace[i]
        if v > trace[i - 1][1] and v >= trace[i + 1][1] \
                and v >= rel_threshold * mx:
            out.append((rt, v))
    return out


def uv_apex_spectra(scans, apex_rts):
    """Full wavelength spectrum of the scan nearest each apex rt."""
    out = []
    for rt_a in apex_rts:
        best = min(scans, key=lambda sc: abs(sc[0] - rt_a))
        out.append((best[0], dict(best[1])))
    return out


def measure_parabene_anchors(specs, polarity):
    """Retention-index standards measured on a converted 'parabens' run.

    For each standard of the series, the most intense recalibrated MS1
    peak within PARABEN_PPM of the expected ion gives the retention-time
    anchor: [M-H]- in negative, [M+H]+ in positive. Returns
    [{name, rt_sec, mz, ri, ion}].
    """
    anchors = []
    for name, mneg, mh, ri in PARABEN_STANDARDS:
        target, ion = (mneg, '[M-H]-') if polarity == 'negative' \
            else (mh, '[M+H]+')
        best = None  # (intensity, rt, mz, ion)
        for sp in specs:
            if sp['ms_level'] != 1:
                continue
            for mz, inten in sp['peaks']:
                if inten > 0 and abs(mz - target) / target * 1e6 \
                        <= PARABEN_PPM:
                    if best is None or inten > best[0]:
                        best = (inten, sp['rt_sec'], mz, ion)
        if best is not None:
            anchors.append(dict(name=name, rt_sec=best[1], mz=best[2],
                                ri=ri, target=target, ion=ion))
    return anchors


def average_parabene_anchors(runs):
    """Mean retention anchors over ALL the parabens runs of one mode:
    several 'parabens' files in the
    folder are all processed, and each anchor name is AVERAGED (rt and
    mz) over the runs where it was found. Returns an anchor list with
    source = 'mean of <runs>'.
    runs: [(file_name, [anchor, ...]), ...]"""
    by_name = {}
    for fname, found in runs:
        for a in found:
            e = by_name.setdefault(a['name'], dict(
                name=a['name'], ri=a['ri'], target=a['target'],
                ion=a['ion'], rts=[], mzs=[], sources=[]))
            e['rts'].append(a['rt_sec'])
            e['mzs'].append(a['mz'])
            e['sources'].append(fname)
    out = []
    for name in sorted(by_name, key=lambda n: by_name[n]['ri']):
        e = by_name[name]
        out.append(dict(
            name=name,
            rt_sec=sum(e['rts']) / len(e['rts']),
            mz=sum(e['mzs']) / len(e['mzs']),
            ri=e['ri'], target=e['target'], ion=e['ion'],
            source=('mean of %s' % ', '.join(e['sources'])
                    if len(e['sources']) > 1 else e['sources'][0])))
    return out


def parse_parabenes(vals):
    """--parabene values ('name=rt' or bare 'rt') -> [(name, rt_min), ...]"""
    out = []
    for i, v in enumerate(vals, 1):
        if '=' in v:
            n, r = v.split('=', 1)
            n = n.strip() or 'parabene_%d' % i
        else:
            n, r = 'parabene_%d' % i, v
        out.append((n, float(r)))
    return out


BAF_MD_ALIASES = (
    # logical field -> Bruker acquisition-parameter names (matched
    # case/space/underscore/minus-insensitively, also by prefix)
    ('rt', ('rt', 'retentiontime', 'scantime')),
    ('segment', ('segment', 'segmentnumber', 'mssegment',
                 'segmentid')),
    ('polarity', ('polarity', 'ionpolarity')),
    # the never-opened DDA baf carries the per-scan precursor
    # in the NUMERIC row 'MSMS_IsolationMass_Act' (probe of
    # 37320: 50.0 on MS1, 181.95/176.97... on MS2) -- the normalized
    # key 'msmsisolationmassact' matches NO old alias (it ends with
    # 'massact', not 'isolationmass'), so every scan decoded MS1.
    ('precursor', ('msmsisolationmassact', 'isolationmass',
                   'isolationtarget', 'targetmass',
                   'msdisolationmass', 'precursor')),
    ('ce', ('collisionenergy', 'collisionenergyev')),
    ('iso', ('isolationwidth', 'isolationwindow', 'isolationmz')),
)


def _md_lookup(low, names):
    """First metadata value matching one of the alias names (exact,
    then by prefix, then by SUFFIX so 'SpectraAcquisition_SegmentId'
    still hits 'segmentid'; the real never-opened runs
    carry SegmentId only under the full 'SpectraAcquisition_*' name,
    which a prefix-only match never caught, so every scan decoded
    without a Segment and the calibrant gate never matched)."""
    for n in names:
        if n in low:
            return low[n]
    for k, v in low.items():
        for n in names:
            if k.startswith(n) or k.endswith(n):
                return v
    return None


def _md_num_from_string(raw, lo=0.0, hi=1e9):
    """first float in [lo, hi] decoded from a STRING/list
    metadata row. The never-opened baf carries the per-scan MS/MS values under
    list-typed names (e.g. 'MSMS_ListIsolationMass') whose value is a
    string like '641.1715' or '641.1715;656.1891' -- the
    numeric-only rule (which fixed an int() crash on binary
    blobs) also silently dropped these, so EVERY scan of a
    never-opened DDA run decoded as MS1 and the report showed n_ms2=0
    (batch 261106: isolates A1/A2/B2, blanks and most 100-200 eV
    standards lost their MSn this way). Decode instead of dropping;
    returns None when no positive number is embedded."""
    if raw is None or isinstance(raw, (int, float, bool)):
        return None
    try:
        txt = raw.decode('utf-8', 'replace') if isinstance(
            raw, (bytes, bytearray)) else str(raw)
    except Exception:
        return None
    for tok in txt.replace(',', ';').replace(' ', ';').split(';'):
        tok = tok.strip()
        if ':' in tok:                       # Bruker 'code:label' rows
            tok = tok.split(':', 1)[0].strip()
        try:
            v = float(tok)
        except ValueError:
            continue
        if lo <= v <= hi:
            return v
    return None


def load_spectra_table_baf(d_folder, log=print):
    """(spectra, steps, psv) rebuilt from the raw analysis.baf for a
    .d folder that was NEVER OPENED in Bruker DataAnalysis (no
    analysis.sqlite -- the fallback: analysis.sqlite
    is only written at the FIRST
    DataAnalysis opening; the never-opened acquisitions carry the
    whole scan list in the baf, which holds one metadata section PER
    SCAN with the values of the acquisition-header parameters -- the
    sqlite tables are only DataAnalysis's rendering of these blocks).

    The alias table is deliberately conservative: a field the baf
    does not name as expected simply degrades (no precursor -> the
    scan decodes as MS1; no Segment -> no calibrant-segment gate ->
    deferred, calibration borrowed from a same-mode neighbour), it
    never corrupts the conversion. When a field is missing the log
    lists the parameter names actually seen so --probe on a real
    never-opened .d can extend the aliases in one step."""
    baf = open_baf(d_folder)
    header_off, sections = read_index(
        os.path.join(d_folder, 'analysis.baf_idx'))
    if not sections:
        raise ValueError('no scan section in the baf index of %s'
                         % d_folder)
    params = _read_baf_params(baf, header_off) if header_off is not None \
        else {}

    def quiet(*_a, **_k):
        pass

    spectra, steps, psv = [], {}, {}
    no_rt, keys_seen, first_block = [], set(), None
    # acquisition-mode COVERAGE signals: d2mz
    # is validated on DDA CID (constant CE or CE ramps); IsCID energy
    # or a non-CID fragmentation mode must produce a WARNING instead
    # of a silently-converted, out-of-scope run.
    _modes = {'iscid': 0, 'fragmode': set(), 'iso_act': 0}
    str_prec = {}          # scans whose precursor came
                            # from a string/list metadata row
    n_meta = n_rt_head = 0
    for flag in sorted(sections):
        moff = sections[flag].get('metadata')
        if moff is None:
            continue
        n_meta += 1
        # the per-scan Rt of never-opened baf files lives in
        # the metadata-block HEADER, not in the parameter rows. The
        # diff (scan 1 vs scan 3970 of the real Qualitative
        # file) showed the i32 at header offset 24 = 299 vs 965142:
        # retention time in MILLISECONDS (0.299 s -> 965.142 s =
        # 16.09 min, 243.2 ms per scan, consistent with the ~244 ms
        # Digitizer spectra sample time). Guarded by mhs >= 28 (the
        # legacy 12-byte mock headers keep the alias fallback) and a
        # 48 h sanity window.
        rt_sec = None
        try:
            with open(baf, 'rb') as _f:
                _f.seek(moff)
                _mbs, _tag, _mhs = struct.unpack('<i4si', _f.read(12))
                if _mhs >= 28:
                    _f.seek(moff + 24)
                    _rtms = struct.unpack('<i', _f.read(4))[0]
                    if -3600000 <= _rtms <= 172800000:
                        rt_sec = _rtms / 1000.0
        except (OSError, struct.error):
            rt_sec = None
        vals = _read_metadata_values(baf, moff, params, quiet)
        if first_block is None:
            first_block = vals
        low, low_num = {}, {}
        for name, raw in vals.items():
            k = str(name).strip().lower().replace(' ', '').replace(
                '_', '').replace('-', '')
            low[k] = raw
            # numeric fields (rt/segment/precursor/ce/iso) may
            # only match NUMERIC rows -- the old suffix match let
            # string/list rows through (MSMS_ListIsolationMass ends
            # with 'isolationmass') and int()/float() then crashed
            # on the raw binary value ('invalid literal for int()')
            if isinstance(raw, (int, float)) and not isinstance(
                    raw, bool):
                low_num[k] = raw
            keys_seen.add(k)
        _alias = dict(BAF_MD_ALIASES)
        # polarity is a genuine string row ('positive'/'negative')
        got = {'polarity': _md_lookup(low, _alias['polarity'])}
        for f in ('rt', 'segment', 'precursor', 'ce', 'iso'):
            got[f] = _md_lookup(low_num, _alias[f])
        # never-opened DDA files carry the MS/MS values in
        # string/list rows -- decode them when no numeric row matched
        str_src = {}
        for f in ('precursor', 'ce', 'iso'):
            if got[f] is not None:
                continue
            for k, raw in low.items():
                if k in low_num:
                    continue          # numeric rows already tried
                hit = None
                for n in _alias[f]:
                    if k == n or k.startswith(n) or k.endswith(n):
                        hit = k
                        break
                if hit is None:
                    continue
                _rng = {'precursor': (20.0, 5000.0),
                        'ce': (0.5, 1000.0),
                        'iso': (0.001, 500.0)}[f]
                v = _md_num_from_string(raw, *_rng)
                if v is not None:
                    got[f] = v
                    str_src[f] = hit
                    break
        # gate: 'MSMS_Isolation_Act' is a per-scan NUMERIC flag
        # (probe of 37320: 0 on the MS1 scan 1, 1 on the MS2 scan
        # 948). Without this gate, the 50.0 DEFAULT that sits in
        # 'MSMS_IsolationMass_Act' on every MS1 scan would turn the
        # whole run into MS2. When the flag is absent (older files,
        # sqlite-style blocks) the previous behaviour is unchanged.
        _iso_act = _md_lookup(low_num, ('msmsisolationact',))
        if _iso_act is not None:
            try:
                if float(_iso_act) == 0.0:
                    got['precursor'] = None
                else:
                    _modes['iso_act'] += 1
            except (TypeError, ValueError):
                pass
        # mode-coverage signals
        for _k in ('transferiscideenergyact', 'transferiscidenergy',
                   'iscidenergy'):
            _v = low_num.get(_k)
            if _v is not None and float(_v) > 0.0:
                _modes['iscid'] += 1
                break
        _fm = low_num.get('msmsfragmentationmodeact')
        if _fm is None:
            _fm = low_num.get('fragmentationmodeact')
        if _fm is not None and float(_fm) != 0.0:
            _modes['fragmode'].add(float(_fm))
        rt_val = rt_sec if rt_sec is not None else got['rt']
        if rt_val is None:
            no_rt.append(flag)
            continue
        if rt_sec is not None:
            n_rt_head += 1
        # belt and braces -- a bad Rt value degrades to "no
        # Rt for this scan" instead of killing the whole conversion
        try:
            rt_f = float(rt_val)
        except (TypeError, ValueError):
            no_rt.append(flag)
            continue
        row = {'Id': flag, 'Rt': rt_f}
        if got['segment'] is not None:
            try:
                row['Segment'] = int(got['segment'])
            except (TypeError, ValueError):
                pass
        if got['polarity'] is not None:
            row['Polarity'] = got['polarity']
        spectra.append(row)
        if got['precursor'] is not None:
            try:
                m = float(got['precursor'])
            except (TypeError, ValueError):
                m = None
            if m and m > 0:
                steps[flag] = [{'Mass': m}]
                if str_src.get('precursor'):
                    str_prec[flag] = m
        pv = {}
        for key, idx in (('ce', 5), ('iso', 7)):
            if got[key] is not None:
                try:
                    pv[idx] = float(got[key])
                except (TypeError, ValueError):
                    pass
        if pv:
            psv[flag] = pv
    # guard: a method-level MS/MS target list (one identical
    # string precursor on ~EVERY scan, MS1 scans included) must not
    # turn the whole run into MS2. A real DDA acquisition varies the
    # precursor scan to scan and/or leaves the MS1 scans without one.
    if str_prec:
        _cnt = {}
        for v in str_prec.values():
            _cnt[v] = _cnt.get(v, 0) + 1
        _modal, _n = max(_cnt.items(), key=lambda kv: kv[1])
        if _n >= 0.98 * len(spectra) and len(spectra) > 20:
            log('    WARNING: string-row precursor %r identical '
                'on %d/%d scans -- suspected METHOD-LEVEL MS/MS list, '
                'not applied (run d2mz --probe to inspect the rows)'
                % (_modal, _n, len(spectra)))
            for flag in list(str_prec):
                steps.pop(flag, None)
        else:
            log('    %d/%d scan(s) given a precursor by '
                'string/list metadata rows (%s...)'
                % (len(str_prec), len(spectra),
                   sorted(set(str_prec.values()))[:6]))
    if not spectra:
        # SELF-DIAGNOSING failure -- the error (it lands in the
        # per-file report and the NOT converted recap) carries the
        # counts, the acquisition-header status and the parameter
        # names+values of the FIRST metadata block, so one rerun is
        # enough to extend BAF_MD_ALIASES without a --probe round-trip
        # A 12-item sample showed only instrument
        # currents (Source/TOF/Digitizer...) on the real never-opened
        # baf -- the WHOLE first block is now dumped (all metadata
        # blocks share the same parameter table, so this is the
        # complete name list of the file)
        # A full dump CONFIRMED the real metadata blocks
        # carry ONLY instrument status (38 params, none Rt-named), so
        # the scan table must live elsewhere -- the error now also
        # dumps the RAW BYTES of (1) the first baf_idx record, (2)
        # the first spectrum section (1,16), (3) the first metadata
        # block header, so ONE rerun gives the ground truth to
        # locate the per-scan Rt without any probe round-trip
        # A rerun decoded idx0 as the baf_idx FILE header
        # (not a scan record) and md0 as a sane block header (mbs=758,
        # mhs=40, ncnt=35 numeric rows) -- no Rt in either, and spec0
        # was ABSENT (no (1,16) record in the index at all). The diag
        # is extended to the ground-truth diff: (a) census of the
        # section-id pairs actually present in the baf_idx, (b) the
        # RAW (2,32) index records of the FIRST and LAST scan (the
        # unused bytes around flag/off are the prime Rt candidates),
        # (c) the metadata-block headers of the FIRST and LAST scan
        # (a field growing ~0.24 s per scan is the Rt), (d) the whole
        # acquisition-header parameter table with its row TYPE codes
        # (the denormal values = integers read as doubles,
        # so the numeric rows are typed and the type lives here).
        # Values in the sample are truncated so one bogus string row
        # cannot swallow megabytes of error text again.
        _items = list((first_block or {}).items())

        def _short(v):
            s = repr(v)
            return s if len(s) <= 96 else s[:93] + '...'
        _sample = ', '.join('%s=%s' % (k, _short(v))
                            for k, v in _items)
        _diag = []
        try:
            _vectors = get_raw_vectors(
                os.path.join(d_folder, 'analysis.baf_idx'))
            _census = {}
            for p in _vectors:
                _census[(p[0], p[1])] = _census.get((p[0], p[1]), 0) + 1
            _diag.append('ids=%s' % ','.join(
                '%02x%02x:%d' % (a, b, n)
                for (a, b), n in sorted(_census.items())))
            if _vectors:
                _diag.append('idx0=%s' % _vectors[0][:32].hex())
            _mrec = {}
            for p in _vectors:
                if len(p) >= 24 and (p[0], p[1]) == (2, 32):
                    _mrec[struct.unpack_from('<i', p, 8)[0]] = p
            if _mrec:
                _f0 = min(_mrec)
                _f1 = max(_mrec)
                _diag.append('idxrec%d=%s' % (_f0, _mrec[_f0][:32].hex()))
                _diag.append('idxrec%d=%s' % (_f1, _mrec[_f1][:32].hex()))
        except (OSError, EOFError, struct.error):
            pass
        try:
            _so = sections[sorted(sections)[0]].get('spectrum')
            if _so is not None:
                with open(baf, 'rb') as _f:
                    _f.seek(_so)
                    _diag.append('spec0=%s' % _f.read(64).hex())
        except (OSError, KeyError, IndexError):
            pass
        try:
            _moffs = [sections[_k]['metadata'] for _k in sorted(sections)
                      if 'metadata' in sections[_k]]
            if _moffs:
                with open(baf, 'rb') as _f:
                    _f.seek(_moffs[0])
                    _diag.append('md0=%s' % _f.read(64).hex())
                    _f.seek(_moffs[-1])
                    _diag.append('mdN=%s' % _f.read(64).hex())
        except (OSError, KeyError, IndexError):
            pass
        try:
            if params:
                _diag.append('ptab=%s' % ';'.join(
                    '%d/%d/%s' % (ri, t,
                                  (f[0].decode('utf-8', 'replace')
                                   if f[0] else '')[:28])
                    for ri, (t, f) in sorted(params.items())))
        except Exception:
            pass
        raise ValueError(
            'analysis.baf metadata: no Rt-named parameter (%d scan '
            'section(s), %d metadata block(s), acquisition header %s; '
            'first block (%d params): %s; diag %s)'
            % (len(sections), n_meta,
               'found' if header_off is not None else 'MISSING',
               len(_items), _sample or 'empty',
               ' '.join(_diag) or 'unavailable'))
    log('    spectra table rebuilt from the analysis.baf metadata '
        '(never opened in DataAnalysis): %d scan(s), %d with a '
        'declared precursor' % (len(spectra), len(steps)))
    if n_rt_head:
        log('    Rt source: metadata-block header word at offset +24 '
            '(milliseconds) for %d scan(s)' % n_rt_head)
    if no_rt:
        log('    NOTE: %d scan section(s) without an Rt metadata value '
            'skipped' % len(no_rt))
    if not steps:
        log('    NOTE: no precursor parameter matched in the baf '
            'metadata (names seen: %s) -> every scan decodes as MS1; '
            'run --probe on this folder and send the "header params + '
            'metadata" section so the alias table can be extended'
            % ', '.join(sorted(keys_seen)[:12]))
    psv['_modes'] = {
        'iscid': _modes['iscid'],
        'fragmode': sorted(_modes['fragmode']),
        'iso_act': _modes['iso_act']}
    return spectra, steps, psv


def load_spectra_table(d, log=print):
    sq = os.path.join(d, 'analysis.sqlite')
    if not os.path.exists(sq):
        # never opened in DataAnalysis: rebuild from the raw baf
        return load_spectra_table_baf(d, log)
    c = sqlite3.connect(sq)
    c.row_factory = sqlite3.Row
    spectra = [dict(r) for r in c.execute('select * from Spectra order by Id')]
    cols = {r[1] for r in c.execute('pragma table_info(Steps)')}
    tgt = 'TargetSpectrum' if 'TargetSpectrum' in cols else \
          ('Spectrum' if 'Spectrum' in cols else None)
    if tgt is None:
        c.close()
        raise ValueError('Steps table: no TargetSpectrum/Spectrum column')
    steps = {}
    for r in c.execute('select * from Steps'):
        steps.setdefault(r[tgt], []).append(dict(r))
    psv = {}
    try:
        for r in c.execute(
                'select Spectrum, Variable, Value from PerSpectrumVariables'):
            psv.setdefault(r['Spectrum'], {})[r['Variable']] = r['Value']
    except sqlite3.Error:
        pass
    c.close()
    return spectra, steps, psv

# ===================================================== calibration & re-phasing

def ppm_error(observed, theoretical):
    return (observed - theoretical) / theoretical * 1e6


LOCK_MIN_CO = 5   # min cluster masses that must CO-MATCH in the same
                  # MS1 scan for it to count as a calibrant observation:
                  # the calibrant segment co-matches 9-14 masses, while a
                  # dense analyte scan only hits 1-3 by coincidence inside
                  # the wide 30 ppm lock window (Impact pollen extracts
                  # saturated the matcher at 1019/1019 scans before this)


class LockMassCalibrator:
    """Quadratic ppm-drift model vs retention time, fitted on lock-mass
    peaks (2nd-order polynomial, like the qTOF mass calibration).

    For every MS1 scan the most intense peak within --lock-ppm of each
    theoretical lock mass is measured; the per-scan median ppm error is
    regressed against RT and removed from every spectrum (MS1 and MS2).
    """

    def __init__(self, lock_masses, lock_ppm, min_points=MIN_LOCK_POINTS):
        self.lock_masses = list(lock_masses)
        self.lock_ppm = lock_ppm
        self.min_points = min_points
        self.points = []      # (rt_sec, median_ppm)
        self.model = None
        self.borrowed_from = None

    def observe_scan(self, rt, peaks):
        errs, tol = [], self.lock_ppm * 1e-6
        for lm in self.lock_masses:
            best = None
            for mz, inten in peaks:
                if abs(mz - lm) <= lm * tol:
                    if best is None or inten > best[1]:
                        best = (mz, inten)
            if best is not None:
                errs.append(ppm_error(best[0], lm))
        if len(errs) >= LOCK_MIN_CO:
            # calibrant scan: SEVERAL cluster masses co-match (comb);
            # single-mass coincidences from dense analyte scans are
            # rejected so the ppm(rt) model is fitted on the calibrant
            self.points.append((rt, sorted(errs)[len(errs) // 2]))

    def fit(self):
        """Least-squares fit of a QUADRATIC drift model vs RT:

            ppm(rt) = c0 + c1*rt + c2*rt^2

        The regressed quantity is the deviation from the expected
        lock-mass position -- the same quantity Bruker schematizes with
        a quadratic regression on its qTOFs (on a TOF, m/z scales with
        the square of the flight time, so mass deviations are largest
        at low masses).
        """
        pts = sorted(self.points)
        if len(pts) < self.min_points:
            self.model = None
            return None
        rts = [p[0] for p in pts]
        ppms = [p[1] for p in pts]
        # normal equations for the degree-2 polynomial (pure stdlib, 3x3)
        s = [sum(x ** k for x in rts) for k in range(5)]
        t = [sum(y * x ** k for x, y in zip(rts, ppms)) for k in range(3)]
        A = [[s[0], s[1], s[2], t[0]],
             [s[1], s[2], s[3], t[1]],
             [s[2], s[3], s[4], t[2]]]
        # Gaussian elimination with partial pivoting
        for col in range(3):
            piv = max(range(col, 3), key=lambda r: abs(A[r][col]))
            if abs(A[piv][col]) < 1e-30:
                # degenerate: fall back to the mean
                self.model = (sum(ppms) / len(ppms), 0.0, 0.0)
                return self.model
            A[col], A[piv] = A[piv], A[col]
            for r in range(col + 1, 3):
                f = A[r][col] / A[col][col]
                for c in range(col, 4):
                    A[r][c] -= f * A[col][c]
        coeffs = [0.0, 0.0, 0.0]
        for r in (2, 1, 0):
            acc = A[r][3] - sum(A[r][c] * coeffs[c] for c in range(r + 1, 3))
            coeffs[r] = acc / A[r][r]
        self.model = tuple(coeffs)      # (c0, c1, c2)
        return self.model

    def rms_before(self):
        if not self.points:
            return None
        v = [p[1] for p in self.points]
        return math.sqrt(sum(x * x for x in v) / len(v))

    def rms_after(self):
        if not self.points or self.model is None:
            return None
        c0, c1, c2 = self.model
        r = [p[1] - (c0 + c1 * p[0] + c2 * p[0] ** 2) for p in self.points]
        return math.sqrt(sum(x * x for x in r) / len(r))

    def mean_ppm_after(self):
        """Mean lock-mass error (ppm) after calibration."""
        if not self.points or self.model is None:
            return None
        c0, c1, c2 = self.model
        r = [p[1] - (c0 + c1 * p[0] + c2 * p[0] ** 2) for p in self.points]
        return sum(r) / len(r)

    def stdev_ppm_after(self):
        """Standard deviation of the lock-mass errors (ppm) after
        calibration — the quantity Bruker's mass-accuracy specification
        is expressed on (mean standard deviation < 2 ppm)."""
        if not self.points or self.model is None:
            return None
        c0, c1, c2 = self.model
        r = [p[1] - (c0 + c1 * p[0] + c2 * p[0] ** 2) for p in self.points]
        mu = sum(r) / len(r)
        return math.sqrt(sum((x - mu) ** 2 for x in r) / len(r))

    def correction_ppm(self, rt):
        if self.model is None:
            return 0.0
        c0, c1, c2 = self.model
        return c0 + c1 * rt + c2 * rt * rt

    def apply(self, rt, peaks):
        c = self.correction_ppm(rt) * 1e-6
        if c == 0.0:
            return list(peaks)
        return [(mz * (1.0 - c), inten) for mz, inten in peaks]


def reassign_precursor(declared_mz, ms1_peaks, window_da):
    """Re-phase an MS2 precursor on the recalibrated MS1 apex.

    The search window is the isolation window (in Da) — a physical bound,
    not a user-tuned tolerance. Within it, the most intense recalibrated
    MS1 apex is taken (deterministic tie-break on |mz - setpoint|).

    Returns (new_mz, intensity, delta_ppm) or (None, None, None) when no
    MS1 peak falls within the isolation window.
    """
    if not ms1_peaks:
        return None, None, None
    cand = sorted(((mz, inten) for mz, inten in ms1_peaks
                   if abs(mz - declared_mz) <= window_da),
                  key=lambda t: (-t[1], abs(t[0] - declared_mz)))
    if not cand:
        return None, None, None
    new_mz, inten = cand[0]
    return new_mz, inten, ppm_error(new_mz, declared_mz)

# ============================================================ binary encoding

def encode_binary(values, precision=32, compress=True, big_endian=False):
    end = '>' if big_endian else '<'
    fmt = '%s%d%s' % (end, len(values), 'f' if precision == 32 else 'd')
    raw = struct.pack(fmt, *values)
    if compress:
        raw = zlib.compress(raw, 6)
    return base64.b64encode(raw).decode('ascii')

# ================================================================== mzML writer

def _spectrum_mzml(sp, opts):
    peaks = sp['peaks']
    n = len(peaks)
    mz = [p[0] for p in peaks]
    it = [p[1] for p in peaks]
    tic = float(sum(it))
    bpm = max(it) if it else 0.0
    bpv = mz[it.index(bpm)] if it and bpm > 0 else 0.0
    o = []
    o.append('<spectrum index="%d" id="scan=%d" defaultArrayLength="%d">\n'
             % (sp['scan_index'], sp['scan_index'], n))
    o.append('<cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="%d"/>\n'
             % sp['ms_level'])
    o.append('<cvParam cvRef="MS" accession="MS:1000127" '
             'name="centroid spectrum"/>\n')
    o.append('<cvParam cvRef="MS" accession="MS:%s" name="%s"/>\n'
             % ('1000130' if sp['polarity'] == 'positive' else '1000129',
                'positive scan' if sp['polarity'] == 'positive' else 'negative scan'))
    o.append('<cvParam cvRef="MS" accession="MS:1000285" name="total ion current" value="%.6g"/>\n'
             % tic)
    if bpm:
        o.append('<cvParam cvRef="MS" accession="MS:1000504" name="base peak m/z" value="%.6g" '
                 'unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>\n' % bpv)
        o.append('<cvParam cvRef="MS" accession="MS:1000505" name="base peak intensity" value="%.6g"/>\n'
                 % bpm)
    o.append('<scanList count="1"><scan>')
    o.append('<cvParam cvRef="MS" accession="MS:1000795" name="no combination"/>')
    o.append('<cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value="%.6f" '
             'unitCvRef="UO" unitAccession="UO:0000031" unitName="unit minute"/>'
             '</scan></scanList>\n' % (sp['rt_sec'] / 60.0))
    if sp['ms_level'] >= 2 and sp.get('precursor_mz') is not None:
        o.append('<precursorList count="1">')
        o.append('<precursor spectrumRef="scan=%d">' % (sp.get('precursor_scan') or 0))
        o.append('<isolationWindow><cvParam cvRef="MS" accession="MS:1000827" '
                 'name="isolation target m/z" value="%.10g" '
                 'unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/></isolationWindow>' % sp['precursor_mz'])
        o.append('<selectedIonList count="1"><selectedIon>')
        o.append('<cvParam cvRef="MS" accession="MS:1000744" name="selected ion m/z" value="%.10g" '
                 'unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>' % sp['precursor_mz'])
        if sp.get('precursor_intensity'):
            o.append('<cvParam cvRef="MS" accession="MS:1000427" name="peak intensity" value="%.6g"/>'
                     % sp['precursor_intensity'])
        o.append('</selectedIon></selectedIonList><activation>')
        o.append('<cvParam cvRef="MS" accession="MS:1000044" name="collision-induced dissociation"/>')
        if sp.get('ce') is not None:
            o.append('<cvParam cvRef="MS" accession="MS:1000045" name="collision energy" value="%.4g" '
                     'unitCvRef="UO" unitAccession="UO:0000266" unitName="electronvolt"/>' % sp['ce'])
            w = ce_window(sp['ce'])
            if w:
                o.append('<userParam name="d2mz:ce_window" value=%s/>'
                         % quoteattr('%s eV' % w))
            # ramp vs constant, from the .method sweep switch
            if sp.get('ce_label'):
                o.append('<userParam name="d2mz:collision_energy" value=%s/>'
                         % quoteattr(sp['ce_label']))
        o.append('</activation></precursor></precursorList>\n')
    mz_enc = encode_binary(mz, opts['precision'], opts['compress'])
    it_enc = encode_binary(it, opts['precision'], opts['compress'])
    dtype = '1000521" name="32-bit float' if opts['precision'] == 32 else '1000523" name="64-bit float'
    comp = '1000574" name="zlib compression' if opts['compress'] else '1000576" name="no compression'
    o.append('<binaryDataArrayList count="2">')
    # units on the array-type cvParams (msconvert-style;
    # strict readers such as MZmine's chromatogram path expect them)
    for acc, unit, enc in (
            ('1000514" name="m/z array', 'MS:1000040" unitName="m/z', mz_enc),
            ('1000515" name="intensity array', 'MS:1000131" unitName="intensity unit', it_enc)):
        o.append('<binaryDataArray encodedLength="%d">'
                 '<cvParam cvRef="MS" accession="MS:%s"/>'
                 '<cvParam cvRef="MS" accession="MS:%s"/>'
                 '<cvParam cvRef="MS" accession="MS:%s" unitCvRef="MS" unitAccession="%s"/>'
                 '<binary>%s</binary></binaryDataArray>'
                 % (len(enc), comp, dtype, acc, unit, enc))
    o.append('</binaryDataArrayList>\n</spectrum>\n')
    return ''.join(o)


def write_mzml(path, spectra, opts, source_name, start_stamp=None,
               uv_traces=None, uv_apexes=None, retention_standards=None,
               extra_traces=None, lc_method=None):
    parts = []
    add = parts.append
    add('<?xml version="1.0" encoding="utf-8"?>\n')
    add('<indexedmzML xmlns="http://psi.hupo.org/ms/mzml">\n')
    add('<mzML xmlns="http://psi.hupo.org/ms/mzml" version="1.1.0">\n')
    add('<cvList><cv id="MS" fullName="Proteomics Standards Initiative Mass '
        'Spectrometry Ontology" '
        'URI="http://psidev.cvs.sourceforge.net/psidev/psi/psi-ms/ms-MsOntology.obo" '
        'version="4.1.0"/></cvList>\n')
    add('<fileDescription><fileContent>'
        '<cvParam cvRef="MS" accession="MS:1000580" name="MS1 spectrum"/>'
        '</fileContent><sourceFileList count="1">'
        '<sourceFile id="sf1" name=%s location="file:///">'
        '<cvParam cvRef="MS" accession="MS:1000768" name="Bruker BAF file"/>'
        '<cvParam cvRef="MS" accession="MS:1000824" name="no combination"/>'
        '</sourceFile></sourceFileList></fileDescription>\n'
        % quoteattr(source_name))
    add('<softwareList count="2">'
        '<software id="d2mz" version="%s"><cvParam cvRef="MS" accession="MS:1000799" '
        'name="custom software"/></software>'
        '<software id="cvt" version="1.0"><cvParam cvRef="MS" accession="MS:1000799" '
        'name="custom software"/></software></softwareList>\n' % __version__)
    add('<instrumentConfigurationList count="1">'
        '<instrumentConfiguration id="ic1">'
        '<cvParam cvRef="MS" accession="MS:1000443" name="mass analyzer type"/>'
        '<componentList>'
        '<source order="1"><cvParam cvRef="MS" accession="MS:1000073" '
        'name="electrospray ionization"/></source>'
        '<analyzer order="2"><cvParam cvRef="MS" accession="MS:1000084" '
        'name="time-of-flight"/></analyzer>'
        '<detector order="3"><cvParam cvRef="MS" accession="MS:1000251" '
        'name="detector type"/></detector></componentList>'
        '<softwareRef ref="cvt"/></instrumentConfiguration>'
        '</instrumentConfigurationList>\n')
    dp = ['<dataProcessingList count="1"><dataProcessing id="dp1">']
    # the centroids are Bruker's own peak picking (firmware 'maxima'
    # blocks, both MS levels) — declare it like msconvert's
    # 'PeakPicking: vendor, MS-levels 1-2'
    dp.append('<processingMethod order="0" softwareRef="d2mz">'
              '<cvParam cvRef="MS" accession="MS:1000035" '
              'name="peak picking"/>'
              '<userParam name="PeakPicking" value="vendor, '
              'MS-levels 1-2 (Bruker firmware maxima blocks)"/>'
              '</processingMethod>')
    lock_series = ('built-in Bruker sodium formate/acetate cluster series '
                   '(38 positive / 36 negative ions; the series matching the '
                   'detected polarity is locked on)')
    dp.append('<processingMethod order="1" softwareRef="d2mz">'
              '<cvParam cvRef="MS" accession="MS:1000544" name="Conversion '
              'to mzML"/>'
              '<userParam name="d2mz:lock_masses" value=%s/>'
              '<userParam name="d2mz:precursor_rephasing" value=%s/>'
              '</processingMethod></dataProcessing></dataProcessingList>\n')
    add(''.join(dp)
        % (quoteattr(lock_series),
           quoteattr('on' if opts['prec_correct'] else 'off')))
    add('<run id=%s defaultInstrumentConfigurationRef="ic1" '
        'defaultDataProcessingRef="dp1" startTimeStamp=%s>\n'
        % (quoteattr(os.path.splitext(source_name)[0]),
           quoteattr(start_stamp or '1970-01-01T00:00:00Z')))
    # citation one-liner for distribution -- the script
    # version was only buried in softwareList and there was NO
    # conversion date anywhere (startTimeStamp is the ACQUISITION
    # start). One explicit userParam, readable by anyone.
    add('<userParam name="d2mz:conversion" value=%s/>\n'
        % quoteattr('converted with d2mz.py %s on %s UTC '
                    '(Stéphane GREFF, IMBE, Marseille -- BEEMED project, '
                    'ANSES PNR EST 22-EST-063)'
                    % (__version__,
                       time.strftime('%Y-%m-%dT%H:%M:%S', time.gmtime()))))
    # Validated acquisition-mode declaration --
    # the d2mz pipeline is validated on DDA CID (constant CE or CE
    # ramps); anything else is flagged here AND in the report.
    if opts.get('acquisition_mode'):
        add('<userParam name="d2mz:acquisition_mode" value=%s/>\n'
            % quoteattr(opts['acquisition_mode']))
        if opts.get('mode_warning'):
            add('<userParam name="d2mz:mode_validation" value=%s/>\n'
                % quoteattr('OUTSIDE the validated scope: %s '
                            '(validated: DDA CID constant-CE / '
                            'CE-ramp)' % opts['mode_warning']))
    if lc_method:
        # LC program provenance (from hystar.method) --
        # solvents, flow, gradient table, oven setpoint; one readable
        # userParam each, ignored transparently by MS software.
        if lc_method.get('solvents'):
            add('<userParam name="d2mz:lc_solvents" value=%s/>\n'
                % quoteattr('; '.join(
                    '%s = %s' % kv
                    for kv in sorted(lc_method['solvents'].items()))))
        if lc_method.get('flow'):
            add('<userParam name="d2mz:lc_flow" value=%s/>\n'
                % quoteattr(', '.join('%g min: %g ml/min' % fv
                                      for fv in lc_method['flow'])))
        add('<userParam name="d2mz:lc_gradient" value=%s/>\n'
            % quoteattr('%s%s' % (
                lc_gradient_str(lc_method),
                ' (curve %g)' % lc_method['curve']
                if lc_method.get('curve') else '')))
        if lc_method.get('oven_c') is not None:
            add('<userParam name="d2mz:lc_column_oven" value=%s/>\n'
                % quoteattr('setpoint %.1f C' % lc_method['oven_c']))
    for a in (uv_apexes or []):
        # full UV spectrum at a chromatographic apex, as plain text the
        # mzML carries for ever -- ignored by MS software, readable by
        # anyone (wl:intensity pairs, comma-separated)
        data = ','.join('%s:%.6f' % (
                           wl if isinstance(wl, str) else '%g' % wl, v)
                        for wl, v in sorted(
                            a['spectrum'].items(),
                            key=lambda kv: str(kv[0])))
        add('<userParam name="d2mz:uv_spectrum" value=%s/>\n'
            % quoteattr('trace=%s; rt=%.4f min; data=%s'
                         % (a['trace'], a['rt_sec'] / 60.0, data)))
    stds = {}   # name -> value string; manual entries win by name
    for a in (retention_standards or []):
        stds[a['name']] = (
            'rt %.4f min; ri %d; mz %.4f (%s); source %s'
            % (a['rt_sec'] / 60.0, a['ri'], a['mz'], a['ion'], a['source']))
    for pname, prt in (opts.get('parabenes') or []):
        stds[pname] = 'rt %.4f min' % prt
    for pname in sorted(stds):
        # run-level userParams: MS software transparently ignores them
        add('<userParam name="d2mz:retention_standard:%s" value=%s/>\n'
            % (escape(pname), quoteattr(stds[pname])))
    add('<spectrumList count="%d" defaultDataProcessingRef="dp1">\n' % len(spectra))
    offsets = []
    for k, sp in enumerate(spectra, 1):
        _progress(k, len(spectra), 'writing mzML')
        offsets.append((sp['scan_index'],
                        sum(len(p) for p in parts)))
        add(_spectrum_mzml(sp, opts))
    _progress_end()
    add('</spectrumList>\n')
    # UV/DAD traces as mzML chromatograms: readers that do not request
    # chromatograms skip this list entirely. The OTHER existing traces
    # of the acquisition (2nd UV channel, temperatures, pressure, MS
    # base peak...) are appended after the analytical UV trace
    # (every existing trace goes into the mzML)
    uv = sorted((uv_traces or {}).items())
    extras = [(ds, wl, pts) for ds, wl, pts in (extra_traces or [])
              if pts]
    chrom_offsets = []
    if uv or extras:
        add('<chromatogramList count="%d" defaultDataProcessingRef="dp1">\n'
            % (len(uv) + len(extras)))

        def _chrom_arrays(pts):
            # MZmine's multi-detector import mis-scales the
            # chromatogram time axis by x1000 when the unit is 'second'
            # (validated empirically 2026-11-07: 60-720 s displayed as
            # 0-720,000; the same values in MINUTES with UO:0000031
            # display correctly). mzML allows any time unit; we now
            # write chromatogram times in minutes, like scan start time.
            times = [p[0] / 60.0 for p in pts]
            ints = [p[1] for p in pts]
            end = 'f' if opts['precision'] == 32 else 'd'
            add('<binaryDataArrayList count="2">\n')
            # MZmine's strict chromatogram parser
            # (MZMLParser.java:643) requires a unitAccession on the
            # array-type cvParams -- msconvert writes MS:1000131
            # "intensity unit" on the intensity array; ours had none
            # -> "Tag cvParam must provide an 'unitAccession'
            # attribute". Both arrays now carry full units.
            for acc, nm, unit, unit_cref, unit_acc, arr in (
                    ('MS:1000595', 'time array', 'minute', 'UO',
                     'UO:0000031', times),
                    ('MS:1000515', 'intensity array', 'intensity unit',
                     'MS', 'MS:1000131', ints)):
                raw = struct.pack('<%d%s' % (len(arr), end), *arr)
                comp = ('<cvParam cvRef="MS" accession="MS:1000574" '
                        'name="zlib compression"/>' if opts['compress'] else
                        '<cvParam cvRef="MS" accession="MS:1000576" '
                        'name="no compression"/>')
                prec = ('<cvParam cvRef="MS" accession="MS:1000521" '
                        'name="32-bit float"/>' if end == 'f' else
                        '<cvParam cvRef="MS" accession="MS:1000523" '
                        'name="64-bit float"/>')
                unit = ((' unitCvRef="%s" unitAccession="%s" '
                         'unitName="%s"' % (unit_cref, unit_acc, unit))
                        if unit_acc else '')
                if opts['compress']:
                    raw = zlib.compress(raw, 6)
                b64 = base64.b64encode(raw).decode('ascii')
                # the base64 payload was computed but never
                # written -- SeeMS/pwiz read ZERO points from these
                # chromatograms (NullReference in
                # seems.Chromatogram.get_Points). <binary> is required
                # by the mzML schema, and the time array unit is UO
                # (UO:0000010 second), not MS:1000004.
                add('<binaryDataArray encodedLength="%d">%s%s'
                    '<cvParam cvRef="MS" accession="%s" name="%s"%s/>'
                    '<binary>%s</binary>'
                    '</binaryDataArray>\n'
                    % (len(b64), comp, prec, acc, nm, unit, b64))
            add('</binaryDataArrayList>\n</chromatogram>\n')

        def _chrom_id(raw_label, prefix=''):
            # human-readable chromatogram ids -- MZmine and
            # SeeMS display the id; 'chromatogram=N' carries no
            # information (every trace also showed as 'TIC' from the
            # cvParam type). Keep [A-Za-z0-9 .+_-] only.
            s = ('%s%s' % (prefix, raw_label)).replace(' ', '_')
            s = ''.join(c if (c.isalnum() or c in '.+-_') else '_'
                        for c in s)
            return s

        def _chrom_header(cid, i, npts):
            # mzML schema -- defaultArrayLength + a
            # chromatogram-type cvParam (readers use the id/userParams
            # for the label; TIC is the generic type for intensity
            # chromatograms -- MZmine's strict parser requires one)
            add('<chromatogram id=%s index="%d" '
                'defaultArrayLength="%d">\n'
                % (quoteattr(cid), i, npts))
            add('<cvParam cvRef="MS" accession="MS:1000235" '
                'name="total ion current chromatogram"/>\n')

        # Display order: BPC first, then
        # the analytical UV trace, UV_VIS_2, then Pump_Pressure,
        # ColumnOven_Temp, Temp_Lamphouse (lamp LAST), then anything
        # else alphabetical. BPC was previously last.
        def _extra_rank(ds):
            if ds.startswith('BPC'):
                return 0
            if ds.startswith('UV_VIS'):
                return 1
            if ds.startswith('Gradient'):
                return 2
            if 'Pump' in ds or 'Pressure' in ds:
                return 3
            if 'ColumnOven' in ds or 'Column' in ds:
                return 4
            if 'Lamphouse' in ds or 'Lamp' in ds:
                return 5
            return 6

        ordered = sorted(extras, key=lambda e: (_extra_rank(e[0]), e[0]))
        bpc = [e for e in ordered if _extra_rank(e[0]) == 0]
        rest = [e for e in ordered if _extra_rank(e[0]) != 0]

        def _write_extra(i, ds, wl, pts):
            chrom_offsets.append((i, sum(len(p) for p in parts)))
            # UV_VIS_2..n traces get the same readable id as
            # UV1 (e.g. UV_254_nm) when the HyStar method holds the
            # channel wavelength; raw UV_VIS_n only as a fallback.
            m = re.match(r'(?i)UV_VIS_\d+$', ds.strip())
            if m and wl is not None:
                cid = _chrom_id('%g nm' % wl, 'UV_')
            else:
                cid = _chrom_id(ds)
            _chrom_header(cid, i, len(pts))
            add('<userParam name="trace" value=%s/>\n' % quoteattr(ds))
            if wl is not None:
                add('<userParam name="wavelength" value="%g nm"/>\n' % wl)
            _chrom_arrays(pts)
            return i + 1

        i = 0
        for ds, wl, pts in bpc:
            i = _write_extra(i, ds, wl, pts)
        for label, pts in uv:
            chrom_offsets.append((i, sum(len(p) for p in parts)))
            _chrom_header(_chrom_id(label, 'UV_'), i, len(pts))
            add('<userParam name="wavelength" value=%s/>\n'
                % quoteattr(label))
            _chrom_arrays(pts)
            i += 1
        for ds, wl, pts in rest:
            i = _write_extra(i, ds, wl, pts)
        add('</chromatogramList>\n')
    add('</run>\n</mzML>\n')
    add('<indexList count="%d"><index name="spectrum">\n'
        % (2 if chrom_offsets else 1))
    indexlist_open = '<indexList count="%d"><index name="spectrum">\n' \
        % (2 if chrom_offsets else 1)
    indexlist_open_len = len(indexlist_open)
    for scan_id, off in offsets:
        add('<offset idRef="scan=%d">%d</offset>\n' % (scan_id, off))
    add('</index>\n')
    if chrom_offsets:
        add('<index name="chromatogram">\n')
        for i, off in chrom_offsets:
            add('<offset idRef="chromatogram=%d">%d</offset>\n' % (i, off))
        add('</index>\n')
    add('</indexList>\n')
    add('<indexListOffset>%d</indexListOffset>\n'
        % (sum(len(p) for p in parts) - indexlist_open_len))
    body = ''.join(parts)
    checksum = hashlib.sha1(body.encode('utf-8')).hexdigest()
    if path.endswith('.gz'):
        f = gzip.open(path, 'wt', encoding='utf-8')
    else:
        f = open(path, 'w', encoding='utf-8')
    with f:
        f.write(body)
        f.write('<fileChecksum>%s</fileChecksum>\n' % checksum)
        f.write('</indexedmzML>\n')

# ================================================================== mzXML writer

def _spectrum_mzxml(sp, opts):
    peaks = sp['peaks']
    n = len(peaks)
    mz = [p[0] for p in peaks]
    it = [p[1] for p in peaks]
    tic = float(sum(it))
    bpm = max(it) if it else 0.0
    bpv = mz[it.index(bpm)] if it and bpm > 0 else 0.0
    end = 'f' if opts['precision'] == 32 else 'd'
    raw = struct.pack('>%d%s' % (n, end), *mz)
    raw += struct.pack('>%d%s' % (n, end), *it)
    if opts['compress']:
        raw = zlib.compress(raw, 6)
    b64 = base64.b64encode(raw).decode('ascii')
    cent = ' centroided="1"'
    o = ['<scan num="%d" msLevel="%d" peaksCount="%d" retentionTime="PT%.6fM" '
         'polarity="%s" lowMz="%.6g" highMz="%.6g" basePeakMz="%.6g" '
         'basePeakIntensity="%.6g" totIonCurrent="%.6g" scanType="Full"%s>\n'
         % (sp['scan_index'], sp['ms_level'], n, sp['rt_sec'] / 60.0,
            '+' if sp['polarity'] == 'positive' else '-',
            mz[0] if mz else 0.0, mz[-1] if mz else 0.0, bpv, bpm, tic,
            cent)]
    if sp['ms_level'] >= 2 and sp.get('precursor_mz') is not None:
        ce = ' collisionEnergy="%.4g"' % sp['ce'] if sp.get('ce') is not None else ''
        o.append('<precursorMz precursorIntensity="%.6g"%s>%.6g</precursorMz>\n'
                 % (sp.get('precursor_intensity') or 0, ce, sp['precursor_mz']))
    o.append('<peaks compressionType="%s" compressedLen="%d" precision="%d" '
             'byteOrder="network" contentType="m/z-intensity">%s</peaks>\n'
             % ('zlib' if opts['compress'] else 'none',
                len(raw), opts['precision'], b64))
    o.append('</scan>\n')
    return ''.join(o)


def write_mzxml(path, spectra, opts, source_name, start_stamp=None):
    parts = []
    add = parts.append
    add('<?xml version="1.0" encoding="utf-8"?>\n')
    add('<mzXML xmlns="http://sashimi.sourceforge.net/schema_revision/mzXML_3-2" '
        'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" '
        'xsi:schemaLocation="http://sashimi.sourceforge.net/schema_revision/'
        'mzXML_3-2 http://sashimi.sourceforge.net/schema_revision/mzXML_3-2/'
        'mzXML_idx_3.2.xsd">\n')
    # mzXML_3.2 has NO standard
    # place for the acquisition start (no startTimeStamp equivalent;
    # mzXML retentionTime is relative to run start). The acquisition
    # date is stamped as an XML comment just before <msRun> opens --
    # comments are tolerated by every mzXML parser and keep the run
    # date traceable by eye, without pretending to the standard.
    if start_stamp:
        add('<!-- d2mz: acquisition start %s -->\n' % start_stamp)
    add('<msRun>\n')
    offsets = []
    for k, sp in enumerate(spectra, 1):
        _progress(k, len(spectra), 'writing mzXML')
        offsets.append((sp['scan_index'], sum(len(p) for p in parts)))
        add(_spectrum_mzxml(sp, opts))
    _progress_end()
    add('</msRun>\n')
    add('<index name="scan">\n')
    for scan_id, off in offsets:
        add('<offset id="%d">%d</offset>\n' % (scan_id, off))
    add('</index>\n')
    add('<indexOffset>%d</indexOffset>\n' % sum(len(p) for p in parts))
    add('<instrument manufacturer="d2mz" model=%s ionisation="ESI" msType="TOF"/>\n'
        % quoteattr(source_name))
    add('<dataProcessing software="d2mz %s" conversionDate="%s" '
        'lockMasses=%s precursorRephasing=%s/>\n'
        % (__version__,
           time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
           quoteattr('built-in Bruker FA/Ac/Na cluster series '
                     '(38 positive / 36 negative ions)'),
           quoteattr('on' if opts['prec_correct'] else 'off')))
    add('</mzXML>\n')
    f = gzip.open(path, 'wt', encoding='utf-8') if path.endswith('.gz') \
        else open(path, 'w', encoding='utf-8')
    with f:
        f.write(''.join(parts))

# ================================================================== citation

CITATION_TEXT = (
    "d2mz was written by Vibe (AI coding agent powered by GLM, developed "
    "by Mistral AI) at the request of, and under the scientific direction "
    "of, Stephane GREFF (IMBE, Aix-Marseille Universite).\n\n"
    "GREFF S. (2026). d2mz {version}: batch Bruker .d to mzML/mzXML "
    "converter with lock-mass recalibration and precursor re-phasing. "
    "https://github.com/stephanegreff/d2mz (GPL-3.0 licence).\n\n"
    "Funding: BeeMed project, ANSES PNR EST 22-EST-063.")


def write_citation(folder):
    """Write CITATION.cff (GitHub-native), CITATION.bib and the plain
    citation string into the annexes/ subfolder of the output ."""
    cff = (
        "cff-version: 1.2.0\n"
        'message: "If you use d2mz to convert or re-phase your data, please cite it as below."\n'
        "type: software\n"
        "title: d2mz\n"
        "authors:\n"
        "  - family-names: GREFF\n"
        "    given-names: Stephane\n"
        "    affiliation: IMBE, Aix-Marseille Universite\n"
        "acknowledgements: >-\n"
        "  The d2mz software was written by Vibe (AI coding agent\n"
        "  powered by GLM, developed by Mistral AI) at the request of,\n"
        "  and under the scientific direction of, the author.\n"
        "  Funding: BeeMed project, ANSES PNR EST 22-EST-063.\n"
        "version: \"%s\"\n"
        "date-released: \"%s\"\n"
        "license: GPL-3.0\n"
        "repository-code: \"https://github.com/stephanegreff/d2mz\"\n"
        "abstract: >\n"
        "  Batch converter for Bruker qTOF .d folders (BAF format) to\n"
        "  mzML 1.1.0 / mzXML 3.2, with lock-mass mass recalibration and\n"
        "  precursor re-phasing onto the recalibrated MS1 apex. Pure\n"
        "  Python, standard library only.\n"
        "keywords:\n"
        "  - mass spectrometry\n"
        "  - Bruker\n"
        "  - mzML\n"
        "  - mzXML\n"
        "  - recalibration\n"
        "  - precursor correction\n"
        % (__version__, __date__))
    bib = (
        "@misc{greff_d2mz_%s,\n"
        "  author       = {GREFF, St{\\'e}phane},\n"
        "  title        = {{d2mz}: batch {Bruker} .d to mzML/mzXML converter "
        "with lock-mass recalibration and precursor re-phasing},\n"
        "  year         = {2026},\n"
        "  version      = {%s},\n"
        "  publisher    = {GitHub},\n"
        "  note         = {Software written by Vibe (AI coding agent\n"
        "                  powered by GLM, developed by Mistral AI).\n"
        "                  Funding: BeeMed project, ANSES PNR EST\n"
        "                  22-EST-063},\n"
        "  url          = {https://github.com/stephanegreff/d2mz},\n"
        "  license      = {GPL-3.0}\n"
        "}\n" % (__date__[:4], __version__))
    with open(os.path.join(folder, 'CITATION.cff'), 'w') as f:
        f.write(cff)
    with open(os.path.join(folder, 'CITATION.bib'), 'w') as f:
        f.write(bib)
    with open(os.path.join(folder, 'CITATION.txt'), 'w') as f:
        f.write(CITATION_TEXT.format(version=__version__) + '\n')


# ================================================================== conversion

def open_baf(d_folder):
    baf = os.path.join(d_folder, 'analysis.baf')
    if not os.path.exists(baf):
        alt = [f for f in sorted(os.listdir(d_folder)) if f.endswith('.baf')]
        if not alt:
            raise FileNotFoundError('no .baf file in %s' % d_folder)
        baf = os.path.join(d_folder, alt[0])
    return baf


def ms1_peaks_iter(baf, sections, cal, spectra_rows, steps):
    """Yield (rt_sec, maxima_peaks) for every MS1 scan."""
    for s in spectra_rows:
        sid, rt = s.get('Id'), s.get('Rt', s.get('RT'))
        if sid is None or rt is None:
            continue
        if any(st.get('Mass') and st['Mass'] > 0 for st in steps.get(sid, [])):
            continue  # MS2
        moff = sections.get(sid, {}).get('maxima')
        if moff is None:
            continue
        peaks = read_maxima(baf, moff, cal)
        if peaks:
            yield float(rt), peaks


def observe_both_series(ms1_iter, opts):
    """Observe BOTH calibrant series (positive and negative clusters) on
    the MS1 scans. Returns the two calibrators (unfitted)."""
    calibs = {}
    for pol, ions in BRUKER_FA_AC_NA_CLUSTERS.items():
        calibs[pol] = LockMassCalibrator([m for m, _f in ions],
                                         opts['lock_ppm'], opts['min_lock'])
    for rt, peaks in ms1_iter:
        for c in calibs.values():
            c.observe_scan(rt, peaks)
    return calibs


def pick_polarity(calibs, opts):
    """Auto-detect the acquisition polarity from the calibrant matches.

    Physically sound: the calibrant's negative cluster ions are not
    transmitted in positive mode (and vice versa), so the correct series
    matches many MS1 scans and the wrong one almost none. The decision
    is a match COUNT, never a fit quality. --polarity overrides it.
    """
    override = opts.get('polarity')
    if override in ('positive', 'negative'):
        return override
    n_pos = len(calibs['positive'].points)
    n_neg = len(calibs['negative'].points)
    return 'positive' if n_pos >= n_neg else 'negative'


def find_method_folder(d_folder):
    """The acquisition method folder (*.m) inside a .d folder: it holds
    the otoControl method files (e.g. microTOFQImpactAcquisition.method,
    hystar.method) written at acquisition time."""
    base_depth = d_folder.rstrip(os.sep).count(os.sep)
    for root, dirs, _files in os.walk(d_folder):
        if root.rstrip(os.sep).count(os.sep) - base_depth > 4:
            dirs[:] = []
            continue
        for dn in dirs:
            if not dn.endswith('.m'):
                continue
            path = os.path.join(root, dn)
            try:
                if any(f.endswith('.method') for f in os.listdir(path)):
                    return path
            except OSError:
                pass
    return None


def method_reference_masses(d_folder):
    """EXACT calibrant reference list DECLARED by the acquisition method
    (.m): the vector <para_vec_double permname=
    "Calibration_LastCalibrationReferenceMass"> (with its
    ...MassName twin) inside the acquisition .method file lists the
    reference masses the instrument itself calibrates on (e.g. 37
    Na(HCOONa)n / Na(CH3COONa)n / Na(NaAc)x(NaF)y clusters on the lab's
    Impact II POS runs). The method KNOWS
    its calibrant -- read the exact list and fit it directly instead of
    trying every adduct family. Returns a list of masses (>= 4) or
    None when the method declares no usable vector."""
    mf = find_method_folder(d_folder)
    if not mf:
        return None
    for fn in sorted(os.listdir(mf)):
        if not fn.endswith('.method') or 'acquisition' not in fn.lower():
            continue
        txt = None
        for enc in ('utf-8', 'latin-1'):
            try:
                with open(os.path.join(mf, fn), encoding=enc) as f:
                    txt = f.read()
                break
            except (UnicodeDecodeError, OSError):
                continue
        if txt is None or 'Calibration_LastCalibrationReferenceMass"' \
                not in txt:
            continue
        i = txt.find('permname='
                     '"Calibration_LastCalibrationReferenceMass"')
        if i < 0:
            continue
        end = txt.find('</para_vec', i)
        chunk = txt[i:end] if end > 0 else ''
        masses = [float(v) for v in re.findall(
            r'entry_double value="([0-9.eE+-]+)"', chunk)]
        if len(masses) >= 4:
            return masses
    return None


def acquisition_start_stamp(d_folder, log=print):
    """Acquisition start time for the mzML <run startTimeStamp>.
    MZmine shows run_date as 1970-01-01T00:00 when the mzML run
    element carries the epoch fallback, i.e. when no acquisition
    date is read. PRIMARY
    source: SampleInfo.xml <AnalysisHeader CreationDateTime="...">
    (written by HyStar at acquisition start, ISO-8601 with local
    UTC offset, UTF-16 or UTF-8 depending on otoControl version).
    Fallback: mtime of the .baf (file written during the run,
    local time -> UTC). Returns 'YYYY-MM-DDTHH:MM:SSZ' | None."""
    import datetime as _dt
    sinfo = os.path.join(d_folder, 'SampleInfo.xml')
    cand = []
    if os.path.exists(sinfo):
        try:
            with open(sinfo, 'rb') as f:
                raw = f.read()
            try:
                raw.decode('utf-8'); txt = raw.decode('utf-8')
            except UnicodeDecodeError:
                txt = raw.decode('utf-16')
            m = re.search(r'CreationDateTime="([^"]+)"', txt)
            if m:
                cand.append(('SampleInfo.xml CreationDateTime', m.group(1)))
        except OSError as e:
            log('    (SampleInfo.xml unreadable: %s)' % e)
    stamp = None
    for src, s in cand:
        s = s.strip()
        if not s:
            continue
        try:
            t = _dt.datetime.fromisoformat(s)
        except ValueError:
            continue
        if t.tzinfo is None:          # no offset: assume local time
            t = t.astimezone()
        t = t.astimezone(_dt.timezone.utc)
        stamp = t.strftime('%Y-%m-%dT%H:%M:%SZ')
        break
    if stamp is None:
        baf = os.path.join(d_folder, 'analysis.baf')
        if not os.path.exists(baf):
            baf = ([os.path.join(d_folder, f)
                    for f in sorted(os.listdir(d_folder))
                    if f.endswith('.baf')] + [None])[0]
        if baf:
            stamp = time.strftime('%Y-%m-%dT%H:%M:%SZ',
                                  time.gmtime(os.path.getmtime(baf)))
            log('    run start (fallback, baf mtime UTC): %s' % stamp)
    else:
        log('    run start: %s' % stamp)
    return stamp


def method_polarity(d_folder):
    """Acquisition polarity DECLARED in the .m method folder (upstream
    check). PRIMARY source (2026-11-07_z): otoControl writes the
    acquisition polarity as <para_int value="N" permname=
    "Mode_IonPolarity"/> inside the acquisition .method file
    (microTOFQImpactAcquisition.method): 0=positive, 1=negative --
    validated 27/27 on the lab's real .m folders. The old word-count
    heuristic (kept as fallback only) ALWAYS returned None on Impact II
    methods because the acquisition file mentions both polarities
    (calibrant lists, switching settings: 61 'positive' + 58
    'negative' on a pure POS run). Returns 'positive' | 'negative' |
    None (None = unknown or switching)."""
    mf = find_method_folder(d_folder)
    if not mf:
        return None
    vals = set()
    for fn in sorted(os.listdir(mf)):
        if not fn.endswith('.method') or 'acquisition' not in fn.lower():
            continue
        txt = None
        for enc in ('utf-8', 'latin-1'):
            try:
                with open(os.path.join(mf, fn), encoding=enc) as f:
                    txt = f.read()
                break
            except (UnicodeDecodeError, OSError):
                continue
        if txt is None:
            continue
        for m in re.finditer(r'<para_int\s+value="(\d+)"\s+'
                             r'permname="Mode_IonPolarity"', txt):
            vals.add(m.group(1))
    if vals == {'0'}:
        return 'positive'
    if vals == {'1'}:
        return 'negative'
    if vals:
        return None  # mixed values: polarity switching or unknown code
    npos = nneg = 0
    for fn in sorted(os.listdir(mf)):
        if not fn.endswith('.method'):
            continue
        txt = None
        for enc in ('utf-8', 'latin-1'):
            try:
                with open(os.path.join(mf, fn), encoding=enc) as f:
                    txt = f.read()
                break
            except (UnicodeDecodeError, OSError):
                continue
        if txt is None:
            continue
        low = txt.lower()
        if 'polarity' not in low and 'ionmode' not in low:
            continue
        npos += low.count('positive')
        nneg += low.count('negative')
    if npos and not nneg:
        return 'positive'
    if nneg and not npos:
        return 'negative'
    return None


def fit_lock_model(d_folder, opts, log=print):
    """Light pre-pass: decode MS1 scans only, detect the polarity and fit
    the lock-mass model. Returns (polarity, model, cal, calib, reason) —
    model is None when the calibrant is too weak (below the minimum lock
    observations), reason is the CalibrantAbsentError message of a
    deferred file, None otherwise. Used to decide, before conversion, which file can
    lend its calibration to a same-polarity neighbour; cal and the
    fitted calib are REUSED by the conversion pass (each file is
    calibrated once, never twice).
    """
    baf = open_baf(d_folder)
    _, sections = read_index(os.path.join(d_folder, 'analysis.baf_idx'))
    cal = read_calibration(os.path.join(d_folder, 'analysis.baf_xtr'))
    spectra_rows, steps, _ = load_spectra_table(d_folder, log)
    if cal.get('mode') not in (4, 5):
        try:
            cal = fit_impact_calibration(
                baf, sections, spectra_rows, steps, opts,
                method_polarity(d_folder), log,
                method_reference_masses(d_folder))
        except CalibrantAbsentError as e:
            log('    calibrant absent (%s)' % e)
            log('    -> file DEFERRED: its calibration will be borrowed '
                'at the end from a same-mode (same polarity) file of '
                'the same folder')
            return 'deferred', None, None, None, str(e)
    if cal.get('segment_fit'):
        # The calibrant lives in segment 2
        # ONLY -- there are no calibrant ions in the analysis segment, so
        # NO lock-mass observation over the run (it was scanning every
        # MS1 scan for nothing). The single quadratic fitted on the
        # calibrant segment is applied unchanged to the whole file, the
        # same way as the reference extraction script.
        cal.pop('_impact_scans', None)
        m_pol = method_polarity(d_folder)
        polarity = m_pol or cal.get('polarity')
        if m_pol and cal.get('polarity') and m_pol != cal['polarity']:
            log('    WARNING: ladder polarity (%s) differs from the .m '
                'method folder (%s)' % (cal['polarity'], m_pol))
        calib = LockMassCalibrator([], opts['lock_ppm'], opts['min_lock'])
        calib.polarity = polarity
        calib.method_polarity = m_pol
        calib.segment_only = True
        calib.n_lock_pos = calib.n_lock_neg = 0
        log('    polarity: %s%s' % (
            polarity,
            ' (declared in .m method folder)' if m_pol == polarity
            else ' (calibrant ladder)'))
        log('    no lock-mass observation (calibrant in segment 2 '
            'only): single quadratic applied to the whole run')
        return polarity, None, cal, calib, None
    recs = cal.pop('_impact_scans', None)
    if recs is not None:
        # impact mode: observe the lock masses on the SAME raw records
        # that were read once for the quad fit (no second pass)
        def _ms1_from_recs():
            for rt, _xs, pairs, _c in recs:
                yield rt, [(index_to_mz(cal, i), v)
                           for i, v in pairs]
        calibs = observe_both_series(_ms1_from_recs(), opts)
    else:
        calibs = observe_both_series(
            ms1_peaks_iter(baf, sections, cal, spectra_rows, steps), opts)
    polarity = pick_polarity(calibs, opts)
    calib = calibs[polarity]
    calib.fit()
    calib.polarity = polarity
    calib.method_polarity = method_polarity(d_folder)
    calib.n_lock_pos = len(calibs['positive'].points)
    calib.n_lock_neg = len(calibs['negative'].points)
    # (polarity, model, cal, fitted-calib) -- the caller reuses cal and
    # calib during the conversion pass instead of fitting them twice
    return polarity, calib.model, cal, calib, None


def convert_one(d_folder, opts, log=print, borrow_model=None,
                borrow_from=None, prefit=None):
    name = os.path.basename(os.path.normpath(d_folder))
    log('  reading %s' % name)
    method_pol = method_polarity(d_folder)
    # CE sweep mode (ramp vs constant) from the .method file --
    # feeds the per-MS2-scan d2mz:collision_energy userParam
    opts['ce_info'] = method_ce_sweep(d_folder)
    if method_pol:
        log('    acquisition polarity in .m method folder: %s' % method_pol)
    baf = os.path.join(d_folder, 'analysis.baf')
    if not os.path.exists(baf):
        alt = [f for f in sorted(os.listdir(d_folder)) if f.endswith('.baf')]
        if not alt:
            raise FileNotFoundError('no .baf file in %s' % d_folder)
        baf = os.path.join(d_folder, alt[0])
    _, sections = read_index(os.path.join(d_folder, 'analysis.baf_idx'))
    cal = read_calibration(os.path.join(d_folder, 'analysis.baf_xtr'))
    spectra_rows, steps, psv = load_spectra_table(d_folder, log)
    # acquisition-mode COVERAGE: d2mz is
    # validated on DDA CID with constant CE or CE ramps (batch
    # 2026-11-06, A1/B1 Impact II). Anything else gets a WARNING,
    # never a silent conversion.
    _modes = psv.get('_modes', {}) if isinstance(psv, dict) else {}
    _cei = opts.get('ce_info') or {}
    _acq_mode = 'DDA, CE %s' % ('ramp/sweep' if _cei.get('sweep')
                               else 'constant')
    _mode_notes = []
    if _modes.get('iscid'):
        _mode_notes.append('IsCID active on %d scan(s) '
                           '(Transfer_ISCIDEnergy_Act > 0)'
                           % _modes['iscid'])
    if _modes.get('fragmode'):
        _mode_notes.append('non-CID fragmentation mode value(s) %s '
                           '(MSMS_FragmentationMode_Act)'
                           % _modes['fragmode'])
    if _cei.get('manual_list_ce'):
        _mode_notes.append('targeted MS/MS list (MSMSManual_List*)')
    # v_co: for the residual-list note the console warning is
    # printed ONCE, AFTER the acquired precursors have confirmed or
    # disproved it (see the retraction after the scan loop): a pure
    # auto-DDA file whose .method carries a leftover
    # MSMSManual_List* shows the calm withdrawal note instead of a
    # scary !! banner -- no systematic uncertainty on the
    # conversion (request 25/11/2026).
    opts['acquisition_mode'] = _acq_mode          # mzML userParam
    opts['mode_warning'] = '; '.join(_mode_notes) or None
    if prefit is not None and prefit[0] is not None:
        cal = prefit[0]     # quad fitted once in the pre-pass: reused
    elif cal.get('mode') not in (4, 5):
        cal = fit_impact_calibration(baf, sections, spectra_rows, steps,
                                     opts, method_pol, log,
                                     method_targets=method_reference_masses(
                                         d_folder))
    cal.pop('_impact_scans', None)   # raw records freed after the pre-pass
    if cal.get('quad') is not None:
        log('    %d spectra rows (impact quadratic calibration)' %
            len(spectra_rows))
    else:
        log('    %d spectra rows (calibration mode %d)'
            % (len(spectra_rows), cal['mode']))
    uv_wanted = opts['uv_traces']
    uv_traces, uv_apex = {}, []
    extra_chroms = []          # other existing traces -> mzML too
    named_traces = read_all_trace_chunks(d_folder)
    hystar_wls = hystar_uv_wavelengths(d_folder)
    uv_srcs = [(ds, pts) for ds, _ins, pts in named_traces
               if re.search(r'(?i)uv', ds)]
    uv1_used = None            # the UV1 points used as the main trace
    if not uv_wanted:
        # Record the FIRST acquired UV
        # chromatogram (UV1) of the file -- NOT a hardcoded 254 nm
        # (the lab's UV1 is e.g. 280 nm)
        # The analytical UV channel
        # is the FIRST acquired one ABOVE 250 nm -- UV1 when it runs at
        # 254/280 nm, but UV2 (e.g. 254 nm) when UV1 sits at 210 nm.
        # Falls back to the first channel when no hystar wavelength is
        # known or none is above 250 nm.
        if uv_srcs:
            pick = None
            for desc, pts in uv_srcs:
                m = re.match(r'(?i)UV_VIS_(\d+)', desc.strip())
                if not m:
                    continue
                wl = hystar_wls.get(int(m.group(1)))
                if wl is not None and wl > 250.0:
                    pick = (desc, pts, wl)
                    break
            if pick is None:
                desc, pts = uv_srcs[0]
                m = re.match(r'(?i)UV_VIS_(\d+)', desc.strip())
                wl = hystar_wls.get(int(m.group(1))) if m else None
                pick = (desc, pts, wl)
            desc, uv1_pts, wl1 = pick
            label = ('%g nm' % wl1) if wl1 else 'UV1'
            log('    UV trace: first acquired channel above 250 nm '
                '(%s)%s' % (desc, ' at %g nm' % wl1 if wl1 else ''))
            uv_wanted = [label]
            key = wl1 if wl1 else label    # numeric wl key when known
            uv_scans = [(rt, {key: inten}) for rt, inten in uv1_pts]
            uv_db = None
            uv_traces = {label: uv1_pts}
            uv1_used = uv1_pts
        else:
            # no named trace: DAD cube / legacy tables, first wavelength
            uv_wanted = uv_first_wavelength(d_folder, spectra_rows)
            if uv_wanted:
                log('    UV trace: first acquired channel (UV1) at %g nm'
                    % uv_wanted[0])
            uv_scans, uv_db = read_uv_scans(d_folder, spectra_rows,
                                             uv_wanted)
            if uv_scans:
                uv_traces = uv_traces_from_scans(uv_scans, uv_wanted)
    else:
        uv_scans, uv_db = read_uv_scans(d_folder, spectra_rows, uv_wanted)
        if uv_scans:
            uv_traces = uv_traces_from_scans(uv_scans, uv_wanted)
    # Every OTHER existing trace goes into the mzML as well (cheap in
    # space): the second UV channel
    # (with its hystar wavelength), lamp/column temperatures, pump
    # pressure, MS base peak... The UV1 trace is skipped (already the
    # analytical chromatogram above)
    _skip = set()
    if uv1_used is not None:
        _skip.add(id(uv1_used))
    for ds, _ins, pts in named_traces:
        if id(pts) in _skip:
            continue
        m = re.match(r'(?i)UV_VIS_(\d+)', ds)
        wl = hystar_wls.get(int(m.group(1))) if m else None
        if wl is not None and uv_wanted and any(
                abs(wl - w) < 1e-6 for w in uv_wanted
                if isinstance(w, float)):
            continue        # duplicates a requested --uv-trace channel
        extra_chroms.append((ds, wl, pts))
    if extra_chroms:
        log('    extra traces -> mzML: %s' % ', '.join(
            '%s (%d pts)' % (ds, len(pts))
            for ds, _wl, pts in extra_chroms))
    # The PROGRAMMED gradient from the
    # HyStar method replaces the lamp temperature trace (the %B
    # channel is not recorded as data by the instrument).
    lc = hystar_lc_method(d_folder)
    if lc:
        extra_chroms = [e for e in extra_chroms
                        if 'lamp' not in e[0].lower()]
        gpts = lc_gradient_points(lc)
        if gpts:
            extra_chroms.append(('Gradient_percent_B', None, gpts))
        log('    LC program (hystar.method): solvents %s; flow %s; '
            'gradient %s; oven %s'
            % ('; '.join('%s=%s' % kv for kv in sorted(
                    lc['solvents'].items())) or 'n/a',
               ', '.join('%g min: %g ml/min' % fv
                         for fv in lc['flow']) or 'n/a',
               lc_gradient_str(lc),
               ('%.1f C' % lc['oven_c']) if lc['oven_c'] else 'n/a'))
    if uv_scans:
        apex_rts = []
        for label, pts in sorted(uv_traces.items()):
            for rt_a, inten_a in uv_apices(pts):
                uv_apex.append(dict(trace=label, rt_sec=rt_a,
                                    intensity=inten_a))
                apex_rts.append(rt_a)
        apex_rts = sorted(set(apex_rts))
        if uv_db:
            # cube source: traces were read lambda-targeted; pull the
            # full wavelength spectrum only for the few apex scans
            uv_apex_spectra_ = _read_dad_full_at(uv_db, apex_rts)
        else:
            uv_apex_spectra_ = uv_apex_spectra(uv_scans, apex_rts)
        for spec in uv_apex:
            rt_s, wl_int = min(uv_apex_spectra_,
                               key=lambda sc: abs(sc[0] - spec['rt_sec']))
            spec['spectrum'] = wl_int
        log('    UV traces: %s' % ', '.join(
            '%s (%d pts)' % (l, len(p)) for l, p in sorted(uv_traces.items())))
        log('    UV apex spectra: %d (at %s)' % (
            len(uv_apex),
            ', '.join('%.2f min' % (r / 60.0) for r in apex_rts) or 'none'))
    else:
        log('    UV traces: none found')

    # pass 1: decode
    raw_specs = []
    n_rows = len(spectra_rows)
    for k, s in enumerate(spectra_rows, 1):
        _progress(k, n_rows, 'decoding', step=500)
        sid, rt = s.get('Id'), s.get('Rt', s.get('RT'))
        if sid is None or rt is None:
            continue
        sec = sections.get(sid, {})
        sts = steps.get(sid, [])
        parent = None
        for st in sts:
            m = st.get('Mass')
            if m and m > 0:
                parent = m
                break
        ms_level = 2 if parent is not None else 1
        peaks = []
        moff, soff = sec.get('maxima'), sec.get('spectrum')
        if moff is not None:
            peaks = read_maxima(baf, moff, cal)
        raw_specs.append(dict(
            id=sid, rt=float(rt), ms_level=ms_level, peaks=peaks,
            declared_parent=parent,
            ce=psv.get(sid, {}).get(5), iso=psv.get(sid, {}).get(7)))
    raw_specs.sort(key=lambda r: r['rt'])
    _progress_end()

    # pass 2: polarity auto-detection + lock-mass calibration
    if prefit is not None and prefit[1] is not None:
        calib = prefit[1]       # fitted once in the pre-pass: reused
        polarity = calib.polarity
        if getattr(calib, 'segment_only', False):
            log('    polarity: %s (calibrant-segment calibration%s)'
                % (polarity,
                   ', borrowed from %s' % calib.calibration_borrowed_from
                   if getattr(calib, 'calibration_borrowed_from', None)
                   else ''))
        if method_pol and polarity != method_pol:
            log('    WARNING: detected polarity (%s) differs from the '
                'acquisition method folder (.m: %s)'
                % (polarity, method_pol))
    elif cal.get('segment_fit'):
        # segment-calibrated impact file reached without the pre-pass
        # (e.g. the parabens anchor measurement): same rule -- no
        # lock-mass observation over the analysis segment
        polarity = method_pol or cal.get('polarity')
        calib = LockMassCalibrator([], opts['lock_ppm'], opts['min_lock'])
        calib.polarity = polarity
        calib.method_polarity = method_pol
        calib.segment_only = True
        calib.n_lock_pos = calib.n_lock_neg = 0
        log('    polarity: %s (calibrant-segment calibration)'
            % polarity)
    else:
        calibs = observe_both_series(
            ((r['rt'], r['peaks']) for r in raw_specs
             if r['ms_level'] == 1 and r['peaks']), opts)
        polarity = pick_polarity(calibs, opts)
        calib = calibs[polarity]
        other = calibs['negative' if polarity == 'positive' else 'positive']
        log('    polarity detected: %s (calibrant matches %d vs %d)'
            % (polarity, len(calib.points), len(other.points)))
        if method_pol and polarity != method_pol:
            log('    WARNING: detected polarity (%s) differs from the '
                'acquisition method folder (.m: %s)' % (polarity, method_pol))
        calib.fit()
        calib.polarity = polarity
        calib.method_polarity = method_pol
        calib.n_lock_pos = len(calibs['positive'].points)
        calib.n_lock_neg = len(calibs['negative'].points)
    # 2026-11-25_cp: durable per-file calibration trace -- the fit
    # stats (n_targets, residual ppm) used to reach the console log
    # only; the report needs them (port of the 2026-11-16_ci surface
    # onto the GitHub _co lineage; deferred-and-borrowed files carry
    # the donor's fit plus calibration_borrowed_from)
    _cfit = cal.get('fit') or {}
    _cborrowed = getattr(calib, 'calibration_borrowed_from', None) \
        or getattr(calib, 'borrowed_from', None)
    calib.calibration_trace = dict(
        mode=cal.get('mode'),
        applied=bool(_cfit or cal.get('quad') is not None
                     or cal.get('mode') in (4, 5) or _cborrowed),
        n_targets=_cfit.get('n_targets'),
        ppm_median=_cfit.get('med_ppm'),
        ppm_max=_cfit.get('max_ppm'))

    if getattr(calib, 'segment_only', False):
        # single quadratic from the calibrant segment: no lock-mass
        # observation, nothing to report as lock RMS
        log('    single quadratic applied to the whole run '
            '(calibrant-segment fit%s)'
            % (', borrowed from %s' % calib.calibration_borrowed_from
               if getattr(calib, 'calibration_borrowed_from', None)
               else ''))
    else:
        log('    lock masses: %d MS1 scans matched; ppm RMS %.2f -> %.2f'
            % (len(calib.points), calib.rms_before() or 0.0,
               calib.rms_after() or 0.0))
        if calib.model is None:
            log('    WARNING: only %d lock observation(s) (< %d): '
                'calibrant weak or absent in this run'
                % (len(calib.points), opts['min_lock']))
            if borrow_model is not None:
                calib.model = borrow_model
                calib.borrowed_from = borrow_from
                log('    -> lock-mass model borrowed from adjacent file %s '
                    '(same series/polarity)' % borrow_from)
            else:
                log('    WARNING: no adjacent calibrated file either: '
                    'spectra NOT recalibrated')

    # pass 3: apply + rephase precursors
    out_specs, prec_table = [], []
    last_ms1 = []
    scan_index = 0
    n_raw = len(raw_specs)
    for ridx, r in enumerate(raw_specs, 1):
        _progress(ridx, n_raw, 'recalibrating', step=200)
        if not r['peaks']:
            continue
        scan_index += 1
        peaks = calib.apply(r['rt'], r['peaks'])
        prec_mz = prec_int = None
        if r['ms_level'] == 1:
            last_ms1 = peaks
        elif r['declared_parent'] is not None:
            prec_mz = r['declared_parent']
            if opts['prec_correct']:
                new_mz, inten, dppm = reassign_precursor(
                    r['declared_parent'], last_ms1, opts['prec_window'])
                row = dict(scan=scan_index, rt_min=r['rt'] / 60.0,
                           declared=r['declared_parent'], corrected=new_mz,
                           delta_ppm=dppm, ms1_intensity=inten,
                           ce_window=ce_window(r['ce']))
                prec_table.append(row)
                if new_mz is not None:
                    prec_mz, prec_int = new_mz, inten
        out_specs.append(dict(
            scan_index=scan_index, ms_level=r['ms_level'], rt_sec=r['rt'],
            polarity=polarity, peaks=peaks,
            precursor_mz=prec_mz, precursor_intensity=prec_int,
            precursor_scan=None, ce=r['ce'], isolation_mz=r['iso'],
            # Explicit ramp vs constant CE label: a run can be all-
            # constant or ramped, the mzML must tell which
            ce_label=ce_label(r['ce'], opts.get('ce_info'))))
    _progress_end()
    # the .method OFTEN carries a residual MSMSManual_List* even
    # when the run is pure auto-DDA (batch 261106: every file declared
    # one, only the Parabens runs actually used it). A REALLY targeted
    # run fragments the same handful of precursor masses over and over;
    # an auto-DDA run spreads its MS2 over many distinct precursors.
    # Confirm the manual-list note with the acquired data before
    # warning, so the quality gate stays meaningful.
    # v_co: the retraction used to be gated on `not _modes` -- an
    # empty dict on the analysis.sqlite path (retraction worked), but
    # NEVER empty on the baf-metadata fallback of never-opened .d
    # files, where _modes['iso_act'] counts EVERY MS2 scan (the
    # per-scan MSMS_Isolation_Act flag, pure bookkeeping that
    # generates NO warning note). The gate silently disabled the
    # retraction there (real case 38087: never opened, auto-DDA with
    # 2839 jittering precursors, still flagged 'targeted MS/MS
    # list'). `len(_mode_notes) == 1` already guarantees that no
    # IsCID / fragmode note was added, so nothing else may block it.
    if _mode_notes and len(_mode_notes) == 1 \
            and _cei.get('manual_list_ce'):
        # a REAL manual list replays the exact same precursor m/z
        # every cycle (bit-identical declared value, thousands of
        # repeats); auto-DDA recomputes the precursor from each MS1,
        # so the value jitters cycle to cycle (batch 261106: A1
        # cycles at 479.12032 / 479.11965 / 479.12031). Exact-value
        # repetition is therefore the discriminator.
        _targeted = manual_list_targeted(prec_table)
        if not _targeted:
            log('    note: residual MSMSManual_List* in the .method but '
                'precursor m/z vary cycle to cycle (auto-DDA): '
                'coverage warning withdrawn')
            _mode_notes = []
            opts['mode_warning'] = None
    # v_co: console warning printed HERE, once the retraction has
    # decided -- never before (no systematic !! banner on auto-DDA
    # files with a leftover manual list in the .method).
    if _mode_notes:
        log('    !! acquisition outside the validated scope '
            '(validated: DDA CID constant-CE / CE-ramp): %s'
            % '; '.join(_mode_notes))
    return out_specs, prec_table, calib, dict(
        traces=uv_traces, apexes=uv_apex, extras=extra_chroms,
        lc=lc)


def _d_decodable(full):
    """A .d folder d2mz can decode: either DataAnalysis already
    wrote analysis.sqlite, or the folder was never opened and the raw
    analysis.baf (+ its index) is enough for the baf-metadata fallback
    of load_spectra_table."""
    return (os.path.exists(os.path.join(full, 'analysis.sqlite'))
            or (os.path.exists(os.path.join(full, 'analysis.baf'))
                and os.path.exists(os.path.join(full,
                                                'analysis.baf_idx'))))


def find_d_folders(input_path, log=None):
    """Scan an input folder for convertible Bruker .d acquisitions.

    Returns (d_folders, incomplete): d_folders = the decodable .d
    folders (analysis.sqlite present, OR never opened in DataAnalysis
    with the raw analysis.baf + baf_idx -- the fallback);
    incomplete = the *.d folders with NEITHER (interrupted
    acquisition, partial copy or a folder that only looks like a .d)
    -- the incomplete list is written into the report scan section so
    the user can see exactly why the [i/N] counter of the log is
    smaller than the number of .d folders in the Explorer.
    """
    if not os.path.isdir(input_path):
        raise NotADirectoryError(input_path)
    if input_path.endswith('.d') and _d_decodable(input_path):
        return [input_path], []
    out = []
    incomplete = []
    for root, dirs, files in os.walk(input_path):
        for d in dirs:
            if d.endswith('.d'):
                full = os.path.join(root, d)
                if _d_decodable(full):
                    out.append(full)
                else:
                    incomplete.append(full)
        dirs[:] = [x for x in dirs if not x.startswith('.')]
    if log and incomplete:
        for p in incomplete:
            log('  !! %s: no analysis.sqlite and no analysis.baf -- '
                'not a complete Bruker acquisition (interrupted or '
                'partial copy?), ignored'
                % os.path.basename(os.path.normpath(p)))
    return sorted(out), sorted(incomplete)


def _progress(k, n, label, step=200, width=24):
    """One-line ASCII progress bar, redrawn in place with \\r."""
    if n <= 0 or (k % step and k != n):
        return
    frac = k / float(n)
    filled = int(width * frac)
    sys.stdout.write('\r    %s [%s%s] %3d%% (%d/%d)'
                     % (label, '#' * filled, '.' * (width - filled),
                        int(100 * frac), k, n))
    sys.stdout.flush()


def _progress_end():
    sys.stdout.write('\n')
    sys.stdout.flush()


def banner():
    return r"""
  ================================================================
     ____    ____                       _   _   _
    |  _ \  |  _ \    ___     ___    __| | | \ | |
    | | | | | |_) |  / __|   / _ \  / _` | |  \| |
    | |_| | |  _ <  | (__   | (_) || (_| | | |\  |
    |____/  |_| \_\  \___|   \___/  \__,_| |_| \_|
     Bruker .d -> mzML / mzXML converter (pure stdlib)
     lock-mass recalibration - precursor re-phasing - UV traces
     version %s (%s) - ANSES PNR EST 22-EST-063 (BeeMed)
  ================================================================
""" % (__version__, __date__)


def _file_warnings(report):
    """Per-file quality flags mined from the report. 2026-11-25_cp:
    factored out of main() so every per-series report dump carries
    the same warnings as the end-of-run quality gate (new flag:
    a file whose calibration trace says applied=False)."""
    warns = []
    for _f in report['files']:
        _nm = _f.get('name') or '?'
        if _f.get('error'):
            continue                    # already in problems
        if _f.get('ms2_warning'):
            warns.append((_nm, _f['ms2_warning']))
        if _f.get('mode_warning'):
            warns.append((_nm, 'acquisition outside the validated '
                              'scope (%s): %s -- conversion NOT '
                              'validated for this mode'
                          % (_f.get('acquisition_mode') or '?',
                             _f['mode_warning'])))
        if _f.get('calibration_borrowed_from'):
            warns.append((_nm, 'lock calibration BORROWED from %s '
                              '(no own calibrant observation)'
                          % _f['calibration_borrowed_from']))
        if _f.get('calibration_applied') is False:
            warns.append((_nm, 'NO mass calibration applied '
                              '(no own calibrant fit and no borrowed '
                              'model): m/z NOT recalibrated'))
        _dp, _mp = _f.get('detected_polarity'), \
            _f.get('method_polarity')
        if _dp and _mp and _dp != _mp:
            warns.append((_nm, 'polarity mismatch: detected %s, '
                              'method %s' % (_dp, _mp)))
    return warns


def _calibration_summary(files):
    """Batch-level calibration block (2026-11-25_cp, port of the
    2026-11-16_ci surface): every file's calibration trace with the
    fit statistics (n_targets, residual ppm) -- the GO/NO-GO gate
    before any deposit reads ok / n_uncalibrated / files_above_5ppm."""
    rows = [_f for _f in files if _f.get('error') is None
            and _f.get('calibration_applied') is not None]
    medians = [_f['calibration_ppm_median'] for _f in rows
               if _f.get('calibration_ppm_median') is not None]
    n_own = sum(1 for _f in rows if _f['calibration_applied']
                and not _f.get('calibration_borrowed_from'))
    n_borrowed = sum(1 for _f in rows
                     if _f.get('calibration_borrowed_from'))
    n_uncal = sum(1 for _f in rows if not _f['calibration_applied'])
    maxes = [_f['calibration_ppm_max'] for _f in rows
             if _f.get('calibration_ppm_max') is not None]
    above5 = [_f['name'] for _f in rows
              if (_f.get('calibration_ppm_median') or 0) > 5.0]
    max_med = max(medians) if medians else None
    max_max = max(maxes) if maxes else None
    ok = (n_uncal == 0) and (max_med is None or max_med <= 5.0)
    if n_own and max_med is not None:
        summary = ('%d file(s) with own calibration fit '
                  '(worst median residual %.2f ppm, worst max %.2f ppm)'
                  % (n_own, max_med, max_max))
    elif n_own:
        summary = ('%d file(s) calibrated (instrument calibration, '
                   'no per-file fit statistics)' % n_own)
    else:
        summary = 'NO calibration trace found'
    if n_borrowed:
        summary += '; %d borrowed' % n_borrowed
    if n_uncal:
        summary += '; %d NOT calibrated' % n_uncal
    return dict(ok=ok, n_files_with_trace=len(rows), n_own_fit=n_own,
                n_borrowed=n_borrowed, n_uncalibrated=n_uncal,
                ppm_median_max=max_med, ppm_max_max=max_max,
                files_above_5ppm=above5, summary=summary)


def main(argv=None):
    ap = argparse.ArgumentParser(
        description='Batch Bruker .d -> mzML/mzXML converter with lock-mass '
                    'recalibration and precursor re-phasing (pure stdlib). '
                    'Accepts one or several input folders; each gets its own '
                    'd2mz_output_<source-name>/ folder created next to '
                    'it, with the mzML and mzXML files written in its '
                    'analyses/ subfolder and the precursor tables in '
                    'its annexes/TSV/ subfolder, while the series '
                    'report_<source>.json and the CITATION files sit '
                    'in its annexes/ subfolder.')
    ap.add_argument('--version', action='version',
                    version='d2mz %s (pure-Python Bruker .d to '
                    'mzML/mzXML converter)' % __version__)
    ap.add_argument('inputs', nargs='*',
                    help='folder(s) containing .d data, or single .d '
                         'folder(s), or .zip archive(s) of them '
                         '(outputs then go next to the archive)')
    ap.add_argument('--out', default=None,
                    help='explicit output folder: relative -> created '
                         'INSIDE each input folder, absolute -> exact '
                         'path shared by the whole queue. Default '
                         '(omitted): each input folder gets its own '
                         'd2mz_output_<source-name>/ folder created '
                         'NEXT TO it')
    ap.add_argument('--format', default='both', choices=['mzml', 'mzxml', 'both'])
    ap.add_argument('--precision', type=int, default=32, choices=[32, 64],
                   help='binary array float precision (default: 32)')
    ap.add_argument('--no-compress', action='store_true',
                    help='disable zlib compression of binary arrays')
    ap.add_argument('--gzip', action='store_true',
                    help='write .mzML/.mzXML files gzipped (lossless, '
                         'directly readable by MZmine, OpenMS and '
                         'msconvert; no effect on report.json or the '
                         'precursor tables)')
    ap.add_argument('--zip-output', action='store_true',
                    help='after conversion, also write ONE .zip archive '
                         'next to the output folder '
                         '(d2mz_output_<source>.zip), structured for a '
                         'DIRECT Zenodo upload: '
                         'analyses/ = every .mzML/.mzML.gz/'
                         '.mzXML/.mzXML.gz FLAT, annexes/ = TSV/ + the '
                         'series report_<source>.json + the CITATION '
                         'files. The flat files stay in place: MZmine '
                         'reads the .mzML/.mzML.gz files directly and '
                         'cannot open a .zip')
    ap.add_argument('--uv-trace', action='append', type=float,
                    metavar='NM',
                    help='wavelength (nm) of a UV chromatographic trace '
                         '(repeatable; default: the FIRST acquired UV '
                         'channel (UV1) of each file, e.g. 280 nm)')
    ap.add_argument('--parabene', action='append', default=[],
                    metavar='NAME=RT',
                    help='retention-time standard (repeatable), e.g. '
                         '--parabene C4=12.34 (minutes); stored as mzML '
                         'run-level userParams and in the report')
    ap.add_argument('--lock-ppm', type=float, default=30.0,
                    help='tolerance to match a lock mass (ppm, default: 30)')
    ap.add_argument('--min-lock', type=int, default=MIN_LOCK_POINTS,
                    help='minimum number of lock-mass observations required '
                         'to trust a recalibration (default: %d, the same '
                         'number of calibration points Bruker DataAnalysis '
                         'uses); below it, the model is borrowed from the '
                         'adjacent calibrated file of the same series, or '
                         'the run is left uncalibrated'
                         % MIN_LOCK_POINTS)
    ap.add_argument('--no-prec-correct', action='store_true',
                    help='disable precursor re-phasing on MS1')
    ap.add_argument('--prec-window', type=float, default=PREC_WINDOW_DA,
                    help='precursor search window half-width in Da '
                         '(default: %.1f — the isolation window, a physical '
                         'bound; within it the most intense MS1 apex wins '
                         'automatically)' % PREC_WINDOW_DA)
    ap.add_argument('--polarity', default='auto',
                    choices=['auto', 'positive', 'negative'],
                    help='acquisition polarity (default: auto — detected '
                         'per .d file from the calibrant cluster matches; '
                         'positive and negative runs can be mixed in the '
                         'same folder)')
    ap.add_argument('--report', default=None, help='write a JSON processing report')
    ap.add_argument('--probe', metavar='FOLDER.d', default=None,
                    help='diagnostic mode: dump the calibration-record '
                         'candidates found in this .d folder and exit '
                         '(no conversion, nothing written)')
    args = ap.parse_args(argv)

    print(banner())
    if args.probe:
        # hand-typed console paths often carry quotes, trailing
        # backslashes or spaces -- clean them before validating, and
        # give an actionable error instead of the bare "not a .d".
        d = args.probe.strip().strip('"').strip().rstrip('\\/')
        if d.endswith('.d') and os.path.isdir(d):
            xtr = os.path.join(d, 'analysis.baf_xtr')
            if not os.path.exists(xtr):
                print('probe: analysis.baf_xtr not found in %s' % d)
                return 2
            lines = []

            def log(msg):
                print(msg)
                lines.append(str(msg))
            probe_folder(d, log)
            base = os.path.basename(d)[:-2] or 'd'
            rpt = os.path.join(os.getcwd(),
                               'd2mz_probe_%s_%s.txt'
                               % (base, __version__))
            try:
                with open(rpt, 'w', encoding='utf-8') as f:
                    f.write('\n'.join(lines) + '\n')
                print('probe: report written: %s' % rpt)
            except OSError as e:
                print('probe: could not write report file: %s' % e)
        elif os.path.isdir(os.path.dirname(d) or '.'):
            # the PARENT exists -- say what is actually in it so
            # the right folder name is obvious from one run.
            try:
                cand = [x for x in os.listdir(os.path.dirname(d))
                        if x.lower().endswith('.d')
                        and os.path.isdir(
                            os.path.join(os.path.dirname(d), x))]
            except OSError:
                cand = []
            print('probe: %s is not a .d folder -- %d .d folder(s) '
                  'found next to it:%s'
                  % (d, len(cand), ''.join('\n    ' + c for c in cand
                                           [:10])))
        else:
            print('probe: %s is not a .d folder' % d)
            return 2
        return 0
    if not args.inputs:
        ap.error('the following arguments are required: inputs')

    opts = dict(
        lock_ppm=args.lock_ppm, min_lock=args.min_lock,
        prec_correct=not args.no_prec_correct, prec_window=args.prec_window,
        precision=args.precision, compress=not args.no_compress,
        polarity=args.polarity, parabenes=parse_parabenes(args.parabene),
        uv_traces=args.uv_trace or None)
    report = dict(tool='d2mz', version=__version__, options=vars(args),
                  files=[], scan=[])
    # .zip archives are accepted as inputs -- each archive is
    # extracted to a temporary folder (cleaned up on exit via atexit),
    # scanned for .d folders, and the outputs go NEXT TO the archive
    # (unless --out gives an absolute path).
    resolved_inputs = []
    for p in args.inputs:
        if os.path.isfile(p) and p.lower().endswith('.zip'):
            import atexit
            import shutil
            import tempfile
            import zipfile
            tmp = tempfile.mkdtemp(prefix='d2mz_unzip_')
            atexit.register(shutil.rmtree, tmp, True)   # ignore errors
            print('%s: extracting archive...' % p)
            try:
                with zipfile.ZipFile(p) as zf:
                    unsafe = [i for i in zf.namelist()
                              if i.startswith(('/', '\\')) or '..' in i]
                    if unsafe:
                        print('  !! archive contains unsafe path(s) '
                              '(absolute or ".."): %s (skipped)'
                              % os.path.basename(p))
                        problems.append((os.path.basename(p),
                                         'unsafe path(s) in archive'))
                        continue
                    zf.extractall(tmp)
            except (zipfile.BadZipFile, OSError) as e:
                print('  !! cannot read archive: %s (skipped)' % e)
                problems.append((os.path.basename(p),
                                 'cannot read archive: %s' % e))
                continue
            resolved_inputs.append(
                (tmp, os.path.dirname(os.path.abspath(p)), p))
        else:
            resolved_inputs.append((p, None, None))
    t_run0 = time.time()
    n_converted = 0
    problems = []     # (name, reason) of everything NOT converted,
                      # recapitulated before "process finished."
    zip_targets = []  # output folders to archive with --zip-output
    for input_path, out_override, display in resolved_inputs:
        d_folders, incomplete_d = find_d_folders(input_path, log=print)
        # the scan result goes into the report (durable evidence
        # of what was seen, ignored and why -- the log alone scrolls
        # away) and the ignored .d enter the NOT converted recap
        report.setdefault('scan', []).append(dict(
            input=(display or input_path),
            d_folders_found=len(d_folders),
            ignored_incomplete=[os.path.basename(os.path.normpath(p))
                               for p in incomplete_d]))
        for _p in incomplete_d:
            problems.append(
                (os.path.basename(os.path.normpath(_p)),
                 'no analysis.sqlite and no analysis.baf -- not a '
                 'complete Bruker acquisition (interrupted or partial '
                 'copy?), ignored at scan'))
        if not d_folders:
            if incomplete_d:
                print('only INCOMPLETE .d folder(s) found under %s '
                      '(no analysis.sqlite) (skipped)'
                      % (display or input_path))
            else:
                print('no .d folder found under %s (skipped)'
                      % (display or input_path))
                problems.append((os.path.basename(display or input_path),
                                 'no .d folder found'))
            continue
        # by DEFAULT (no --out) each
        # input folder gets its OWN d2mz_output_<source-name>/ folder
        # created NEXT TO it, so the outputs of several sources never
        # mix in one place. An explicit --out still overrides:
        # relative -> inside each input folder, absolute -> that exact
        # path for the whole queue (GUI single output folder).
        if args.out:
            base = out_override if out_override else input_path
            out_root = args.out if os.path.isabs(args.out) else \
                os.path.join(base, args.out)
        elif out_override:
            # .zip input: name the output after the archive, next to it
            arc = os.path.basename(display or input_path)
            if arc.lower().endswith('.zip'):
                arc = arc[:-4]
            out_root = os.path.join(out_override,
                                    'd2mz_output_' + arc)
        else:
            src = input_path.rstrip('/\\')
            srcname = os.path.basename(src)
            if srcname.lower().endswith('.d'):
                srcname = srcname[:-2]
            parent = os.path.dirname(os.path.abspath(src)) or '.'
            out_root = os.path.join(parent, 'd2mz_output_' + srcname)
        # The output folder mirrors
        # the Zenodo-ready archive layout -- analyses/ holds every
        # mzML/mzXML FLAT, annexes/ holds TSV/ + report_<source>.json
        # + the CITATION files, NOTHING at the output folder root and
        # NOTHING one level above the source folder anymore.
        analyses_dir = os.path.join(out_root, 'analyses')
        annexes_dir = os.path.join(out_root, 'annexes')
        os.makedirs(analyses_dir, exist_ok=True)
        os.makedirs(annexes_dir, exist_ok=True)
        write_citation(annexes_dir)
        # ONE report PER SERIES
        # (report_<source>.json). The old single report.json lived at
        # the same place for every series sharing the parent folder,
        # so each run OVERWROTE the report of the previous series and
        # destroyed its diagnostics.
        _srcname = os.path.basename((display or input_path).rstrip('/\\'))
        if _srcname.lower().endswith('.d'):
            _srcname = _srcname[:-2]
        report_path = os.path.join(annexes_dir,
                                  'report_%s.json' % _srcname)
        print('%s: %d .d folder(s) -> %s (analyses/ + annexes/ with '
              'report: %s)'
              % (display or input_path, len(d_folders), out_root,
                 os.path.basename(report_path)))
        # pre-pass: fit each file's lock-mass model so that a run with a
        # weak or absent calibrant can borrow the model of the nearest
        # calibrated file of the same series (same folder/polarity)
        fits = []
        n_df = len(d_folders)
        for di, d in enumerate(d_folders, 1):
            print('  [%d/%d] calibrating %s'
                  % (di, n_df, os.path.basename(os.path.normpath(d))))
            try:
                fits.append(fit_lock_model(d, opts))
            except Exception as e:
                # Never swallow the reason again -- an old version printed
                # "-> None, lock model weak/absent" for these files,
                # which read like a SUCCESS (misread by everyone)
                fits.append((None, None, None, None, str(e)))
            pol, model = fits[-1][0], fits[-1][1]
            if pol == 'deferred':
                _why = fits[-1][4] or 'calibrant absent'
                print('      -> calibrant absent: DEFERRED (%s). Its '
                      'calibration will be borrowed at the end'
                      % _why)
            elif pol is None:
                print('      -> ERROR at calibration: %s (the file will '
                      'fail again at conversion -- see NOT converted)'
                      % (fits[-1][4] or 'unknown'))
            elif getattr(fits[-1][3], 'segment_only', False):
                print('      -> OK, calibrant-segment calibration '
                      'fitted (no lock-mass observation: normal on '
                      'Impact segment runs)')
            else:
                print('      -> %s, lock model %s'
                      % (pol, 'ok' if model else
                         'weak/absent (may borrow a neighbour\'s)'))
        borrowed = []
        for i, fit_ in enumerate(fits):
            pol, m = fit_[0], fit_[1]
            if m is not None or (fit_[3] is not None and
                                 getattr(fit_[3], 'segment_only', False)):
                # segment-only impact files carry no lock model BY
                # DESIGN (single quadratic): nothing to borrow
                borrowed.append((None, None))
                continue
            donor = None
            for dist in range(1, len(fits)):
                for j in (i - dist, i + dist):
                    if 0 <= j < len(fits) and fits[j][1] is not None \
                            and fits[j][0] == pol:
                        donor = j
                        break
                if donor is not None:
                    break
            borrowed.append((fits[donor][1], os.path.basename(
                os.path.normpath(d_folders[donor]))) if donor is not None
                else (None, None))
        # retention-index anchors: when a 'parabens' run sits in the
        # folder it is CONVERTED FIRST,
        # its retention anchors are measured on its own conversion, and
        # every analysis file converted after it receives them in its
        # mzML. Anchors of each polarity go to that polarity's files;
        # when only one polarity has a parabens run, its anchors serve
        # the opposite-polarity files as well (same series, same
        # injection sequence) -- the ion label in the userParam says
        # which ion they were measured on.
        rstds_by_pol = {}
        rstds_runs = {}    # mode -> [(parabens file name, [anchors])]
        # conversion order: files with their OWN calibration first (in
        # folder order), DEFERRED files LAST so their calibration can be
        # borrowed from a same-mode (same polarity) file of the same
        # folder
        def _is_paraben(i):
            return 'paraben' in os.path.basename(
                os.path.normpath(d_folders[i])).lower()
        order = ([i for i, f_ in enumerate(fits)
                  if f_[0] != 'deferred' and _is_paraben(i)] +
                 [i for i, f_ in enumerate(fits)
                  if f_[0] == 'deferred' and _is_paraben(i)] +
                 [i for i, f_ in enumerate(fits)
                  if f_[0] != 'deferred' and not _is_paraben(i)] +
                 [i for i, f_ in enumerate(fits)
                  if f_[0] == 'deferred' and not _is_paraben(i)])
        cal_donors = {}
        for i, f_ in enumerate(fits):
            if f_[0] in ('deferred', None) or f_[2] is None:
                continue
            if f_[2].get('segment_fit'):
                cal_donors.setdefault(f_[0], []).append(i)
        for di, i in enumerate(order, 1):
            d = d_folders[i]
            fit_ = fits[i]
            bm, bf = borrowed[i]
            name = os.path.basename(os.path.normpath(d))
            t_file0 = time.time()
            prefit = (fit_[2], fit_[3]) if fit_[2] is not None else None
            if fit_[0] == 'deferred':
                m_pol = method_polarity(d)
                donor = None
                if m_pol and cal_donors.get(m_pol):
                    donor = cal_donors[m_pol][0]
                elif len(cal_donors) == 1:
                    donor = next(iter(cal_donors.values()))[0]
                if donor is not None:
                    dname = os.path.basename(os.path.normpath(
                        d_folders[donor]))
                    dcal = dict(fits[donor][2])
                    dcal.pop('_impact_scans', None)
                    dcalib = LockMassCalibrator(
                        [], opts['lock_ppm'], opts['min_lock'])
                    dcalib.polarity = m_pol or fits[donor][0]
                    dcalib.method_polarity = m_pol
                    dcalib.segment_only = True
                    dcalib.calibration_borrowed_from = dname
                    dcalib.n_lock_pos = dcalib.n_lock_neg = 0
                    prefit = (dcal, dcalib)
                    print('  [%d/%d] converting %s (calibration borrowed '
                          'from %s, same mode/polarity)'
                          % (di, n_df, name, dname))
                else:
                    # 1. segment 2, 2. another file. NO other fit
                    # path -- a deferred
                    # file with no same-mode donor is NOT converted and
                    # the failure is reported.
                    print('  [%d/%d] %s: no calibrant segment and no '
                          'same-mode donor -> NOT converted (reported)'
                          % (di, n_df, name))
                    _why = fit_[4] if len(fit_) > 4 else None
                    problems.append(
                        (name, 'no calibrant segment (Segment==2) and '
                         'no same-mode donor%s'
                         % (': ' + _why if _why else '')))
                    report['files'].append(
                        dict(name=name,
                             error='no calibrant segment (Segment==2) '
                                   'and no same-mode donor: not '
                                   'converted'))
                    continue
            else:
                print('  [%d/%d] converting %s' % (di, n_df, name))
            try:
                specs, prec_table, calib, uv = convert_one(
                    d, opts, borrow_model=bm, borrow_from=bf,
                    prefit=(prefit[0], prefit[1]) if prefit else None)
            except Exception as e:
                print('  ERROR in %s: %s' % (name, e))
                report['files'].append(dict(name=name, error=str(e)))
                problems.append((name, 'ERROR: %s' % e))
                continue
            if 'paraben' in name.lower():
                # parabens run converted first: measure the retention
                # anchors ON THIS CONVERSION (no second pass) so every
                # analysis file converted after receives them
                pol_a = getattr(calib, 'polarity', None) or 'positive'
                found = measure_parabene_anchors(specs, pol_a)
                if found:
                    # every parabens run of the mode contributes; the
                    # AVERAGE per anchor name is what gets stamped
                    rstds_runs.setdefault(pol_a, []).append(
                        (name, found))
                    rstds_by_pol[pol_a] = average_parabene_anchors(
                        rstds_runs[pol_a])
                    nruns = len(rstds_runs[pol_a])
                    print('  parabens run (%s): %d anchors -> mode %s '
                          'now averaged over %d run(s) (%s)'
                          % (name, len(found), pol_a, nruns,
                             ', '.join('%s %.3f min' % (
                                 a['name'], a['rt_sec'] / 60.0)
                                 for a in rstds_by_pol[pol_a])))
                else:
                    print('  WARNING: %s: no parabene anchor matched '
                          'in this run' % name)
            pol_f = getattr(calib, 'polarity', None) or 'positive'
            rstds = rstds_by_pol.get(pol_f)
            if rstds is None:
                # only one polarity had a parabens run: its anchors are
                # used for this opposite-polarity file as well
                other = [v for p, v in rstds_by_pol.items() if p != pol_f]
                rstds = other[0] if other else []
            # mzML and mzXML are
            # written FLAT in the analyses/ subfolder of the output
            # folder, and the precursor tables live in its
            # annexes/TSV/ subfolder -- the folder layout mirrors the
            # Zenodo-ready archive.
            tsv_dir = os.path.join(out_root, 'annexes', 'TSV')
            os.makedirs(tsv_dir, exist_ok=True)
            base = os.path.join(out_root, 'analyses',
                                name.replace('.d', ''))
            written = []
            gzsfx = '.gz' if args.gzip else ''
            if args.format in ('mzml', 'both'):
                write_mzml(base + '.mzML' + gzsfx, specs, opts, name,
                           start_stamp=acquisition_start_stamp(d),
                           uv_traces=uv['traces'], uv_apexes=uv['apexes'],
                           retention_standards=rstds,
                           extra_traces=uv['extras'],
                           lc_method=uv.get('lc'))
                written.append(base + '.mzML' + gzsfx)
            if args.format in ('mzxml', 'both'):
                write_mzxml(base + '.mzXML' + gzsfx, specs, opts, name,
                            start_stamp=acquisition_start_stamp(d))
                written.append(base + '.mzXML' + gzsfx)
            if prec_table:
                with open(os.path.join(
                        tsv_dir,
                        name.replace('.d', '') + '.tsv'), 'w') as f:
                    # NO ms1_intensity
                    # column in the report files (only scan, RT, m/z)
                    f.write('scan\trt_min\tdeclared_mz\tcorrected_mz\t'
                            'delta_ppm\tce_window\n')
                    for row in prec_table:
                        f.write('%d\t%.3f\t%.5f\t%s\t%s\t%s\n' % (
                            row['scan'], row['rt_min'], row['declared'],
                            '%.5f' % row['corrected'] if row['corrected'] is not None else 'NA',
                            '%.2f' % row['delta_ppm'] if row['delta_ppm'] is not None else 'NA',
                            (row.get('ce_window') or 'NA') + ' eV' if row.get('ce_window') else 'NA'))
            n_ms2 = sum(1 for s in specs if s['ms_level'] == 2)
            n_fix = sum(1 for r in prec_table if r['corrected'] is not None)
            print('    %d scans (%d MS2), %d precursors re-phased -> %s'
                  % (len(specs), n_ms2, n_fix,
                     ', '.join(os.path.basename(w) for w in written)))
            _ms2_source = ('analysis.sqlite (opened in DataAnalysis)'
                           if os.path.exists(os.path.join(
                               d, 'analysis.sqlite'))
                           else 'baf fallback (never opened)')
            _ms2_warning = None
            if _ms2_source.startswith('baf') and n_ms2 == 0:
                _ms2_warning = ('no precursor metadata matched -- '
                                'the baf fallback decoded every scan '
                                'as MS1; run d2mz --probe on this .d')
            _ctr = getattr(calib, 'calibration_trace', None) or {}
            report['files'].append(dict(
                name=name, input=display or input_path, n_scans=len(specs), n_ms2=n_ms2,
                ms2_source=_ms2_source, ms2_warning=_ms2_warning,
                acquisition_mode=opts.get('acquisition_mode'),
                mode_warning=opts.get('mode_warning'),
                duration_s=round(time.time() - t_file0, 1),
                precursors_rephased=n_fix,
                n_lock_points=len(calib.points),
                detected_polarity=getattr(calib, 'polarity', None),
                method_polarity=getattr(calib, 'method_polarity', None),
                lock_points_positive=getattr(calib, 'n_lock_pos', None),
                lock_points_negative=getattr(calib, 'n_lock_neg', None),
                lock_rms_ppm_before=calib.rms_before(),
                lock_rms_ppm_after=calib.rms_after(),
                lock_mean_ppm_after=calib.mean_ppm_after(),
                lock_stdev_ppm_after=calib.stdev_ppm_after(),
                calibration_applied=_ctr.get('applied'),
                calibration_mode=_ctr.get('mode'),
                calibration_n_targets=_ctr.get('n_targets'),
                calibration_ppm_median=_ctr.get('ppm_median'),
                calibration_ppm_max=_ctr.get('ppm_max'),
                calibration_borrowed_from=getattr(
                    calib, 'calibration_borrowed_from', None)
                or getattr(calib, 'borrowed_from', None),
                uv_traces={l: len(p) for l, p in sorted(uv['traces'].items())},
                extra_traces=[ds for ds, _wl, _p in uv['extras']],
                lc_flow=(', '.join('%g min: %g ml/min' % fv
                                   for fv in (uv.get('lc') or {}).get(
                                       'flow', []))
                          if uv.get('lc') else None),
                lc_gradient=(lc_gradient_str(uv['lc'])
                             if uv.get('lc') else None),
                lc_solvents=(uv['lc'].get('solvents')
                             if uv.get('lc') else None),
                uv_apex_spectra=len(uv['apexes']),
                n_retention_standards=len(rstds),
                written=[os.path.basename(w) for w in written]))
            n_converted += 1
        # retention standards: per mode, the AVERAGED anchors actually
        # stamped into the mzMLs + their RT variation across the
        # parabens runs of that mode
        report['retention_standards'] = {}
        for pol_a, runs in sorted(rstds_runs.items()):
            per_name = {}
            for fname, found in runs:
                for a in found:
                    per_name.setdefault(a['name'], dict(
                        name=a['name'], ion=a['ion'], ri=a['ri'],
                        rts=[], sources=[], per_run=[]))
                    per_name[a['name']]['rts'].append(a['rt_sec'])
                    per_name[a['name']]['sources'].append(fname)
                    # Every individual measurement kept: a drifted run
                    # can be
                    # identified and removed later, run by run
                    per_name[a['name']]['per_run'].append(dict(
                        file=fname, rt_min=round(a['rt_sec'] / 60.0, 4)))
            anchors = []
            for name in sorted(per_name, key=lambda n: per_name[n]['ri']):
                e = per_name[name]
                rts = e['rts']
                mean = sum(rts) / len(rts)
                spread = (max(rts) - min(rts)) if len(rts) > 1 else 0.0
                stdev = (math.sqrt(sum((x - mean) ** 2 for x in rts) /
                                   len(rts)) if len(rts) > 1 else 0.0)
                anchors.append(dict(
                    name=name, ion=e['ion'], ri=e['ri'],
                    n_runs=len(rts), sources=e['sources'],
                    per_run=e['per_run'],
                    rt_mean_min=round(mean / 60.0, 4),
                    rt_min=round(min(rts) / 60.0, 4),
                    rt_max=round(max(rts) / 60.0, 4),
                    rt_spread_s=round(spread, 2),
                    rt_stdev_s=round(stdev, 2)))
            report['retention_standards'][pol_a] = anchors
        # 2026-11-25_cp: stamp the batch blocks BEFORE the per-series
        # dump -- report_<series>.json used to be written before
        # quality_gate/calibration were computed, so it carried neither
        report['calibration'] = _calibration_summary(report['files'])
        _gwarns = _file_warnings(report)
        report['quality_gate'] = dict(
            n_converted=n_converted, n_errors=len(problems),
            n_warnings=len(_gwarns),
            errors=[{'file': _nm, 'reason': _rs} for _nm, _rs in problems],
            warnings=[{'file': _nm, 'reason': _rs}
                      for _nm, _rs in _gwarns])
        report['duration_s'] = round(time.time() - t_run0, 1)
        # the per-series report sits in the annexes/ subfolder
        # of the output folder, next to the CITATION files and TSV/
        with open(report_path, 'w') as f:
            json.dump(report, f, indent=2)
        if args.zip_output:
            # archive AFTER report.json is written so the .zip
            # contains the complete output. Shared --out: the same
            # folder is appended once per input and deduplicated below,
            # so the single archive is written at the very end with the
            # whole queue inside. The archive is structured for a
            # DIRECT Zenodo upload and the output folder on disk now
            # MIRRORS it exactly:
            # <root>/analyses/ = every mzML/mzXML FLAT, <root>/annexes/
            # = TSV/ + report_<source>.json + CITATION files of every
            # series that landed in this output folder. The walk below
            # therefore archives the folder AS IS.
            zip_targets.append(out_root)
    if args.zip_output:
        import zipfile as _zipf

        for _oroot in sorted(set(zip_targets)):
            _zpath = os.path.join(os.path.dirname(_oroot),
                                  os.path.basename(_oroot) + '.zip')
            _zroot = os.path.basename(_oroot)
            try:
                _n_an, _n_ann = 0, 0
                with _zipf.ZipFile(_zpath, 'w',
                                   _zipf.ZIP_DEFLATED) as _zf:
                    for _dp, _dirs, _fns in os.walk(_oroot):
                        _dirs.sort()
                        for _fn in sorted(_fns):
                            _fp = os.path.join(_dp, _fn)
                            _rel = os.path.relpath(_fp, _oroot)
                            _arc = '%s/%s' % (_zroot,
                                              _rel.replace('\\', '/'))
                            _zf.write(_fp, _arc)
                            if _arc.split('/')[1] == 'analyses':
                                _n_an += 1
                            else:
                                _n_ann += 1
                print('archived: %s (%d analyses + %d annexes) -- flat '
                      'files kept for MZmine, Zenodo-ready (analyses/ '
                      ' + annexes/)' % (_zpath, _n_an, _n_ann))
            except OSError as e:
                print('  !! cannot write %s: %s (skipped)' % (_zpath, e))
                problems.append((os.path.basename(_zpath),
                                 'cannot write archive: %s' % e))
    # QUALITY GATE -- a run with
    # errors or silent-quality flags must end UNMISTAKABLY, never on a
    # bare "process finished." Recapitulate BOTH the NOT-converted
    # and the per-file quality flags mined from the report:
    # borrowed calibration, polarity mismatch, and the MS2-loss
    # warning. Errors fail the run (exit code 1); warnings recap but
    # keep exit 0.
    # 2026-11-25_cp: warn list factored (_file_warnings), shared
    # between the per-series report dumps and this final gate
    warns = _file_warnings(report)
    report['calibration'] = _calibration_summary(report['files'])
    report['quality_gate'] = dict(
        n_converted=n_converted, n_errors=len(problems),
        n_warnings=len(warns),
        errors=[{'file': _nm, 'reason': _rs} for _nm, _rs in problems],
        warnings=[{'file': _nm, 'reason': _rs} for _nm, _rs in warns])
    if args.report:
        report['duration_s'] = round(time.time() - t_run0, 1)
        with open(args.report, 'w') as f:
            json.dump(report, f, indent=2)
    # Expose the outcome so GUI wrappers can react without
    # scraping stdout (LAST_RUN is reset at every main() call).
    global LAST_RUN
    LAST_RUN = {'n_warnings': len(warns),
                'warnings': list(warns),
                'n_errors': len(problems),
                'errors': list(problems)}
    print('done. %d .d folder(s) converted in %.1f s.'
          % (n_converted, time.time() - t_run0))
    # The run used to end on a
    # bare "done." even when files had been skipped along the way, so
    # there was no way to see WHY some requested files were
    # missing from the output folder. Recapitulate every skip with
    # its reason, then mark the end unambiguously.
    if problems:
        print('NOT converted (%d):' % len(problems))
        for _nm, _rs in problems:
            print('    %s: %s' % (_nm, _rs))
    if warns:
        print('QUALITY WARNINGS (%d):' % len(warns))
        for _nm, _rs in warns:
            print('    %s: %s' % (_nm, _rs))
        # warnings used to end in a
        # quiet list; the run then closed on a bare "process finished",
        # and a degraded batch (e.g. 0-MS2 fallback files) looked
        # identical to a clean one at a glance. Close on a loud banner.
        print('=' * 66)
        print('!! %d QUALITY WARNING(S) -- CONVERSION INCOMPLETE OR '
              'DEGRADED: inspect the list above. !!' % len(warns))
        print('=' * 66)
    if problems:
        print('=' * 66)
        print('!! %d ERROR(S) -- OUTPUT INCOMPLETE: %d/%d converted. '
              'Fix and rerun. !!'
              % (len(problems), n_converted, n_converted + len(problems)))
        print('=' * 66)
        print('process finished WITH ERRORS.')
        return 1
    print('process finished. %d file(s) converted, no error '
          '(%d quality warning(s)).' % (n_converted, len(warns)))
    return 0


if __name__ == '__main__':
    # No argument at all (typically a double-click on Windows): open an
    # interactive session: the user drags folder(s) into the console
    # window (this pastes their paths) and presses Enter to convert.
    if len(sys.argv) == 1:
        print(banner())
        print('  Interactive mode: DRAG one or several folders (or a\n'
              '  single .d folder) into THIS window, then press Enter.\n'
              '  Converted to mzML 32-bit, same options as the batch run.\n'
              '  Output goes to <each input folder>\\d2mz_output.\n'
              '  Empty line = quit.')
        while True:
            try:
                raw = input('\n  Folder(s) to convert (drag here) > ').strip()
            except (EOFError, KeyboardInterrupt):
                break
            if not raw:
                break
            # accept several dragged paths at once: "p1" "p2" or bare
            paths = re.findall(r'"([^"]+)"', raw)
            paths += [t for t in re.sub(r'"[^"]+"', ' ', raw).split()]
            paths = [p.rstrip('\\/') for p in paths if p]
            if not paths:
                continue
            print('  Converting %d input(s): %s'
                  % (len(paths), ', '.join(paths)))
            try:
                main(paths + ['--format', 'mzml', '--precision', '32'])
            except SystemExit as e:
                if e.code:
                    print('  (exit code %s)' % e.code)
            except Exception as e:
                print('  ERROR: %s' % e)
            print('  Done. You can drag another folder, or press Enter to quit.')
        print('  Bye.')
        sys.exit(0)
    try:
        sys.exit(main())
    finally:
        # In a PyInstaller --onefile exe, the console window closes at
        # the end of the run (or on an argparse error). On Windows,
        # keep it open so the final recap / errors stay readable.
        if getattr(sys, 'frozen', False) and os.name == 'nt':
            try:
                input('\nPress Enter to close... ')
            except (EOFError, KeyboardInterrupt):
                pass


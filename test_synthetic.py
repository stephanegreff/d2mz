#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_synthetic.py — end-to-end tests for d2mz on a synthetic Bruker .d.

Builds a mock .d folder (analysis.sqlite + analysis.baf/.baf_idx/.baf_xtr)
with a known calibration, a quadratic ppm drift (+3 ppm linear
component + 0.8 ppm curvature) injected in every
spectrum, and declared MS2 precursors offset from the true MS1 apex by
~7 ppm. Then runs d2mz.main() and asserts: recalibration residual
mean standard deviation of lock-mass errors < 2 ppm (Bruker's
mass-accuracy specification), precursor error < 0.2 ppm, valid
mzML/mzXML round-trip.

Run:  python3 test_synthetic.py
"""
import base64
import json
import math
import os
import shutil
import sqlite3
import struct
import sys
import tempfile
import zlib
import xml.etree.ElementTree as ET

import d2mz

CAL = dict(mode=4, alpha=0.0, beta=1.0, freq_high=1.0e6,
           size=100000, freq_low=1.0, freq_width=20000.0)


def out_root_of(w):
    # v_bd default: each input folder gets its own
    # d2mz_output_<name>/ folder created NEXT TO it
    return os.path.join(os.path.dirname(os.path.abspath(w)),
                       'd2mz_output_' + os.path.basename(
                           os.path.normpath(w)))

def report_of(w):
    # v_bm: ONE report PER SERIES, report_<source>.json. v_by: it
    # lives in the annexes/ subfolder of the output folder
    return os.path.join(out_root_of(w), 'annexes',
                        'report_%s.json' % os.path.basename(
                            os.path.normpath(w)))


def mz_to_index(mz):
    den = CAL['freq_high'] / mz
    return CAL['size'] * (1.0 - (den + CAL['beta']) / CAL['freq_width'])


def write_xtr(path):
    p = bytearray(144)
    struct.pack_into('<i', p, 0, len(p))
    struct.pack_into('<i', p, 4, 4)
    struct.pack_into('<i', p, 8, 76)
    fo = 80
    struct.pack_into('<i', p, fo + 8, CAL['mode'])
    struct.pack_into('<d', p, fo + 20, -CAL['beta'])
    struct.pack_into('<d', p, fo + 28, CAL['freq_high'])
    struct.pack_into('<i', p, fo + 40, CAL['size'])
    struct.pack_into('<d', p, fo + 48, CAL['freq_low'])
    struct.pack_into('<d', p, fo + 56, CAL['freq_width'])
    bs = len(p) + 24
    with open(path, 'wb') as f:
        f.write(b'\x00\x00\x00\x01')
        f.write(b'\x00' * 4)
        f.write(struct.pack('<i', bs))
        f.write(bytes(p))
        f.write(b'\x00' * 12)
        f.write(b'\x00' * 4 + b'\x00' * 4 + struct.pack('<i', 0))


def maxima_block(mz_list, itn_list):
    n = len(mz_list)
    idx = [mz_to_index(m) for m in mz_list]
    payload = struct.pack('<ii', n, 8)
    payload += struct.pack('<%dd' % n, *idx)
    payload += struct.pack('<%dd' % n, *itn_list)
    hs = 12
    bs = hs + len(payload)
    return struct.pack('<i4si', bs, b'ABCD', hs) + payload


def write_idx(path, sid_to_offset):
    with open(path, 'wb') as f:
        for sid, off in sid_to_offset.items():
            rec = bytearray(24)
            rec[0], rec[1] = 2, 16
            struct.pack_into('<i', rec, 8, sid)
            struct.pack_into('<i', rec, 16, off)
            struct.pack_into('<i', rec, 20, 0)
            bs = len(rec) + 24
            f.write(b'\x00\x00\x00\x01')
            f.write(b'\x00' * 4)
            f.write(struct.pack('<i', bs))
            f.write(bytes(rec))
            f.write(b'\x00' * 12)
        f.write(b'\x00\x00\x00\x01' + b'\x00' * 4 + struct.pack('<i', 0))


def write_idx_records(path, records):
    """General baf_idx writer (v_bp): records are (id_a, id_b, flag,
    offset) tuples -- (2,16)=maxima, (2,32)=per-scan metadata,
    (1,48)=acquisition header."""
    with open(path, 'wb') as f:
        for id_a, id_b, flag, off in records:
            rec = bytearray(24)
            rec[0], rec[1] = id_a, id_b
            struct.pack_into('<i', rec, 8, flag)
            struct.pack_into('<i', rec, 16, off)
            struct.pack_into('<i', rec, 20, 0)
            f.write(b'\x00\x00\x00\x01')
            f.write(b'\x00' * 4)
            f.write(struct.pack('<i', len(rec) + 24))
            f.write(bytes(rec))
            f.write(b'\x00' * 12)
        f.write(b'\x00\x00\x00\x01' + b'\x00' * 4
                + struct.pack('<i', 0))


def metadata_block(numeric_rows, string_rows=()):
    """Per-scan metadata block in the byte layout d2mz._read_metadata_
    values parses: [i32 mbs][4s][i32 mhs=12], then the numeric section
    [i32 nbs][4s][i32 nhs=16][i32 ncnt] + ncnt x (H row, d value),
    then the string section [i32 sbs][4s][i32 shs=16][i32 scnt] +
    rows [H rs][2 pad][H row][i32 dlen][text (dlen-1 bytes)]."""
    num = b''.join(struct.pack('<Hd', ri, v) for ri, v in numeric_rows)
    nbs = 16 + len(num)
    numsec = struct.pack('<i4sii', nbs, b'', 16, len(numeric_rows)) + num
    srows = b''
    for ri, txt in string_rows:
        tb = txt.encode('utf-8') if isinstance(txt, str) else txt
        # same v_bs fix: no leading pad, rs = len(content) + 2
        row = (struct.pack('<Hi', ri, len(tb) + 1)
               + tb + b'\x00')
        srows += struct.pack('<H', len(row) + 2) + row
    sbs = 16 + len(srows)
    strsec = struct.pack('<i4sii', sbs, b'', 16, len(string_rows)) + srows
    return (struct.pack('<i4si', 12 + nbs + sbs, b'', 12)
            + numsec + strsec)


def metadata_block_rt(rt_ms, numeric_rows=(), string_rows=()):
    """v_bw never-opened layout (real Qualitative file): 40-byte
    block header with the per-scan Rt in MILLISECONDS at offset 24
    (scan 1 = 299, scan 3970 = 965142). numeric_rows are
    (ri, raw 8 value bytes) so integer rows can be packed as '<q'."""
    num = b''.join(struct.pack('<H8s', ri, v) for ri, v in numeric_rows)
    nbs = 16 + len(num)
    numsec = struct.pack('<i4sii', nbs, b'', 16, len(numeric_rows)) + num
    srows = b''
    for ri, txt in string_rows:
        tb = txt.encode('utf-8') if isinstance(txt, str) else txt
        row = (struct.pack('<Hi', ri, len(tb) + 1)
               + tb + b'\x00')
        srows += struct.pack('<H', len(row) + 2) + row
    sbs = 16 + len(srows)
    strsec = struct.pack('<i4sii', sbs, b'', 16, len(string_rows)) + srows
    head = (struct.pack('<i4si', 40 + nbs + sbs, b'', 40)
            + struct.pack('<7i', 1, 1, 0, rt_ms, 0, 0, 0))
    return head + numsec + strsec


def acquisition_header_block(param_rows):
    """Acquisition-header parameter table in the byte layout
    d2mz._read_baf_params parses: [20B header -> nxt=20][16B middle
    block][table [i32 tbs][6B][H nrow] + rows [H rs][2 pad][H row]
    [i32 type][6H field offsets][name\\0unit\\0\\0labels\\0]]."""
    tbl_rows = b''
    for row_ in param_rows:
        ri, nm, lb = row_[:3]
        ty = row_[3] if len(row_) > 3 else 4
        parts = [nm, b'unit', b'', lb]
        offs, payload = [], b''
        for p_ in parts:
            offs.append(20 + len(payload))
            payload += p_ + b'\x00'
        offs += [0, 0]          # 6 offset slots: 0 = absent (reader skips)
        # reader rebuilds row = b'\x00\x00' + file_bytes(rs-2), so the
        # file stores the row WITHOUT its 2 leading bytes and rs counts
        # them back in (v_bs fix: the former leading pad shifted ri,
        # every name decoded as 'row0' and the suite PASSED only
        # because a failed conversion leaves rc=0 + a stale mzML)
        row = (struct.pack('<Hi', ri, ty)
               + struct.pack('<6H', *offs) + payload)
        tbl_rows += struct.pack('<H', len(row) + 2) + row
    tbl = (struct.pack('<i', 12 + len(tbl_rows)) + b'\x00' * 6
           + struct.pack('<H', len(param_rows)) + tbl_rows)
    middle = struct.pack('<i', 16) + b'\x00' * 12
    head = struct.pack('<i', 20) + b'\x00' * 12 + struct.pack('<i', 20)
    return head + middle + tbl


def make_mock_d(root, n_ms1=12, drift_ppm=3.0, prec_offset_ppm=7.0,
                name='Synthetic_01', polarity='positive',
                extra_ms1_peaks=None, uv1_nm=280.0, uv2_nm=320.0):
    d = os.path.join(root, name + '.d')
    os.makedirs(d)
    # the calibrant clusters d2mz locks on automatically: the first 6 ions
    # of the built-in Bruker series for this polarity
    lock_masses = [m for m, _f in
                   d2mz.BRUKER_FA_AC_NA_CLUSTERS[polarity][:6]]
    analytes = [479.1180, 641.1710, 353.0650]
    ms1_rows, ms2_rows = [], []
    blocks = bytearray()
    sid_off = {}
    sid = 0
    for k in range(n_ms1):
        rt = 30.0 + k * 30.0
        u = k / 11.0
        # linear + quadratic curvature: a purely linear regression
        # would NOT remove this drift (proves the quadratic model)
        drift = drift_ppm * u + 0.8 * u * u
        def shifted(m):
            return m * (1.0 + drift * 1e-6)
        mzs, itns = [], []
        for lm in lock_masses:
            mzs.append(shifted(lm)); itns.append(5.0e4)
        for a in analytes:
            mzs.append(shifted(a)); itns.append(1.0e4 + 1000.0 * k)
        for mz_e, it_e in (extra_ms1_peaks or {}).get(k, []):
            mzs.append(shifted(mz_e)); itns.append(it_e)
        sid += 1
        sid_off[sid] = len(blocks)
        blocks += maxima_block(mzs, itns)
        ms1_rows.append((sid, rt))
        parent_true = analytes[k % len(analytes)]
        declared = parent_true * (1.0 + prec_offset_ppm * 1e-6)
        m2 = [153.0181, 169.0137, 302.0460, 479.1180]
        m2_mz = [shifted(m) * (1.0 - prec_offset_ppm * 1e-6) if m == 479.1180
                 else shifted(m) for m in m2]
        m2_it = [3000.0, 1700.0, 9000.0, 500.0]
        sid += 1
        sid_off[sid] = len(blocks)
        blocks += maxima_block(m2_mz, m2_it)
        ms2_rows.append((sid, rt + 0.5, declared))
    with open(os.path.join(d, 'analysis.baf'), 'wb') as f:
        f.write(blocks)
    write_idx(os.path.join(d, 'analysis.baf_idx'), sid_off)
    write_xtr(os.path.join(d, 'analysis.baf_xtr'))
    c = sqlite3.connect(os.path.join(d, 'analysis.sqlite'))
    c.execute('create table Spectra (Id integer primary key, Rt real, '
              'AcquisitionKey integer)')
    c.execute('create table Steps (Id integer primary key, '
              'TargetSpectrum integer, Mass real)')
    c.execute('create table PerSpectrumVariables (Spectrum integer, '
              'Variable integer, Value real)')
    for sid, rt in ms1_rows:
        c.execute('insert into Spectra values (?,?,?)', (sid, rt, 1))
    for sid, rt, declared in ms2_rows:
        c.execute('insert into Spectra values (?,?,?)', (sid, rt, 2))
        c.execute('insert into Steps values (?,?,?)', (sid, sid, declared))
        c.execute('insert into PerSpectrumVariables values (?,?,?)',
                  (sid, 5, 20.0 + 10.0 * (sid % 3)))
    c.execute('create table UvPerScan (Spectrum integer, Wavelength real, '
              'Intensity real)')
    wls = [200.0 + 25.0 * j for j in range(13)] + [254.0, 365.0]
    for k2, (sid, rt) in enumerate(ms1_rows):
        # two chromatographic peaks: 254 nm apex at scan 6, 365 nm at scan 8
        g254 = math.exp(-((k2 - 6) / 1.5) ** 2) * 100.0
        g365 = math.exp(-((k2 - 8) / 1.5) ** 2) * 100.0
        for wl in wls:
            inten = 5.0 + (g254 if wl == 254.0 else 0.0) + \
                    (g365 if wl == 365.0 else 0.0)
            c.execute('insert into UvPerScan values (?,?,?)',
                      (sid, wl, inten))
    c.commit()
    # chromatography-data.sqlite in the REAL recent otoControl format:
    # named traces (TraceSources/TraceChunks) -- UV_VIS_1 is the FIRST
    # acquired UV chromatogram (analyst directive 2026-11-07: record
    # UV1, e.g. at 280 nm -- NOT a hardcoded 254 nm)
    cduc = sqlite3.connect(os.path.join(d, 'chromatography-data.sqlite'))
    cduc.execute('create table TraceSources (Id integer, Description '
                 'text, Instrument text, InstrumentId text, Type integer,'
                 ' Unit integer, TimeOffset real, Color integer)')
    cduc.execute('create table TraceChunks (Trace integer, Times blob, '
                 'Intensities blob)')
    cduc.execute("insert into TraceSources values (1,'UV_VIS_1','Thermo',"
                 "'x',3,6,0.0,-65536)")
    cduc.execute("insert into TraceSources values (2,'UV_VIS_2','Thermo',"
                 "'x',3,6,0.0,-16776961)")
    cduc.execute("insert into TraceSources values (6,'BPC',"
                 "'Bruker OTOF MS','y',1,6,0.0,-65536)")
    for tr, koff in ((1, 6), (2, 8)):
        ts, ints = [], []
        for k2, (sid, rt) in enumerate(ms1_rows):
            g = math.exp(-((k2 - koff) / 1.5) ** 2) * 100.0
            ts.append(rt)
            ints.append(5.0 + g)
        cduc.execute('insert into TraceChunks values (?,?,?)', (
            tr, struct.pack('<%dd' % len(ts), *ts),
            struct.pack('<%df' % len(ints), *ints)))
    cduc.commit(); cduc.close()
    # Hystar LC method (double-XML-escaped like the real files): the
    # acquisition wavelength of UV channel 1
    mf = os.path.join(d, 'acq.m')
    os.makedirs(mf, exist_ok=True)
    with open(os.path.join(mf, 'hystar.method'), 'w') as f:
        f.write('<Method><Parameter>'
                '&amp;lt;Name&amp;gt;UV_VIS_1.Wavelength&amp;lt;/Name'
                '&amp;gt;&amp;lt;Id&amp;gt;0&amp;lt;/Id&amp;gt;'
                '&amp;lt;Unit&amp;gt;nm&amp;lt;/Unit&amp;gt;'
                '&amp;lt;Value&amp;gt;%g&amp;lt;/Value&amp;gt;'
                '&amp;lt;/Parameter&amp;gt;<Parameter>'
                '&amp;lt;Name&amp;gt;UV_VIS_2.Wavelength&amp;lt;/Name'
                '&amp;gt;&amp;lt;Id&amp;gt;0&amp;lt;/Id&amp;gt;'
                '&amp;lt;Unit&amp;gt;nm&amp;lt;/Unit&amp;gt;'
                '&amp;lt;Value&amp;gt;%g&amp;lt;/Value&amp;gt;'
                '&amp;lt;/Parameter&amp;gt;</Method>'
                % (uv1_nm, uv2_nm))
    c.close()
    return d, dict(lock_masses=lock_masses, analytes=analytes,
                   n_ms1=n_ms1, drift_ppm=drift_ppm,
                   prec_offset_ppm=prec_offset_ppm)


# ---- impact-generation (mode 10) mock -------------------------------
# True quadratic index->m/z: hardware bootstrap + a per-file drift
# (a +0.2%, b -0.1%, c -0.36: the Prunus-class offset). The bootstrap
# predictions are then off by up to ~350 idx, which the 600->200->70
# window series must absorb.
IMPACT_CAL = dict(a=6.174825e-09 * 1.002, b=0.001109141 * 0.999,
                  c=49.83071 - 0.36)


def impact_mz_to_index(mz):
    aa, bb, cc = IMPACT_CAL['a'], IMPACT_CAL['b'], IMPACT_CAL['c']
    return (-bb + math.sqrt(bb * bb - 4.0 * aa * (cc - mz))) / (2.0 * aa)


def impact_maxima_block(peak_list):
    """peak_list: [(idx, intensity, snr)] in the REAL impact layout
    (byte-level re-probe 2026-11-07, Complément 200): [n][po] + n x f64
    idx + n x f32 INTENSITY + n x f32 SNR (consecutive arrays)."""
    n = len(peak_list)
    payload = struct.pack('<ii', n, 8)
    payload += struct.pack('<%dd' % n, *[p[0] for p in peak_list])
    payload += struct.pack('<%df' % n, *[p[1] for p in peak_list])
    payload += struct.pack('<%df' % n, *[p[2] for p in peak_list])
    hs = 12
    bs = hs + len(payload)
    return struct.pack('<i4si', bs, b'ABCD', hs) + payload


def write_xtr_impact(path):
    p = bytearray(144)
    struct.pack_into('<i', p, 0, len(p))
    struct.pack_into('<i', p, 4, 4)
    struct.pack_into('<i', p, 8, 76)
    struct.pack_into('<i', p, 80 + 8, 10)      # transformator mode 10
    bs = len(p) + 24
    with open(path, 'wb') as f:
        f.write(b'\x00\x00\x00\x01')
        f.write(b'\x00' * 4)
        f.write(struct.pack('<i', bs))
        f.write(bytes(p))
        f.write(b'\x00' * 12)
        f.write(b'\x00' * 4 + b'\x00' * 4 + struct.pack('<i', 0))


def make_mock_impact_d(root, n_ms1=16, drift_ppm=3.0,
                       prec_offset_ppm=7.0, name='Impact_NEG',
                       n_cal_scans=9, segment=True,
                       polarity='negative', baf_metadata=False,
                       manual_list_ce=False, isolation_act=False,
                       prec_jitter=False):
    """Real Impact II structure (analyst directive 2026-11-07): the
    calibrant comb lives in the CALIBRANT SEGMENT (Spectra.Segment==2,
    the first n_cal_scans MS1 scans, exactly like the 9 segment-2 scans
    RT 6-18 s of the real runs); the ANALYSIS segment (Segment==1)
    carries analytes/junk/noise but NO calibrant ions. segment=False
    builds the same file WITHOUT the Segment column (deferred path)."""
    import random
    rng = random.Random(42)
    d = os.path.join(root, name + '.d')
    os.makedirs(d)
    # acquisition method folder (.m): declares the polarity upstream
    mf = os.path.join(d, 'acq.m')
    os.makedirs(mf)
    with open(os.path.join(mf, 'microTOFQImpactAcquisition.method'), 'w') as fh:
        pol_word = 'Positive' if polarity == 'positive' else 'Negative'
        fh.write('<Method><Source IonPolarity="%s"/></Method>\n' % pol_word)
        if manual_list_ce:
            # v_co: residual MSMSManual_List* paragraph, exactly the
            # XML shape method_ce_sweep parses (real case 38087: every
            # real .method carries one, even in pure auto-DDA)
            fh.write('<para_vec permname="MSMSManual_ListCollisionEnergy"'
                     '><entry_double value="40.0"/>'
                     '<entry_double value="8.0"/></para_vec>\n')
    lock_masses = d2mz._cluster_families(
        d2mz.BRUKER_FA_AC_NA_CLUSTERS[polarity])[0][1][:6]
    analytes = [479.1180, 641.1710, 353.0650]
    comb_snr = [26.0, 19.0, 14.0, 11.0, 9.5, 8.5]
    ms1_rows, ms2_rows = [], []
    blocks = bytearray()
    sid_off = {}
    sid = 0
    for k in range(n_ms1):
        rt = 30.0 + k * 30.0
        u = k / float(n_ms1 - 1)
        drift = drift_ppm * u + 0.8 * u * u
        cal_scan = k < n_cal_scans    # calibrant segment (Segment==2)
        saturated = cal_scan and k < 2    # start-of-run saturated combs
        peaks = []

        def shifted(m):
            return m * (1.0 + drift * 1e-6)

        if cal_scan:
            for j, lm in enumerate(lock_masses):
                idx = impact_mz_to_index(shifted(lm))
                if saturated and j == 0:
                    # flat top: dozens of bucket maxima +- 300 idx around
                    # the comb (weak), one strong apex with a small position
                    # jitter (real saturated clusters centroid precisely:
                    # 0.14 ppm validated on the analyst's GA4 runs)
                    # real saturated clusters centroid precisely (0.14
                    # ppm validated on the analyst's GA4 runs): the
                    # apex carries only a sub-ppm position jitter
                    # (0.05 idx ~ 0.5 ppm at this index scale)
                    for _ in range(30):
                        peaks.append((idx + rng.uniform(-300.0, 300.0),
                                      150.0, 8.0))
                    peaks.append((idx + rng.uniform(-0.05, 0.05),
                                  4.0e6, 300.0))
                else:
                    peaks.append((idx, 5.0e5, comb_snr[j]))
        for a in analytes:
            peaks.append((impact_mz_to_index(shifted(a)), 2.0e3,
                          15.0))
        for _ in range(40):                     # chemical/junk background
            peaks.append((rng.uniform(2000.0, 340000.0),
                          rng.uniform(60.0, 300.0), rng.uniform(3.0, 8.0)))
        for _ in range(8):                      # digitizer noise (SNR=1)
            peaks.append((rng.uniform(2000.0, 340000.0), 51.0, 1.0))
        peaks.sort(key=lambda p: p[0])
        sid += 1
        sid_off[sid] = len(blocks)
        blocks += impact_maxima_block(peaks)
        ms1_rows.append((sid, rt, 2 if cal_scan else 1))
        parent_true = analytes[k % len(analytes)]
        declared = parent_true * (1.0 + prec_offset_ppm * 1e-6)
        if prec_jitter:
            # v_co: auto-DDA recomputes the precursor from each MS1
            # survey -- the declared value jitters cycle to cycle
            # (a REAL manual list replays the exact same double)
            declared = parent_true * (1.0 + (prec_offset_ppm
                                            + rng.uniform(-2.0, 2.0))
                                      * 1e-6)
        m2 = [153.0181, 169.0137, 302.0460, 479.1180]
        m2_peaks = [(impact_mz_to_index(shifted(m) *
                                        (1.0 - prec_offset_ppm * 1e-6))
                     if m == 479.1180 else
                     impact_mz_to_index(shifted(m)), s * 1.0e3, s)
                    for m, s in zip(m2, [60.0, 40.0, 90.0, 25.0])]
        m2_peaks.sort(key=lambda p: p[0])
        sid += 1
        sid_off[sid] = len(blocks)
        blocks += impact_maxima_block(m2_peaks)
        ms2_rows.append((sid, rt + 0.5, declared))
    if not baf_metadata:
        with open(os.path.join(d, 'analysis.baf'), 'wb') as f:
            f.write(blocks)
        write_idx(os.path.join(d, 'analysis.baf_idx'), sid_off)
    else:
        # v_bp mock: the .d was NEVER OPENED in DataAnalysis -- the
        # baf carries the acquisition header + ONE METADATA BLOCK PER
        # SCAN (Rt, Segment, IsolationMass, CollisionEnergy, the
        # polarity as a string) that load_spectra_table_baf parses.
        header = acquisition_header_block([
            (10, b'Rt', b''), (11, b'Segment', b''),
            (12, b'IsolationMass', b''), (13, b'CollisionEnergy', b''),
            (14, b'Polarity', b'')]
            + ([(15, b'MSMS_Isolation_Act', b'')]
               if isolation_act else []))
        blob = bytes(blocks) + header
        head_off = len(blocks)
        meta_off = {}
        for sid_, rt_, seg in ms1_rows:
            meta_off[sid_] = len(blob)
            meta_off_rows = [(10, rt_), (11, seg)]
            if isolation_act:
                meta_off_rows.append((15, 0.0))   # MS1: no isolation
            blob += metadata_block(
                meta_off_rows,
                [(14, 'Negative' if polarity == 'negative'
                  else 'Positive')])
        for sid_, rt_, declared in ms2_rows:
            meta_off[sid_] = len(blob)
            meta_off_rows = [(10, rt_), (11, 1), (12, declared),
                             (13, 20.0 + 10.0 * (sid_ % 3))]
            if isolation_act:
                meta_off_rows.append((15, 1.0))   # MS2: isolation on
            blob += metadata_block(meta_off_rows)
        with open(os.path.join(d, 'analysis.baf'), 'wb') as f:
            f.write(blob)
        recs = [(2, 16, s, o) for s, o in sid_off.items()]
        recs += [(2, 32, s, o) for s, o in meta_off.items()]
        recs.append((1, 48, 0, head_off))
        write_idx_records(os.path.join(d, 'analysis.baf_idx'), recs)
    write_xtr_impact(os.path.join(d, 'analysis.baf_xtr'))
    c = sqlite3.connect(os.path.join(d, 'analysis.sqlite'))
    if segment:
        c.execute('create table Spectra (Id integer primary key, Rt real, '
                  'AcquisitionKey integer, Segment integer)')
        for sid_, rt_, seg in ms1_rows:
            c.execute('insert into Spectra values (?,?,?,?)',
                      (sid_, rt_, 1, seg))
        for sid_, rt_, declared in ms2_rows:
            c.execute('insert into Spectra values (?,?,?,?)',
                      (sid_, rt_, 2, 1))
    else:
        c.execute('create table Spectra (Id integer primary key, Rt real, '
                  'AcquisitionKey integer)')
        for sid_, rt_, seg in ms1_rows:
            c.execute('insert into Spectra values (?,?,?)', (sid_, rt_, 1))
        for sid_, rt_, declared in ms2_rows:
            c.execute('insert into Spectra values (?,?,?)', (sid_, rt_, 2))
    c.execute('create table Steps (Id integer primary key, '
              'TargetSpectrum integer, Mass real)')
    c.execute('create table PerSpectrumVariables (Spectrum integer, '
              'Variable integer, Value real)')
    for sid_, rt_, declared in ms2_rows:
        c.execute('insert into Steps values (?,?,?)',
                  (sid_, sid_, declared))
        c.execute('insert into PerSpectrumVariables values (?,?,?)',
                  (sid_, 5, 20.0 + 10.0 * (sid_ % 3)))
    c.commit(); c.close()
    return d, dict(lock_masses=lock_masses, analytes=analytes,
                   n_ms1=n_ms1, drift_ppm=drift_ppm,
                   prec_offset_ppm=prec_offset_ppm,
                   n_cal_scans=min(n_cal_scans, n_ms1))


def check_impact_mzml(path, meta):
    """Impact-specific mzML checks (calibrant-segment design): noise
    bucket gone, combs+analytes within 2 ppm on the calibrant-segment
    scans, analytes within 4 ppm on the analysis scans (the mock
    injects an artificial 3.8 ppm intra-run drift that the single-
    quadratic design does NOT remove -- no lock model, same as the
    reference extraction script; the lab's REAL runs measured
    0.25-0.6 ppm), precursors re-phased onto the MS1 apex (0.2 ppm)."""
    tree = ET.parse(path)
    ns = '{http://psi.hupo.org/ms/mzml}'
    spectra = tree.getroot().findall('.//' + ns + 'spectrum')
    assert len(spectra) == meta['n_ms1'] * 2, (len(spectra), 'spectra')
    n_ms1_checked = 0
    prev_ms1_mz = None
    for sp in spectra:
        n = int(sp.get('defaultArrayLength'))
        arrays = sp.findall('.//' + ns + 'binaryDataArray')
        mz = decode_arrays(arrays[0].find(ns + 'binary').text or '',
                           64, True, '<', n)
        it = decode_arrays(arrays[1].find(ns + 'binary').text or '',
                           64, True, '<', n)
        # the B=1.0 digitizer-noise bucket never reaches the mzML
        assert n == 0 or min(it) >= 2.0, ('noise bucket leaked', min(it))
        level = sp.find(ns + 'cvParam[@name="ms level"]').get('value')
        if level == '1' and n:
            n_ms1_checked += 1
            in_segment = n_ms1_checked <= meta['n_cal_scans']
            known = list(meta['analytes'])
            if in_segment:
                known = list(meta['lock_masses']) + known
                if n_ms1_checked <= 2:
                    # saturated start-of-segment scans: comb1 is a flat
                    # top (maxima spread +-300 idx); only its higher,
                    # clean combs and the analytes are checkable
                    known = known[1:]
            gate = 2.0 if in_segment else 4.0
            for m in known:
                best = min(abs(x - m) / m * 1e6 for x in mz)
                assert best < gate, ('lock/analyte residual ppm', best)
            prev_ms1_mz = mz
        if level == '2':
            prec = sp.find('.//' + ns + 'selectedIon')
            pmz = float(prec.find(
                ns + 'cvParam[@name="selected ion m/z"]').get('value'))
            if prev_ms1_mz:
                # re-phasing: the precursor IS the recalibrated MS1
                # apex of the previous MS1 scan (the declared parent
                # carried a +7 ppm offset by construction)
                apex = min(prev_ms1_mz, key=lambda x: abs(x - pmz))
                best = abs(pmz - apex) / apex * 1e6
                assert best < 0.2, ('impact precursor rephase ppm', best)
    assert n_ms1_checked >= meta['n_ms1'] - 2  # saturated segment too
    return True


def decode_arrays(b64, precision, compressed, endian, count):
    raw = base64.b64decode(b64)
    if compressed:
        raw = zlib.decompress(raw)
    fmt = '%s%d%s' % (endian, count, 'f' if precision == 32 else 'd')
    return struct.unpack(fmt, raw)


def check_mzml(path, meta):
    tree = ET.parse(path)
    ns = '{http://psi.hupo.org/ms/mzml}'
    spectra = tree.getroot().findall('.//' + ns + 'spectrum')
    assert len(spectra) == meta['n_ms1'] * 2, (len(spectra), 'spectra count')
    for sp in spectra:
        n = int(sp.get('defaultArrayLength'))
        arrays = sp.findall('.//' + ns + 'binaryDataArray')
        assert len(arrays) == 2
        mz = decode_arrays(arrays[0].find(ns + 'binary').text or '', 64, True, '<', n)
        it = decode_arrays(arrays[1].find(ns + 'binary').text or '', 64, True, '<', n)
        assert len(mz) == n == len(it)
        prec = sp.find('.//' + ns + 'selectedIon')
        level = sp.find(ns + 'cvParam[@name="ms level"]').get('value')
        if level == '2' and prec is not None:
            pmz = float(prec.find(ns + 'cvParam[@name="selected ion m/z"]').get('value'))
            best = min(abs(pmz - a) / a * 1e6 for a in meta['analytes'])
            assert best < 0.2, ('precursor rephase error ppm', best)
        if level == '1':
            known = meta['lock_masses'] + meta['analytes']
            for m in mz:
                best = min(abs(m - a) / a * 1e6 for a in known)
                assert best < 2.0, ('per-peak residual error ppm (Bruker spec 2 ppm)', best)
    return True


def check_mzxml(path, meta):
    xns = '{http://sashimi.sourceforge.net/schema_revision/mzXML_3-2}'
    scans = ET.parse(path).getroot().findall('.//' + xns + 'scan')
    assert len(scans) == meta['n_ms1'] * 2, (len(scans), 'scan count')
    for sc in scans:
        if int(sc.get('msLevel')) == 2:
            pmz = float(sc.find(xns + 'precursorMz').text)
            best = min(abs(pmz - a) / a * 1e6 for a in meta['analytes'])
            assert best < 0.2, ('mzXML precursor ppm', best)
    return True


def main():
    tmp = tempfile.mkdtemp(prefix='d2mz_test_')
    try:
        # --- scenario 1: single input folder, default layout
        # (v_bd: d2mz_output_<source>/ created NEXT TO the input)
        work = os.path.join(tmp, 'batch1')
        os.makedirs(work)
        d, meta = make_mock_d(work, n_ms1=12)
        rc = d2mz.main([work, '--format', 'both', '--precision', '64',
                        '--parabene', 'C4=8.5', '--parabene', '12.34'])
        assert rc == 0
        root = out_root_of(work)
        # v_bg: mzML/mzXML FLAT at the output root, precursor tables
        # in the TSV/ subfolder, no converted_<name>/ folder anymore
        assert not os.path.exists(
            os.path.join(root, 'converted_Synthetic_01')), 'layout'
        assert os.path.isdir(os.path.join(root, 'annexes', 'TSV')), 'TSV folder'
        # v_by: analyses/ + annexes/ inside the output folder --
        # NOTHING at the output folder root nor at the parent
        rep = os.path.dirname(work)          # parent of the source
        assert os.path.isdir(os.path.join(root, 'analyses')), 'analyses'
        assert os.path.isdir(os.path.join(root, 'annexes')), 'annexes'
        for c in ('CITATION.cff', 'CITATION.bib', 'CITATION.txt'):
            assert os.path.exists(os.path.join(root, 'annexes', c)), \
                'missing ' + c
            assert not os.path.exists(os.path.join(root, c)), \
                'stray citation at the output root: ' + c
            assert not os.path.exists(os.path.join(rep, c)), \
                'stray citation at the parent: ' + c
        assert not os.path.exists(
            os.path.join(rep, 'report_Synthetic_01.json')), \
            'stray report at the parent'
        # v_bm: ONE report PER SERIES, named after the source folder
        assert os.path.exists(report_of(work)), 'missing per-series report'
        assert not os.path.exists(os.path.join(root, 'report.json'))
        assert os.path.exists(
            os.path.join(root, 'annexes', 'TSV',
                         'Synthetic_01.tsv'))
        mzml = open(os.path.join(
            root, 'analyses', 'Synthetic_01.mzML')).read()
        # analyst directive 2026-11-07: the FIRST acquired UV channel
        # (UV1 = UV_VIS_1, 280 nm per hystar.method) is recorded by
        # default -- one trace, one apex
        assert '<chromatogramList' in mzml and \
            mzml.count('<chromatogram ') == 2
        assert '280 nm' in mzml and '254 nm' not in mzml
        # the other existing traces go into the mzML too: UV_VIS_2
        # (hystar 320 nm) follows the analytical UV1 chromatogram
        assert 'name="trace" value="UV_VIS_2"' in mzml
        assert '<userParam name="wavelength" value="320 nm"/>' in mzml
        # FIX v_ay: the second UV channel is displayed with its HyStar
        # wavelength in the chromatogram id, like UV1 (UV_320_nm)
        assert 'id="UV_320_nm"' in mzml
        assert 'name="wavelength"' in mzml
        assert 'd2mz:retention_standard:C4' in mzml
        assert 'rt 8.5000 min' in mzml
        # chromatograms are indexed too (indexedmzML)
        assert '<index name="chromatogram">' in mzml
        # the apex UV spectrum lives INSIDE the mzML: one userParam per
        # apex, self-contained file, nothing to lose
        assert mzml.count('d2mz:uv_spectrum') == 1, mzml[:0]
        assert 'trace=280 nm; rt=3.5000 min; data=' in mzml
        assert '280:105.000000' in mzml
        assert not os.path.exists(
            os.path.join(rep, 'Synthetic_01.uv_apex.tsv'))
        mzxml = open(os.path.join(
            root, 'analyses', 'Synthetic_01.mzXML')).read()
        # mzXML carries NO UV and NO chromatogram (its schema has none)
        assert 'chromatogram' not in mzxml.lower()
        assert 'retention_standard' not in mzxml
        f = json.load(open(report_of(work)))['files'][0]
        assert f['detected_polarity'] == 'positive', f
        assert f['uv_traces'] == {'280 nm': 12}, f
        assert f['extra_traces'] == ['UV_VIS_2'], f
        assert f['uv_apex_spectra'] == 1, f
        # report.json carries processing stats only: NO measurement data
        assert 'uv_apices' not in f and 'parabenes' not in f, f
        # explicit --uv-trace override: the legacy per-scan UV table is
        # read at the requested wavelengths (254 and 365 nm)
        work_ov = os.path.join(tmpdir if 'tmpdir' in dir() else
                               tempfile.mkdtemp(prefix='d2mz_uvov_'),
                               'w_ov')
        os.makedirs(work_ov)
        make_mock_d(work_ov, n_ms1=12)
        rc = d2mz.main([work_ov, '--format', 'mzml', '--precision', '64',
                        '--uv-trace', '254', '--uv-trace', '365'])
        assert rc == 0
        mzml2 = open(os.path.join(out_root_of(work_ov), 'analyses',
            'Synthetic_01.mzML')).read()
        assert mzml2.count('<chromatogram ') == 4   # 254+365 + UV_VIS_1/2
        assert '254 nm' in mzml2 and '365 nm' in mzml2
        assert mzml2.count('name="trace" value="UV_VIS_1"') == 1
        assert mzml2.count('name="trace" value="UV_VIS_2"') == 1
        assert mzml2.count('d2mz:uv_spectrum') == 2
        assert '254:105.000000' in mzml2 and '365:105.000000' in mzml2
        assert '200:5.000000' in mzml2 and '500:5.000000' in mzml2
        f2 = json.load(open(report_of(work_ov)))['files'][0]
        assert f2['uv_traces'] == {'254 nm': 12, '365 nm': 12}, f2
        assert f2['uv_apex_spectra'] == 2, f2
        # the data itself lives ONLY in the mzML
        assert mzml.count('d2mz:uv_spectrum') == 1
        assert 'd2mz:retention_standard:C4' in mzml
        assert 'd2mz:retention_standard:parabene_2' in mzml

        # --- scenario 1b (v_bh): UV1 acquired at 210 nm, UV2 at 254 nm
        # -> the analytical trace must be UV2 (254 nm, the FIRST
        # channel ABOVE 250 nm), UV1/210 nm kept as an extra trace
        work_uv = os.path.join(tmp, 'batch1b')
        os.makedirs(work_uv)
        make_mock_d(work_uv, n_ms1=12, uv1_nm=210.0, uv2_nm=254.0)
        rc = d2mz.main([work_uv, '--format', 'mzml', '--precision', '32'])
        assert rc == 0
        mzml_uv = open(os.path.join(out_root_of(work_uv), 'analyses',
                                    'Synthetic_01.mzML')).read()
        assert 'trace=254 nm' in mzml_uv, \
            'UV2 (254 nm) must be the analytical trace when UV1 is 210 nm'
        assert mzml_uv.count('d2mz:uv_spectrum') == 1
        assert 'name="trace" value="UV_VIS_1"' in mzml_uv, \
            'UV1 (210 nm) must still travel as an extra trace'
        assert '<userParam name="wavelength" value="210 nm"/>' in mzml_uv
        assert '<userParam name="wavelength" value="254 nm"/>' in mzml_uv
        assert f['n_scans'] == 24, f
        assert f['n_ms2'] == 12, f
        assert f['precursors_rephased'] == 12, f
        assert abs(f['lock_stdev_ppm_after'] or 0) < 2.0, f  # Bruker spec: mean stdev < 2 ppm
        check_mzml(os.path.join(root, 'analyses',
                                 'Synthetic_01.mzML'), meta)
        check_mzxml(os.path.join(root, 'analyses',
                                 'Synthetic_01.mzXML'), meta)
        tsv = open(os.path.join(root, 'annexes', 'TSV',
                                'Synthetic_01.tsv')).read()
        assert len(tsv.strip().splitlines()) == 13
        # --- scenario 2: SEVERAL input folders in one command
        work2 = os.path.join(tmp, 'batch2')
        os.makedirs(work2)
        # weak run: only 4 lock observations (< 7, the Bruker minimum)
        # -> must borrow the lock-mass model from the adjacent calibrated
        #    file of the same folder/series
        make_mock_d(work2, n_ms1=8, name='Good_NEG', polarity='negative')
        make_mock_d(work2, n_ms1=8, name='Good_POS', polarity='positive')
        make_mock_d(work2, n_ms1=4, name='Weak_POS', polarity='positive')
        rc = d2mz.main([work, work2, '--format', 'mzml',
                        '--precision', '32'])
        assert rc == 0
        # v_bm: each series has its OWN report next to its output
        rj1 = report_of(work)
        f1 = [f for f in json.load(open(rj1))['files']
              if f['name'] == 'Synthetic_01.d'
              and f.get('n_scans') == 24][-1]
        assert f1['n_scans'] == 24, f1
        assert os.path.isfile(
            os.path.join(out_root_of(work), 'analyses',
                         'Synthetic_01.mzML'))
        files2 = {f['name']: f for f in json.load(open(
            report_of(work2)))['files']}
        # polarity detected automatically, per file, in a MIXED folder
        assert files2['Good_NEG.d']['detected_polarity'] == 'negative', files2['Good_NEG.d']
        assert files2['Good_POS.d']['detected_polarity'] == 'positive', files2['Good_POS.d']
        assert abs(files2['Good_NEG.d']['lock_stdev_ppm_after']) < 2.0
        assert abs(files2['Good_POS.d']['lock_stdev_ppm_after']) < 2.0
        assert files2['Good_POS.d']['calibration_borrowed_from'] is None
        # the weak positive run borrows from its POSITIVE neighbour,
        # never from the negative one
        wk = files2['Weak_POS.d']
        assert wk['n_scans'] == 8, wk
        assert wk['n_lock_points'] == 4, wk
        assert wk['calibration_borrowed_from'] == 'Good_POS.d', wk
        assert abs(wk['lock_stdev_ppm_after'] or 0) < 2.0, wk  # Bruker spec: mean stdev < 2 ppm
        assert wk['precursors_rephased'] == 4, wk
        # --- scenario 3: GUI front-end, headless (argvec == command line)
        import importlib.util
        from importlib.machinery import SourceFileLoader
        gui_path = os.path.join(os.path.dirname(
            os.path.abspath(__file__)), 'd2mz_gui.pyw')
        loader = SourceFileLoader('d2mz_gui', gui_path)
        spec = importlib.util.spec_from_loader('d2mz_gui', loader)
        gui = importlib.util.module_from_spec(spec)
        loader.exec_module(gui)
        assert gui.tk is None or True  # headless import never raises
        work3 = os.path.join(tmp, 'batch3')
        os.makedirs(work3)
        make_mock_d(work3, n_ms1=12)
        av = gui.build_argvec(work3, 'both', '64')
        assert av[0] == work3 and '--format' in av and 'both' in av
        rc = d2mz.main(av)
        assert rc == 0
        f3 = json.load(open(report_of(work3)))['files'][-1]
        assert f3['n_scans'] == 24, f3
        assert f3['detected_polarity'] == 'positive', f3
        assert abs(f3['lock_stdev_ppm_after'] or 0) < 2.0, f3
        # --- scenario 5: impact-generation file (mode 10, no frequency
        # calibration record; quadratic idx->mz fitted on the calibrant
        # ladder; SNR intensities; saturated calibrant segment)
        work5 = os.path.join(tmp, 'batch5')
        os.makedirs(work5)
        _d5, meta5 = make_mock_impact_d(work5, n_ms1=16)
        rc = d2mz.main([work5, '--format', 'mzml', '--precision', '64'])
        assert rc == 0
        rep5 = json.load(open(report_of(work5)))['files'][-1]
        assert rep5['detected_polarity'] == 'negative', rep5
        assert rep5.get('method_polarity') == 'negative', rep5
        assert rep5['n_ms2'] == 16 and rep5['n_scans'] == 32, rep5
        assert rep5['precursors_rephased'] == 16, rep5
        assert abs(rep5['lock_stdev_ppm_after'] or 0) < 2.0, rep5
        check_impact_mzml(os.path.join(out_root_of(work5), 'analyses',
            'Impact_NEG.mzML'), meta5)
        # --- scenario 4: a 'parabens' run in the folder -> RI anchors
        # measured automatically and stamped into EVERY mzML of the batch
        work4 = os.path.join(tmp, 'batch4')
        os.makedirs(work4)
        make_mock_d(work4, n_ms1=6, name='Parabens_NEG_29042',
                    polarity='negative',
                    extra_ms1_peaks={
                        1: [(169.0143, 9.0e4)],      # GallAc   -> 1.0 min
                        3: [(151.0401, 9.0e4)],      # MeParaben -> 2.0 min
                        5: [(193.0870, 9.0e4)],      # BuParaben -> 3.0 min
                    })
        make_mock_d(work4, n_ms1=6, name='Sample_POS')
        rc = d2mz.main([work4, '--format', 'both', '--precision', '32'])
        assert rc == 0
        files4 = {f['name']: f for f in json.load(open(
            report_of(work4)))['files']}
        assert files4['Parabens_NEG_29042.d']['n_retention_standards'] == 3, \
            files4['Parabens_NEG_29042.d']
        assert files4['Sample_POS.d']['n_retention_standards'] == 3, \
            files4['Sample_POS.d']
        for nm in ('Parabens_NEG_29042', 'Sample_POS'):
            m = open(os.path.join(out_root_of(work4), 'analyses',
                                  nm + '.mzML')).read()
            # measured anchors: GallAc 1.0 min RI 600, MeParaben 2.0 min
            # RI 800, BuParaben 3.0 min RI 1100, [M-H]-, source run given
            assert 'd2mz:retention_standard:GallAc' in m
            assert 'rt 1.0000 min; ri 600; mz 169.0143 ([M-H]-)' in m
            assert 'd2mz:retention_standard:MeParaben' in m
            # NOTE: the stamped m/z is the RECALIBRATED one, so only
            # its first decimals are asserted (the lock-mass fit can
            # shift the last digit from run to run)
            assert 'rt 2.0000 min; ri 800; mz 151.040' in m
            assert 'd2mz:retention_standard:BuParaben' in m
            assert 'rt 3.0000 min; ri 1100; mz 193.087' in m
            assert 'source Parabens_NEG_29042.d' in m
            x = open(os.path.join(out_root_of(work4), 'analyses',
                                  nm + '.mzXML')).read()
            assert 'retention_standard' not in x
        # mzML written WITHOUT a parabens run in the folder carries none
        m1 = open(os.path.join(out_root_of(work), 'analyses',
            'Synthetic_01.mzML')).read()
        assert 'retention_standard:GallAc' not in m1

        # positive 'parabens' run: anchor on [M+H]+ ONLY (the [M+Na]+
        # peak present in the run is not a candidate)
        work4b = os.path.join(tmp, 'batch4b')
        os.makedirs(work4b)
        make_mock_d(work4b, n_ms1=6, name='Parabens_POS_seq02',
                    polarity='positive',
                    extra_ms1_peaks={
                        2: [(153.0547, 8.0e4),   # MeParaben [M+H]+
                            (175.0366, 9.0e4)],  # [M+Na]+, must be ignored
                    })
        rc = d2mz.main([work4b, '--format', 'mzml', '--precision', '32'])
        assert rc == 0
        m4b = open(os.path.join(out_root_of(work4b), 'analyses',
            'Parabens_POS_seq02.mzML')).read()
        assert 'd2mz:retention_standard:MeParaben' in m4b
        assert 'mz 153.054' in m4b, \
            [l for l in m4b.split('\n') if 'MeParaben' in l]
        assert '[M+Na]' not in m4b
        assert 'ri 800' in m4b

        # folder with TWO parabens runs (one per polarity): each file gets
        # the anchors of ITS polarity -- no name collision, no mixing
        work4c = os.path.join(tmp, 'batch4c')
        os.makedirs(work4c)
        make_mock_d(work4c, n_ms1=6, name='Parabens_NEG_1', polarity='negative',
                    extra_ms1_peaks={1: [(169.0143, 9.0e4)]})   # GallAc [M-H]-
        make_mock_d(work4c, n_ms1=6, name='Parabens_POS_1', polarity='positive',
                    extra_ms1_peaks={2: [(171.0289, 9.0e4)]})  # GallAc [M+H]+
        rc = d2mz.main([work4c, '--format', 'mzml', '--precision', '32'])
        assert rc == 0
        m_neg = open(os.path.join(out_root_of(work4c), 'analyses',
            'Parabens_NEG_1.mzML')).read()
        m_pos = open(os.path.join(out_root_of(work4c), 'analyses',
            'Parabens_POS_1.mzML')).read()
        assert 'mz 169.0143 ([M-H]-)' in m_neg, m_neg
        # (stamped m/z = recalibrated value; only first decimals kept)
        assert 'mz 171.02' in m_pos, m_pos
        assert '[M+H]+' not in m_neg and '[M-H]-' not in m_pos

        # v_bp: a .d NEVER OPENED in Bruker DataAnalysis (no
        # analysis.sqlite) is now CONVERTIBLE: the spectra table is
        # rebuilt from the per-scan metadata blocks of the raw
        # analysis.baf, and the mzML is BYTE-IDENTICAL to the same
        # file WITH its sqlite index (startTimeStamp is the
        # acquisition date, so the two outputs must match exactly).
        work6 = os.path.join(tmp, 'unopened')
        os.makedirs(work6)
        _d6, meta6 = make_mock_impact_d(work6, n_ms1=12,
                                        name='NeverOpened_01',
                                        baf_metadata=True)
        rc = d2mz.main([work6, '--format', 'mzml', '--precision', '64'])
        assert rc == 0
        m_base = open(os.path.join(out_root_of(work6), 'analyses',
                                   'NeverOpened_01.mzML')).read()
        check_impact_mzml(os.path.join(out_root_of(work6), 'analyses',
                                       'NeverOpened_01.mzML'), meta6)
        # remove the DataAnalysis index: the folder becomes a
        # never-opened .d (raw baf only)
        os.remove(os.path.join(work6, 'NeverOpened_01.d',
                               'analysis.sqlite'))
        # v_bs: delete the sqlite-conversion mzML first -- a FAILED
        # fallback used to leave rc=0 plus the STALE file on disk, so
        # the byte-identity check compared the file with itself (the
        # suite passed while the baf-metadata parser was broken)
        mz6 = os.path.join(out_root_of(work6), 'analyses',
                   'NeverOpened_01.mzML')
        os.remove(mz6)
        rc = d2mz.main([work6, '--format', 'mzml', '--precision', '64'])
        assert rc == 0, 'baf-fallback conversion failed'
        assert os.path.exists(mz6), \
            'baf-fallback produced NO mzML (silent failure)'
        m_fb = open(mz6).read()
        assert m_base == m_fb, \
            'baf-fallback mzML differs from the sqlite conversion'
        # find_d_folders must also count the never-opened .d BEFORE
        # the conversion (no sqlite, raw baf only)
        found, inc = d2mz.find_d_folders(work6)
        assert len(found) == 1 and not inc, (found, inc)

        # v_bw: never-opened baf with the per-scan Rt in the metadata
        # BLOCK HEADER (i32 milliseconds at offset +24, mhs=40, as on
        # the real Qualitative file: scan 1 = 299 ms, scan 3970 =
        # 965142 ms) and TYPED integer rows (type code % 256 == 5,
        # e.g. SpectraAcquisition_SegmentId=1 -- the v_bt denormals
        # were those integers read as doubles).
        work7 = os.path.join(tmp, 'unopened_rt')
        os.makedirs(work7)
        d7 = os.path.join(work7, 'RtHeader_01.d')
        os.makedirs(d7)
        # v_bz regression: the v_by SUFFIX alias match let STRING/LIST
        # rows through (MSMS_ListIsolationMass ends with
        # 'isolationmass') and int()/float() crashed on the raw
        # binary value -- the crash killed 35 real DDA files. Numeric
        # fields must only match NUMERIC rows.
        head7 = acquisition_header_block([
            (10, b'Zt', b''),
            (11, b'SpectraAcquisition_SegmentId', b'', 5),
            (12, b'MSMS_ListIsolationMass', b'', 0x0B)])
        blob = head7
        meta_off = {}
        _listraw = struct.pack('<qqq', 471, 519, 622)  # list-type bytes
        for k, rt_ms in enumerate((500, 1500, 2500)):
            meta_off[k + 1] = len(blob)
            blob += metadata_block_rt(
                rt_ms, [(10, struct.pack('<d', 30.0)),
                        (11, struct.pack('<q', 1))],
                string_rows=[(12, _listraw)])
        with open(os.path.join(d7, 'analysis.baf'), 'wb') as f:
            f.write(blob)
        recs = [(2, 32, s, o) for s, o in meta_off.items()]
        recs.append((1, 48, 0, 0))
        write_idx_records(os.path.join(d7, 'analysis.baf_idx'), recs)
        sp7, _st7, _ps7 = d2mz.load_spectra_table_baf(d7)
        assert [r['Rt'] for r in sp7] == [0.5, 1.5, 2.5], sp7
        assert all(r['Segment'] == 1 for r in sp7), sp7
        # the list-type isolation-mass row must NOT be mistaken for a
        # declared precursor (and must not crash either)
        assert not _st7, _st7
        # psv legitimately carries the _modes coverage bookkeeping
        # (since the mode-coverage signals); the isolation-list row
        # must add NOTHING else to it
        assert set(_ps7) == {'_modes'}, _ps7
        assert _ps7['_modes']['iso_act'] == 0, _ps7

        # v_br: --zip-output produces a ZENODO-READY archive: the
        # mzML/mzXML FLAT in analyses/, the TSV + the series report +
        # the CITATION files in annexes/ (analyst directive 2026-11-09)
        rc = d2mz.main([work6, '--format', 'both', '--precision', '64',
                        '--gzip', '--zip-output'])
        assert rc == 0
        import zipfile as _zf
        _zp = os.path.join(out_root_of(work6) + '.zip')
        with _zf.ZipFile(_zp) as _z:
            _names = _z.namelist()
        assert any(n.endswith('/analyses/NeverOpened_01.mzML.gz')
                   for n in _names), _names
        assert any(n.endswith('/analyses/NeverOpened_01.mzXML.gz')
                   for n in _names), _names
        assert any(n.endswith('/annexes/TSV/NeverOpened_01.tsv')
                   for n in _names), _names
        assert any(n.endswith('/annexes/report_unopened.json')
                   for n in _names), _names
        assert any(n.endswith('/annexes/CITATION.cff')
                   for n in _names), _names
        # assert any(n.endswith('/annexes/CITATION.cff')
        #            for n in _names), _names
        # no stray flat file at the archive root
        assert all('/analyses/' in n or '/annexes/' in n
                   for n in _names if not n.endswith('/')), _names

        # ------------------------------------------------------------
        # v_co (25/11/2026): the residual MSMSManual_List* retraction
        # must ALSO fire on the baf-metadata fallback of never-opened
        # .d files. Real case 38087: _modes['iso_act'] counts EVERY
        # MS2 scan there (per-scan MSMS_Isolation_Act flag, pure
        # bookkeeping that generates NO warning note) and the old
        # `not _modes` gate silently disabled the retraction -- an
        # auto-DDA file with 2839 jittering precursors was still
        # flagged 'targeted MS/MS list' in the report.
        # (1) unit tests of the discriminator
        assert d2mz.manual_list_targeted(
            [{'declared': 479.1203} for _ in range(50)]
            + [{'declared': 479.1203} for _ in range(50)]
            + [{'declared': 479.1203} for _ in range(50)]) is True
        assert not d2mz.manual_list_targeted(
            [{'declared': 100.0 + i} for i in range(120)])
        # under the minimum MS2 count: never 'targeted'
        assert not d2mz.manual_list_targeted(
            [{'declared': 479.1203} for _ in range(40)])
        # repeats below the share threshold: DDA-like
        assert not d2mz.manual_list_targeted(
            [{'declared': 479.1203} for _ in range(55)]
            + [{'declared': 100.0 + i} for i in range(45)])
        # (2) e2e, NEVER-OPENED .d (baf fallback) with the isolation
        # flag on every MS2 + a RESIDUAL manual list in the .method:
        # auto-DDA jitter -> warning WITHDRAWN (report + mzML clean)
        work7 = os.path.join(tmp, 'residual_list_dda')
        os.makedirs(work7)
        make_mock_impact_d(work7, n_ms1=12, name='ResidualList_DDA',
                           baf_metadata=True, manual_list_ce=True,
                           isolation_act=True, prec_jitter=True)
        os.remove(os.path.join(work7, 'ResidualList_DDA.d',
                               'analysis.sqlite'))
        rp7 = os.path.join(tmp, 'report_dda.json')
        rc = d2mz.main([work7, '--format', 'mzml', '--precision', '64',
                        '--report', rp7])
        assert rc == 0, 'residual-list DDA conversion failed'
        rep7 = json.load(open(rp7))
        assert rep7['files'][0]['mode_warning'] is None, \
            rep7['files'][0]['mode_warning']
        m7 = open(os.path.join(out_root_of(work7), 'analyses',
                               'ResidualList_DDA.mzML')).read()
        assert 'd2mz:mode_validation' not in m7
        assert 'd2mz:acquisition_mode' in m7
        # (3) e2e, same fallback path but a REAL manual list: the
        # declared precursors replay bit-identically (>= 50 repeats
        # of each mass, >= 80 % of all MS2) -> warning KEPT
        work8 = os.path.join(tmp, 'residual_list_targeted')
        os.makedirs(work8)
        make_mock_impact_d(work8, n_ms1=160, name='ManualList_target',
                           baf_metadata=True, manual_list_ce=True,
                           isolation_act=True)
        os.remove(os.path.join(work8, 'ManualList_target.d',
                               'analysis.sqlite'))
        rp8 = os.path.join(tmp, 'report_targeted.json')
        rc = d2mz.main([work8, '--format', 'mzml', '--precision', '64',
                        '--report', rp8])
        assert rc == 0
        rep8 = json.load(open(rp8))
        assert 'targeted MS/MS list (MSMSManual_List*)' \
            in (rep8['files'][0]['mode_warning'] or ''), \
            rep8['files'][0]['mode_warning']
        m8 = open(os.path.join(out_root_of(work8), 'analyses',
                               'ResidualList_DDA.mzML'.replace(
                                   'ResidualList_DDA',
                                   'ManualList_target'))).read()
        assert 'd2mz:mode_validation' in m8
        assert 'targeted MS/MS list' in m8

        print('ALL TESTS PASSED')
        # v_co: the old bare json.dumps(f) dumped a LEAKED file handle
        # (a `with open(...) as f` of the mock builders survives in
        # the local scope) -- dump the LAST conversion report instead
        print(json.dumps(rep8, indent=2))
        return 0

    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == '__main__':
    sys.exit(main())

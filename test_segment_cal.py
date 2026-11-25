#!/usr/bin/env python3
"""Dedicated test for the CALIBRANT-SEGMENT ladder fit (d2mz >= 2026-11-07_q):
the quad calibration must be fitted on the Segment==2 (calibrant segment)
scans only -- same logic as the reference extraction script
(ms_run_extract.calibrate, flag7==51 / segment 2, analyst directive
2026-10-31) -- with a fallback to all MS1 scans when no segment-2 scan
carries a usable ladder.

Builds the impact mock (same builder as the main harness), marks the two
saturated start-of-run MS1 scans as Segment 2 in analysis.sqlite, converts,
and checks:
  1. the log announces the calibrant segment (2 of 16 MS1 scans);
  2. the calibration still fits (lock model, drift removed);
  3. without the Segment column the fallback path still converts.
"""
import io
import contextlib
import json
import os
import shutil
import sqlite3
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import d2mz
import test_synthetic as T


def convert_and_log(d, root):
    import io
    import contextlib
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        rc = d2mz.main([d, '--out', root, '--format', 'mzml'])
    # v_co: the return code is REPORTED, not asserted here -- case 2
    # deliberately CONVERTS NOTHING (no Segment column, no same-mode
    # donor -> deferred file NOT converted, rc != 0); only the
    # callers that expect a successful conversion assert rc == 0.
    return rc, buf.getvalue()


def run_case(tmp, name, mark_segment):
    root = os.path.join(tmp, name)
    os.makedirs(root)
    # segment=True: the mock carries the real Impact II structure
    # (Spectra.Segment column, calibrant combs in the 9 Segment==2
    # scans only); segment=False: no Segment column -> deferred path
    d, meta = T.make_mock_impact_d(root, n_ms1=16, segment=mark_segment)
    rc, log = convert_and_log(d, root)
    mzml = os.path.join(root, 'analyses',
                        os.path.basename(d)[:-2] + '.mzML')
    return rc, log, mzml, meta


def main():
    tmp = tempfile.mkdtemp(prefix='d2mz_segcal_')
    try:
        # case 1: Segment column present, calibrant scans marked
        rc1, log, mzml, meta = run_case(tmp, 'seg', mark_segment=True)
        assert rc1 == 0, 'd2mz.main failed (case 1)'
        assert 'calibrant segment (Segment==2): 9 scans (of 16 MS1)' in log, \
            'calibrant-segment log line missing:\n' + log
        assert os.path.exists(mzml), 'mzML not written (case 1)'
        assert 'using the built-in exact negative list (51 masses)' \
            in log, 'built-in exact-list fit on the calibrant segment ' \
            'missing:\n' + log
        assert 'impact quadratic calibration (negative series)' in log, log
        assert 'falling back to all' not in log, \
            'unexpected fallback with a valid calibrant segment:\n' + log
        with open(mzml) as f:
            x = f.read()
        assert 'd2mz:lock' in x or 'referenceableParamGroup' in x
        print('case 1 (segment-2 marked): OK')

        # case 2: no Segment column -> file DEFERRED, no donor in the
        # folder -> NOT converted (analyst directive: 1. segment 2,
        # 2. another file -- NO other fit path)
        rc2, log2, mzml2, _ = run_case(tmp, 'noseg', mark_segment=False)
        # the lone deferred file converts NOTHING: rc != 0 is the
        # EXPECTED outcome (the error is the product, not a failure
        # of the test)
        assert rc2 != 0, \
            'a deferred file without a donor must NOT convert cleanly'
        assert 'calibrant absent' in log2, log2
        assert 'DEFERRED' in log2, log2
        assert 'no calibrant segment and no same-mode donor -> ' \
            'NOT converted' in log2, log2
        assert 'legacy fit' not in log2, \
            'the legacy all-MS1 fit was removed (analyst directive):' + log2
        assert not os.path.exists(mzml2), \
            'a file without segment 2 and without donor was converted'
        # v_bm: the report is PER SERIES, named report_<source>.json
        # (source = the .d here); the .d name is the mzML basename
        rep2 = json.load(open(os.path.join(
            tmp, 'noseg',
            'annexes',
            'report_%s.json' % os.path.basename(mzml2)[:-5])))['files']
        assert rep2 and rep2[0].get('error'), rep2
        print('case 2 (no Segment column, no donor -> NOT converted, '
              'reported): OK')

        # case 3: deferred file WITH a same-mode donor in the folder
        # -> its calibration is BORROWED from the donor (analyst
        # directive 2026-11-07)
        root3 = os.path.join(tmp, 'borrow')
        os.makedirs(root3)
        da, _ = T.make_mock_impact_d(root3, n_ms1=16, name='Impact_A',
                                     segment=True)
        db, meta_b = T.make_mock_impact_d(root3, n_ms1=16, name='Impact_B',
                                          segment=False)
        # donor A: Segment==2 built in; B has NO Segment column ->
        # calibrant absent -> deferred
        buf3 = io.StringIO()
        with contextlib.redirect_stdout(buf3):
            rc = d2mz.main([root3, '--out', root3, '--format', 'mzml',
                       '--precision', '64'])
        assert rc == 0, 'd2mz.main failed (case 3)'
        log3 = buf3.getvalue()
        assert 'calibration borrowed from Impact_A.d' in log3, log3
        by = {f_['name']: f_ for f_ in json.load(open(os.path.join(
            root3, 'annexes', 'report_borrow.json')))['files']}
        assert by['Impact_B.d']['calibration_borrowed_from'] == \
            'Impact_A.d', by.get('Impact_B.d')
        assert by['Impact_B.d']['detected_polarity'] == 'negative', \
            by['Impact_B.d']
        mbz = os.path.join(root3, 'analyses', 'Impact_B.mzML')
        assert os.path.exists(mbz)
        # borrowed calibration must deliver sound masses on the
        # deferred file. The mock artificially injects a 3.8 ppm
        # intra-run drift that the single-quad design does NOT remove
        # (no lock model, same as the reference extraction script,
        # validated at 0.3-0.6 ppm on the lab's REAL runs): the gate
        # here is therefore the injected drift bound, not Bruker's
        # 2 ppm spec which applies to real data.
        import xml.etree.ElementTree as ET
        ns = '{http://psi.hupo.org/ms/mzml}'
        tree = ET.parse(mbz)
        n_ms1_seen = 0
        for sp in tree.getroot().findall('.//' + ns + 'spectrum'):
            n = int(sp.get('defaultArrayLength'))
            if not n:
                continue
            level = sp.find(ns + 'cvParam[@name="ms level"]').get('value')
            if level != '1':
                continue
            n_ms1_seen += 1
            known = list(meta_b['analytes'])
            if n_ms1_seen <= meta_b['n_cal_scans']:
                # calibrant-segment scans carry the combs (first comb
                # spread +-300 idx on the saturated start-of-segment
                # scans, same exclusion as check_impact_mzml)
                known = list(meta_b['lock_masses']) + known
                if n_ms1_seen <= 2:
                    known = known[1:]
            arrays = sp.findall('.//' + ns + 'binaryDataArray')
            mz = T.decode_arrays(arrays[0].find(ns + 'binary').text,
                                 64, True, '<', n)
            for m in known:
                best = min(abs(x - m) / m * 1e6 for x in mz)
                assert best < 4.0, ('borrowed-cal residual ppm', best)
        print('case 3 (deferred -> calibration borrowed from '
              'same-mode donor, masses < 2 ppm): OK')
        # case 4: a 'parabens' run in the folder is CONVERTED FIRST and
        # its retention anchors are stamped into every mzML that
        # follows (analyst directive 2026-11-07)
        root4 = os.path.join(tmp, 'parabens')
        os.makedirs(root4)
        T.make_mock_d(root4, n_ms1=6, name='Sample_POS')
        T.make_mock_d(root4, n_ms1=6, name='Parabens_NEG_29042',
                      polarity='negative',
                      extra_ms1_peaks={
                          1: [(169.0143, 9.0e4)],      # GallAc  -> 1.0 min
                          3: [(151.0401, 9.0e4)],      # MeParaben -> 2.0
                          5: [(193.0870, 9.0e4)],      # BuParaben -> 3.0
                      })
        buf4 = io.StringIO()
        with contextlib.redirect_stdout(buf4):
            rc = d2mz.main([root4, '--out', root4, '--format', 'mzml'])
        assert rc == 0, 'd2mz.main failed (case 4)'
        log4 = buf4.getvalue()
        conv_lines = [ln for ln in log4.splitlines()
                      if '] converting ' in ln]
        assert len(conv_lines) == 2, conv_lines
        assert 'Parabens_NEG_29042.d' in conv_lines[0], \
            ('parabens run not converted first', conv_lines)
        assert 'Sample_POS.d' in conv_lines[1], conv_lines
        assert 'now averaged over 1 run(s)' in log4, log4
        rep4 = {f_['name']: f_ for f_ in json.load(open(os.path.join(
            root4, 'annexes', 'report_parabens.json')))['files']}
        assert rep4['Parabens_NEG_29042.d']['n_retention_standards'] == 3, \
            rep4['Parabens_NEG_29042.d']
        # anchors also stamped into the analysis file converted after
        assert rep4['Sample_POS.d']['n_retention_standards'] == 3, \
            rep4['Sample_POS.d']
        mz4 = open(os.path.join(root4, 'analyses',
                                'Sample_POS.mzML')).read()
        assert 'd2mz:retention_standard:' in mz4, \
            'retention anchors missing from the analysis mzML'
        print('case 4 (parabens converted first, anchors stamped into '
              'the analysis files): OK')
        # case 5: SEVERAL parabens runs of the same mode -> all are
        # processed and the AVERAGE per anchor is stamped (analyst
        # directive 2026-11-07)
        root5 = os.path.join(tmp, 'parabens_multi')
        os.makedirs(root5)
        T.make_mock_d(root5, n_ms1=6, name='Parabens_NEG_1',
                      polarity='negative',
                      extra_ms1_peaks={1: [(169.0143, 9.0e4)]})
        T.make_mock_d(root5, n_ms1=6, name='Parabens_NEG_2',
                      polarity='negative',
                      extra_ms1_peaks={3: [(169.0143, 9.0e4)]})
        T.make_mock_d(root5, n_ms1=6, name='Sample_NEG')
        buf5 = io.StringIO()
        with contextlib.redirect_stdout(buf5):
            rc = d2mz.main([root5, '--out', root5, '--format', 'mzml'])
        assert rc == 0, 'd2mz.main failed (case 5)'
        log5 = buf5.getvalue()
        conv5 = [ln for ln in log5.splitlines() if '] converting ' in ln]
        assert len(conv5) == 3, conv5
        assert 'Parabens_NEG_1.d' in conv5[0] and \
            'Parabens_NEG_2.d' in conv5[1] and \
            'Sample_NEG.d' in conv5[2], conv5
        assert 'now averaged over 2 run(s)' in log5, log5
        mz5 = open(os.path.join(root5, 'analyses',
                                'Sample_NEG.mzML')).read()
        # GallAc at 1.0 min in run 1, 2.0 min in run 2 -> mean 1.5 min
        assert 'd2mz:retention_standard:GallAc' in mz5, \
            'averaged anchor missing from the analysis mzML'
        assert 'rt 1.5000 min' in mz5, \
            'anchor rt is not the mean of the two parabens runs'
        assert 'mean of Parabens_NEG_1.d, Parabens_NEG_2.d' in mz5, \
            'source does not say the anchor is a mean'
        rep5 = json.load(open(os.path.join(
            root5, 'annexes', 'report_parabens_multi.json')))
        rs5 = rep5['retention_standards']['negative']
        ga5 = [a for a in rs5 if a['name'] == 'GallAc'][0]
        assert ga5['n_runs'] == 2, ga5
        assert ga5['rt_mean_min'] == 1.5, ga5
        assert ga5['rt_min'] == 1.0 and ga5['rt_max'] == 2.0, ga5
        assert ga5['rt_spread_s'] == 60.0, ga5      # 1 min spread
        assert 'Parabens_NEG_1.d' in ga5['sources'] and \
            'Parabens_NEG_2.d' in ga5['sources'], ga5
        # every individual measurement is kept, run by run
        pr5 = sorted(ga5['per_run'], key=lambda x: x['file'])
        assert pr5 == [dict(file='Parabens_NEG_1.d', rt_min=1.0),
                       dict(file='Parabens_NEG_2.d', rt_min=2.0)], pr5
        print('case 5 (two parabens runs -> averaged anchors + RT '
              'variation in the report): OK')

        # case 6: EXACT calibrant reference list declared in the .m
        # method (Calibration_LastCalibrationReferenceMass vector) ->
        # DIRECT fit of that list, no adduct-family scan (analyst
        # directive 2026-11-07: the method KNOWS its calibrant). The
        # vector is polarity-independent (cations even in a NEG
        # method): inject the CATION counterparts of the mock's anion
        # combs; d2mz must fit the anion transform (- 2 Na + 2 e)
        root6 = os.path.join(tmp, 'reflist')
        os.makedirs(root6)
        d6, meta6 = T.make_mock_impact_d(root6, n_ms1=16)
        mfile = os.path.join(d6, 'acq.m',
                             'microTOFQImpactAcquisition.method')
        d_ = 2.0 * d2mz.NA_MONO - 2.0 * d2mz.E_MONO   # cation = anion + 2Na - 2e
        vec = ['<para_vec_double permname='
               '"Calibration_LastCalibrationReferenceMass">']
        vec += ['<entry_double value="%.6f"/>' % (lm + d_)
                for lm in meta6['lock_masses']]
        vec.append('</para_vec_double>')
        with open(mfile, 'a') as f:
            f.write('\n'.join(vec) + '\n')
        buf6 = io.StringIO()
        with contextlib.redirect_stdout(buf6):
            rc = d2mz.main([d6, '--out', root6, '--format', 'mzml',
                            '--precision', '64'])
        assert rc == 0, 'd2mz.main failed (case 6)'
        log6 = buf6.getvalue()
        assert 'calibrant reference list read from .m method: 6 masses' \
            ' -> direct fit' in log6, log6
        assert 'negative run: fitting the anion counterparts' in log6, \
            'the cation vector must be transformed for a NEG run:' + log6
        assert 'fitting calibration ladders' not in log6, \
            'family scan ran although the .m declared its exact list'
        assert 'impact quadratic calibration (negative series)' in log6, \
            log6
        mz6 = os.path.join(root6, 'analyses', 'Impact_NEG.mzML')
        assert os.path.exists(mz6)
        T.check_impact_mzml(mz6, meta6)
        print('case 6 (.m reference list -> direct fit, no family '
              'scan): OK')

        # case 7: the SAME mis-referencing the other way around -- a
        # POSITIVE method whose .m vector is the ANION table (analyst
        # 2026-11-07: "cela peut arriver pour le POS"): d2mz must
        # detect the anionic form of the vector and fit the CATION
        # counterparts (+ 2 Na - 2 e)
        root7 = os.path.join(tmp, 'reflist_pos')
        os.makedirs(root7)
        d7, meta7 = T.make_mock_impact_d(root7, n_ms1=16,
                                         name='Impact_POS',
                                         polarity='positive')
        mfile7 = os.path.join(d7, 'acq.m',
                              'microTOFQImpactAcquisition.method')
        neg6 = d2mz._cluster_families(
            d2mz.BRUKER_FA_AC_NA_CLUSTERS['negative'])[0][1][:6]
        vec7 = ['<para_vec_double permname='
                '"Calibration_LastCalibrationReferenceMass">']
        vec7 += ['<entry_double value="%.6f"/>' % lm for lm in neg6]
        vec7.append('</para_vec_double>')
        with open(mfile7, 'a') as f:
            f.write('\n'.join(vec7) + '\n')
        buf7 = io.StringIO()
        with contextlib.redirect_stdout(buf7):
            rc = d2mz.main([d7, '--out', root7, '--format', 'mzml',
                            '--precision', '64'])
        assert rc == 0, 'd2mz.main failed (case 7)'
        log7 = buf7.getvalue()
        assert 'positive run: reference vector is anionic ' \
            '(mis-referenced calibration table) -> fitting the cation ' \
            'counterparts (6 masses >= 50 Da)' in log7, log7
        assert 'fitting calibration ladders' not in log7, \
            'family scan ran although the .m declared its exact list'
        assert 'impact quadratic calibration (positive series)' in log7, \
            log7
        mz7 = os.path.join(root7, 'analyses', 'Impact_POS.mzML')
        assert os.path.exists(mz7)
        T.check_impact_mzml(mz7, meta7)
        print('case 7 (anion table mis-referenced in a POS .m -> '
              'cation transform): OK')

        # case 8: THEORETICAL cluster extension (analyst directive
        # 2026-11-07: integrate the undetected/unlisted clusters to
        # consolidate the calibration at high masses). The enumeration
        # [Na_{b+c+1}(HCOO)_b(AcO)_c]+ must reproduce the 38-mass .m
        # vector EXACTLY and the observed continuation of the
        # calibrant-segment captures (444.9294 -> 906.8258)
        th = d2mz._theoretical_cluster_targets()
        pos = th['positive']
        neg = th['negative']
        vec38 = d2mz.IMPACT_EXACT_TARGETS['positive']
        for v in vec38:
            assert any(abs(m - v) < 0.0005 for m in pos), \
                'vector mass %.6f not covered by the enumeration' % v
        for obs in (498.9012, 512.9167, 566.8887, 634.8762, 702.8636,
                    770.8510, 838.8382, 906.8258, 648.8919, 716.8792,
                    784.8665, 852.8539):
            assert any(abs(m - obs) < 0.001 for m in pos), \
                'observed continuation %.4f missing' % obs
        assert any(abs(m - 906.825727) < 0.001 for m in pos), '906.8257'
        assert any(abs(m - 996.822140) < 0.001 for m in neg), '996.8221'
        assert any(abs(m - 112.98563) < 0.0005 for m in neg), '112.98563'
        assert max(pos) <= 1100.0 and min(pos) >= 50.0
        assert max(neg) <= 1100.0 and min(neg) >= 50.0
        # the fit path logs the extension (case 6 already ran it)
        assert 'extended with' in log6, log6
        print('case 8 (theoretical cluster enumeration, vector '
              'coverage 38/38 + high-mass continuation): OK')
        print('SEGMENT-CAL TEST PASSED')
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == '__main__':
    main()

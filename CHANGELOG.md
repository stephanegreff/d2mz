# Changelog

Versions are dated `YYYY-MM-DD_xx`. The version embedded in `d2mz.py`
(`__version__`) always matches the git tag of the release.

## 2026-11-25_cp
- Durable calibration trace in the report (port of the 2026-11-16_ci
surface onto this lineage): every converted file now carries
calibration_applied, calibration_mode, calibration_n_targets,
calibration_ppm_median, calibration_ppm_max (fit statistics of
the calibrant-segment quadratic), and the batch report gains a
calibration block (ok, n_own_fit, n_borrowed, n_uncalibrated,
worst median/max residual ppm, files_above_5ppm, plain-language
summary) plus the quality_gate block stamped BEFORE the
per-series report dump (report_<series>.json used to carry neither).
- New quality warning: "NO mass calibration applied" for a file whose
trace says applied=False; the end-of-run warning list is factored
into _file_warnings, shared between the per-series dumps and the
final gate.
- This block is the GO/NO-GO gate read before any deposit (GNPS,
Zenodo, HAL): ok=true, n_uncalibrated=0, files_above_5ppm empty.

## 2026-11-25_co
- Fixed a systematic false-positive quality warning: a residual
  `MSMSManual_List*` entry in the acquisition `.method` made every
  never-opened `.d` file (baf-metadata fallback) report
  "targeted MS/MS list (MSMSManual_List*)" even when the run was pure
  auto-DDA. The per-scan `MSMS_Isolation_Act` bookkeeping counter
  (`_modes['iso_act']`) was silently disabling the existing
  data-based retraction on that code path. The retraction now only
  depends on the acquired precursor m/z (real case: 38087,
  auto-DDA with 2839 jittering precursors, wrongly flagged).
- The "!! acquisition outside the validated scope" console banner is
  now printed once, AFTER the retraction has decided: an auto-DDA
  file with a leftover manual list shows the calm "coverage warning
  withdrawn" note instead of a scary double message.
- New helper `manual_list_targeted()` (unit-tested discriminator:
  bit-identical precursor replays); synthetic suite extended with
  never-opened e2e cases (withdrawn in auto-DDA, kept for a real
  targeted list).


## 2026-11-21_cn
- Licence switched from MIT to GPL-3.0: LICENSE file, README,
  GitHub page, CITATION.cff, and the citation templates embedded in
  d2mz (CITATION.cff / .bib / plain-string generated in reports/).

## 2026-11-21_cm
- Compiled-exe behaviour: the console executable pauses at the end of
  the run (and on errors) so the report stays readable on Windows.
- mzML `<run>` startTimeStamp now carries the acquisition start
  (SampleInfo.xml, UTC); mzXML carries the same date as an XML comment.
- Conversion provenance: "Stéphane GREFF, IMBE, Marseille -- BEEMED
  project, ANSES PNR EST 22-EST-063".
- Comment cleanup: English-only, no internal housekeeping tags.

## 2026-11-12_cg/ch
- Two validated Impact-II maxima-block layouts, auto-detected by
  geometry (MS1 = 16n with SNR, MS2/lock = 12n without SNR).
- Acquisition-mode declaration and warnings outside the validated
  scope (DDA CID, constant CE or ramps).

## 2026-11-09 / 2026-11-10
- Never-opened .d support (scan table rebuilt from the baf metadata),
  self-diagnosing metadata errors, per-series reports, quality gate at
  the end of the batch run, Zenodo-ready --zip-output, output folder
  mirroring the archive layout, UV channel selection above 250 nm.

## 2026-11-07 / 2026-11-08
- Impact II calibration from the calibrant segment (Segment == 2)
  with exact .m reference vectors, anion/cation counterpart handling,
  calibration borrowing from a same-mode file, retention-index
  anchoring on the parabens series, every existing trace written to
  the mzML.

## Earlier (2026-10 / 2026-11-06)
- Pure-Python BAF reader, lock-mass quadratic recalibration, precursor
  re-phasing, mzML/mzXML writers, batch mode, reports.

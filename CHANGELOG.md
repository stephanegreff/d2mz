# Changelog

Versions are dated `YYYY-MM-DD_xx`. The version embedded in `d2mz.py`
(`__version__`) always matches the git tag of the release.

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

# d2mz: convert Bruker .d folders to mzML/mzXML — properly recalibrated, precursors included

## The problem this tool solves

If you analyse Bruker qTOF data (micrOTOF-Q, maXis, Impact series) in open mass-spectrometry formats, you have probably met this: **the precursor m/z stored in your MS2 spectra does not match the peaks you actually see in MS1.**

The "precursor m/z" written by the instrument is not a measurement: it is the **isolation setpoint sent to the quadrupole**, computed from the survey scan at acquisition time, before any mass recalibration:

- the exported value is simply the **centre of the isolation window** (whatever width you set in the acquisition method), not the monoisotopic peak that was actually fragmented;
- it is **never recalibrated**, even when you recalibrate your spectra on an internal standard afterwards;
- the accurate value is **not stored anywhere** in the vendor metadata — not in the .baf files, not in `analysis.sqlite`, and Bruker's own Compass export produces the same imprecise value (see [ProteoWizard issue #483](https://github.com/ProteoWizard/pwiz/issues/483));
- ProteoWizard's `precursorRecalculation` filter explicitly refuses time-of-flight analysers — it was designed for Thermo Orbitrap data.

The concrete case we hit in practice is **MZmine**: MS2 scans are attached to MS1 features by matching the precursor m/z within a small tolerance, and an un-recalibrated setpoint can miss the feature entirely, so the MS2 spectra are simply lost for that MS1 feature. The same inherited error affects, to varying degrees, any tool that uses the precursor m/z quantitatively — for example formula-assignment engines (SIRIUS, MS-FINDER) that enumerate candidate formulas within a ppm window around the declared precursor.

## Calibration: what you don't need to do first

Earlier workflows built around the [Garayev mzXML precursor corrector](https://github.com/elnurgar/mzxml-precursor-corrector) required the `.d` files to be **recalibrated in Bruker DataAnalysis before running the script** — either manually, or through the automatic recalibration applied by Bruker DataAnalysis, which can silently apply a poor calibration when the calibrant is weak or absent. Either way it is an operator-dependent step outside the script, repeated for every run. **d2mz removes it entirely.** Feed it the `.d` folder exactly as written by the acquisition: the lock-mass recalibration is performed *inside* the pipeline, **before** any mzML/mzXML is written, so the output files are never the un-recalibrated ones. The calibration actually applied is the one reported in the JSON report (ppm RMS before/after), and it is the same correction for MS1 and MS2 — hence the "phased" output. A recalibration is only trusted when at least 7 lock-mass observations are found — the same number of calibration points Bruker DataAnalysis itself uses — and when a run's calibrant is too weak or absent, d2mz borrows the lock-mass model from the adjacent calibrated file of the same series and says so in the report, instead of silently applying a poor correction.

One caveat users should know: **opening a `.d` in DataAnalysis triggers Bruker's own recalibration** (internal standard, if configured), so the m/z values you inspect interactively in DataAnalysis may differ from those stored in the raw BAF. If your files have been opened and processed through DataAnalysis before conversion, read the d2mz calibration report before comparing: an already-recalibrated run should show a near-zero drift, a raw export typically shows a small one. The report is your audit trail either way.

And a fair warning about instrument state: a `.d` acquired on a **poorly calibrated instrument** (failing reference nebulizer, stale tune) produces strongly drifting — sometimes bizarre-looking — spectra. d2mz models and removes the ppm drift with a **quadratic** (2nd-order) regression over the run, the same order as the qTOF mass calibration; if the drift is wilder than that, the lock-mass residual in the report stays large, which is precisely the flag telling you the instrument needed attention rather than the software. A real-world example `.d` from a mis-calibrated instrument is planned as a validation case in this repository.

## What d2mz does

`d2mz` is a batch converter that reads Bruker `.d` folders **directly** (no proprietary library, no MSConvert, no admin rights, no Python package beyond the standard library) and produces **mzML 1.1.0 and/or mzXML 3.2** files in which:

1. **every spectrum is recalibrated** on your lock masses: the deviation from the expected lock-mass position is regressed against retention time with a quadratic (2nd-order) model — the same regression Bruker itself uses to schematize mass deviations, which are largest at low masses (on a TOF, m/z scales with the square of the flight time). The modelled deviation is then subtracted from every m/z of every spectrum, MS1 and MS2 alike, so the exported values are the corrected ones (the report gives the ppm RMS before and after, so you can prove the correction in a paper);
2. **every precursor m/z declared in the MS2 spectra is re-phased** — the declared setpoint is replaced by the most intense apex of the recalibrated MS1 peaks inside the isolation window (±0.5 Da). A TSV table records declared vs corrected m/z, the delta in ppm and the matched MS1 intensity, so the correction is fully transparent and auditable.

In other words: after conversion, the precursor list of your DDA run is **consistent with the MS1 trace** — what an annotation workflow sees is what the instrument actually fragmented.

```bash
python3 d2mz.py /path/to/runs --format both --precision 32
```

Nothing to type in: for every `.d` file, d2mz automatically detects the polarity and locks on the built-in sodium formate/acetate cluster series of the standard Bruker calibrant (38 positive / 36 negative cluster ions, exact masses validated against the calibrant tables shipped in the `.baf` acquisition files). A folder may mix positive and negative acquisitions; each file is detected and calibrated on its own. `--polarity` only forces a polarity when auto-detection is not wanted; there is deliberately no option to supply a different calibrant.

## Point-and-click

No command line needed: double-click **`d2mz.bat`** on Windows (a small
tkinter window: pick the folder of `.d` acquisitions, press
*Convert*), or grab the **standalone `d2mz.exe`** built
automatically by the repository's GitHub Actions workflow — ~15 MB, no
Python installation required on the acquisition PC. Both drive exactly
the same conversion code as the command-line tool.

## Why you might prefer it to the usual pipeline

- **Zero dependency**: pure Python ≥ 3.8, standard library only. No `libbaf2sql`, no numpy, no pip install, runs on the acquisition PC, a cluster or a lab laptop.
- **Batch**: point it at a folder containing dozens of `.d` runs, come back to a folder of phased, recalibrated mzML/mzXML.
- **32- or 64-bit arrays, zlib-compressed or not** — whatever your downstream tool prefers.
- **Centroid output only**: Bruker's own peak-picked list (declared in the metadata as `MS:1000035` peak picking, `PeakPicking: vendor, MS-levels 1-2`, and `centroided="1"` in mzXML). The raw profile is deliberately not offered — for downstream processing the vendor centroids are the data you want.
- **UV spectra & retention standards ride along inside the mzML**: the DAD/UV data of the `.d` feeds two chromatographic traces (`--uv-trace`, default 254 and 365 nm) written as mzML `<chromatogram>` elements that MS-only software simply ignores — the mzXML output stays free of them. The full wavelength spectrum is kept only at the chromatographic apexes, as human-readable mzML run-level userParams (`d2mz:uv_spectrum`): the mzML is a self-contained archive of the run, and the heavy vendor `.d` folders never need to be shared. A retention-standard series attaches the same way, automatically: when a `parabens` run (a `.d` whose name contains "parabens") sits in the batch folder, d2mz measures the RT of the RI standard series on it (GallAc, pHOHAc, Me/Et/Pr/BuParaben, RI = 100 × carbons; ion [M−H]⁻ in negative, [M+H]⁺ in positive) and stamps every mzML of the folder with the measured anchors — manual `--parabene` entries add or override. The mzML is the only carrier of these data; the JSON report contains processing statistics only, never measurement data.
- **Reproducible**: a JSON report per batch (calibration RMS, number of precursors re-phased) and a CI workflow that re-runs the test suite on every push, so the code cannot silently degrade.

## How it was validated

`test_synthetic.py` builds a synthetic Bruker `.d` from scratch — with a known mass calibration, a deliberately injected **quadratic drift** (+3 ppm linear component + 0.8 ppm curvature) and precursor setpoints offset by **7 ppm** from the true MS1 apexes — then runs the full conversion and asserts that the recalibration residual stays within Bruker's own 2 ppm mass-accuracy guarantee and every precursor lands within 0.2 ppm of its true apex, in both mzML and mzXML, at 32 and 64 bits, compressed and uncompressed. The same tests run automatically on GitHub Actions at every push.

## Background & credits

The BAF decoding is derived from the public reverse engineering of the Bruker BAF format (**H3Jones/bruker-baf**, itself derived from the R package *rtms* by N. Twarog) and from the **baf2sql / kbaf2mzml** knowledge base (bjonnh). The precursor re-phasing strategy follows the lock-mass-based correction philosophy used in plant metabolomics ([mzxml-precursor-corrector](https://github.com/elnurgar/mzxml-precursor-corrector) by E. Garayev), applied here to the case the vendor tooling leaves unsolved. For the methodological framing of precursor mass accuracy in DDA, see Defossez et al. 2023, *Mass Spectrometry Reviews*.

**Authorship**: the d2mz code was written by **Vibe** (AI coding agent powered by GLM, developed by Mistral AI) at the request of, and under the scientific direction of, Stéphane GREFF (IMBE, Aix-Marseille Université). Every conversion output carries this statement in `reports/CITATION.cff`, `CITATION.bib` and `CITATION.txt`.

**Funding**: This work was developed within the BeeMed project, funded by the French National Research Programme for Environmental and Occupational Health (PNR EST) of ANSES (project 22-EST-063).

**Licence**: GPL-3.0 — free to use, modify and redistribute; any redistributed product must in turn be published under GPL-3.0 with its sources (no closed-source embedding in commercial products); please credit the above.

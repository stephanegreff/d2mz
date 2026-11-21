# d2mz — Bruker .d → mzML/mzXML batch converter

[![tests](https://github.com/stephanegreff/d2mz/actions/workflows/tests.yml/badge.svg)](https://github.com/stephanegreff/d2mz/actions/workflows/tests.yml)

Pure-Python (standard library only) batch converter for **Bruker micrOTOF-Q series
`.d` folders** (BAF format), with **lock-mass mass recalibration** and
**precursor re-phasing**: each precursor m/z declared in the MS2
spectra is re-aligned on the
recalibrated MS1 apex it was isolated from.

No dependency, no proprietary Bruker library, no MSConvert export: the BAF
format is decoded directly (index `.baf_idx`, calibration `.baf_xtr`, peak list
`.baf` maxima blocks, `analysis.sqlite` metadata), so MS1 spectra and MS2
precursors are consistent ("phased") in the converted files.

## Why: the precursor phasing problem

In DDA, the precursor m/z written by the instrument is the **quadrupole
isolation setpoint** (the centre of the isolation window, whatever width you set in the method), not a measurement. On Bruker qTOFs the
exported value is the un-recalibrated window centre — the accurate value is not
even stored in the vendor metadata (ProteoWizard issue #483; Defossez et al.
2023, *Mass Spectrometry Reviews*). d2mz therefore recalibrates every scan on
lock masses and re-phases each precursor onto the recalibrated MS1 apex.

## Features

- **Batch**: converts every `*.d` folder under the input path
- **Lock-mass recalibration**: per-scan ppm error on reference peaks →
  quadratic drift model vs RT (2nd-order polynomial — the same order
  Bruker uses to schematize mass deviations, largest at low masses
  where TOF m/z scales with the square of the flight time) → applied
  to MS1 *and* MS2 (report: ppm RMS
  before/after)
- **Precursor re-phasing**: for each MS2, the declared precursor m/z
  (Bruker `Steps.Mass`) is replaced by the apex of the closest recalibrated
  MS1 peak; a TSV table reports scan, declared vs corrected m/z, delta ppm,
  matched MS1 intensity
- **Outputs**: `indexedmzML` 1.1.0 and/or `mzXML` 3.2, 32- or 64-bit
  binary arrays, optional zlib compression, Bruker centroid data only
- **UV spectra & retention standards, all inside the mzML**: the DAD/UV
  data of the `.d` feeds **two chromatographic traces** (`--uv-trace`,
  default 254 and 365 nm) written as mzML `<chromatogram>` elements —
  MS-only software ignores them; the mzXML output stays free of them.
  The **full wavelength spectrum is kept only at the chromatographic
  apexes**, as human-readable mzML run-level userParams
  (`d2mz:uv_spectrum`): the mzML is self-contained, and the heavy `.d`
  folders never need to be archived or shared. **Retention-index
  standards attach the same way, automatically**: when a `parabens`
  run (a `.d` whose name contains "parabens") sits in the folder,
  d2mz measures the RT of the standard series (GallAc, pHOHAc,
  Me/Et/Pr/BuParaben; RI = 100 × carbons) on that run and stamps every
  mzML of the batch with the measured anchors — ion [M−H]⁻ in negative,
  [M+H]⁺ in positive; each file receives the anchors of its own
  polarity, and when only one polarity has a parabens run, its anchors
  are used for the opposite-polarity files as well (the ion label in
  the userParam says which ion they were measured on). Manual
  `--parabene name=rt` entries add or override. All of it lives in the
  mzML only; `report.json` carries processing statistics, never
  measurement data.
- **Reports**: JSON processing report + `<sample>.precursors.tsv`
- **Zero dependency**: Python ≥ 3.8, standard library only
- **Continuous integration**: tests run automatically on every push
  (Python 3.8 / 3.10 / 3.12) via GitHub Actions

## Point-and-click usage

Two click-button options, no command line:

- **`d2mz.bat`** (Windows): double-click it. A small window opens — pick the
  folder containing your `.d` acquisitions, press *Convert*. Requires
  Python 3.8+ installed; polarity detection and lock-mass recalibration
  are fully automatic (built-in Bruker sodium formate/acetate cluster
  series, positive and negative).

### Windows executables (no Python at all)

Every release (tagged `vYYYY-MM-DD_xx`) carries two standalone
executables built by CI:

- **`d2mz.exe`** — console: drag one or several folders (or `.zip`
  archives, or a single `.d`) onto it, exactly like `d2mz_drop.bat`;
  the console pauses at the end of the run (and on errors) so the
  report stays readable.
- **`d2mz_gui.exe`** — graphical interface, double-click, no console.

`d2mz.exe --version` prints the embedded release tag.

The GUI is plain `tkinter` (standard library, file `d2mz_gui.pyw`);
it just drives the same `d2mz.main()` as the command line.

## Usage

```bash
# convert every .d under one or several folders (batch).
# Nothing to set: polarity is detected per file, and the lock-mass
# recalibration uses the built-in sodium formate/acetate cluster series
# of the standard Bruker calibrant.
python3 d2mz.py /path/to/runs /path/to/other_runs --format both --precision 32

# single .d, 64-bit mzML only (precursor re-phasing is automatic)
python3 d2mz.py /path/to/runs/Sample.d --format mzml --precision 64

# which version is this script?
python3 d2mz.py --version
```

Each input folder receives its own `d2mz_output_<source-name>/` folder,
mirroring the Zenodo-ready archive layout:

```
<input folder>/
└── d2mz_output_<source-name>/
    ├── analyses/                  # every mzML / mzXML, FLAT
    └── annexes/
        ├── TSV/                   # one precursor table per analysis
        ├── report_<source-name>.json   # per-series processing report
        └── CITATION.cff / .bib / .txt
```

| Option | Default | Description |
|---|---|---|
| `--out` | `d2mz_output` | output folder created inside each input folder (or absolute path) |
| `--format` | `both` | `mzml` / `mzxml` / `both` |
| `--precision` | `32` | binary array float precision (32/64) |
| `--no-compress` | off | disable zlib compression |
| `--uv-trace` | `254`, `365` | wavelengths (nm) of the two UV traces (repeatable) |
| `--parabene` | — | retention-time standard, repeatable, e.g. `--parabene C4=12.34` (minutes); stored as mzML userParams + report. Automatic when a `parabens` run is in the folder (measured anchors) |
| `--no-prec-correct` | off | disable precursor re-phasing |
| `--prec-window` | `0.5` | precursor search window half-width (Da) — the isolation window, a physical bound |
| `--polarity` | `auto` | `auto` (detected per file) / `positive` / `negative` |
| `--report` | — | JSON processing report path |
| `--probe FOLDER.d` | — | diagnostic: dump the calibration records found in this .d folder and exit (no conversion) |

## Outputs per `.d` (inside `<input>/d2mz_output/`)

- `converted/<sample>.mzML` — indexedmzML 1.1.0 (little-endian arrays,
  valid for MZmine, GNPS, MS-DIAL, pyteomics, pymzML…)
- `converted/<sample>.mzXML` — mzXML 3.2 (network byte order, as per spec)
- `reports/<sample>.precursors.tsv` — precursor re-phasing table
  (declared vs corrected m/z, delta ppm, MS1 intensity)
- `reports/report.json` — per-file stats incl. lock-mass calibration RMS
  (before/after); an extra copy can be written with `--report`
- `reports/CITATION.cff` / `.bib` / `.txt` — how to cite d2mz

## Post-run hook: convert automatically on the acquisition PC

`d2mz_postrun.vb` is a small **Bruker DataAnalysis script (Visual Basic
.NET)** to save in your acquisition method: at the end of each run it
launches the standalone `d2mz.exe` on the just-written `.d` folder — no
export, no DataAnalysis session, no Python installation.

1. put `d2mz.exe` at a fixed location (e.g. `C:\d2mz\d2mz.exe`) and set
   the `D2MZ_EXE` constant at the top of the script accordingly;
2. in DataAnalysis, open the script editor, paste the script, save it in
   the acquisition method (menus differ slightly between DA 4.x and 5.x);
3. check the two lines marked `--- CHECK` (the analysis-path property
   name) against the *DataAnalysis Script Reference Manual* of your
   version — the candidates are already in the script.

The converted files land in `<run folder>\d2mz_output\` as usual.

## Data embedded in the mzML (self-describing format)

Besides the MS spectra, each mzML written by d2mz carries the non-MS
information of the run, in plain text that any tool ignores and any
reader can extract. No dedicated software is needed — the format is
self-describing:

- **UV traces** — `<chromatogramList>`: two `<chromatogram>` elements
  (ids `chromatogram=0/1`), each with a `<userParam name="wavelength">`
  label (e.g. `254 nm`), a time array (seconds, `MS:1000595`) and an
  intensity array (`MS:1000515`), same precision and zlib compression as
  the spectra. Intensities are **raw** (as acquired); retention times are
  those of the acquisition's scan table (no offset applied).
- **UV spectra at chromatographic apexes** — one mzML run-level
  `<userParam name="d2mz:uv_spectrum">` per apex, value readable in any
  text editor:
  `trace=254 nm; rt=3.5000 min; data=wl:intensity,wl:intensity,…`
  wavelengths in nm, intensities raw, sorted by wavelength. The full
  spectrum is kept only at the apexes (no per-scan bloat).
- **Retention-time standards** — either measured automatically from a
  `parabens` run of the folder, or given with `--parabene`: one
  `<userParam name="d2mz:retention_standard:NAME" value="…"/>` per
  standard. The value of a measured anchor reads
  `rt 3.0000 min; ri 1100; mz 193.0873 ([M-H]-); source Parabens_NEG_29042.d`
  (RI = 100 × carbons; manual entries carry `rt X.XXXX min` only).

The mzML is the only carrier of these data; `report.json` contains
processing statistics only (calibration RMS, scan counts, trace point
counts), never measurement data.

Nothing else is embedded; the mzXML output carries none of it (its
schema has no room for it).

## Method (for a manuscript)

1. **Decode**: BAF index/calibration/peak blocks + `analysis.sqlite`
   (RT, precursor `Steps.Mass`, collision energy `PerSpectrumVariables`).
2. **Recalibrate**: for each MS1 scan, the most intense peak within
   tolerance of each lock mass gives a ppm error; the per-scan median is
   deviation from the expected mass regressed against RT with a
   2nd-order polynomial (Bruker schematizes mass deviations this way,
   most visible at low masses); the model removes the drift from all
   spectra of the run (same correction applied to MS2 product ions).
3. **Re-phase precursors**: each precursor m/z declared in the MS2
   spectra is matched
   against the peaks of the last preceding recalibrated MS1 scan; the
   apex inside the isolation window (+-0.5 Da) replaces the declared
   value — the most intense MS1 apex wins, automatically, and the
correction
   is reported. This makes DDA precursor lists consistent with the
   (recalibrated) MS1 trace — a prerequisite for unambiguous annotation
   of positional isomers.

## Tests

```bash
python3 test_synthetic.py
```

Builds a synthetic `.d` (known calibration, injected quadratic drift —
+3 ppm linear component + 0.8 ppm curvature,
declared precursors offset by 7 ppm), runs the full conversion and asserts:
recalibration residual < 2 ppm (Bruker's own mass-accuracy guarantee),
precursor error < 0.2 ppm, valid
mzML/mzXML round-trip. The same tests run automatically via GitHub Actions
(`.github/workflows/tests.yml`) on Python 3.8 / 3.10 / 3.12 at each push
and pull request.

## Repo layout (GitHub)

```
d2mz/
├── .github/workflows/tests.yml
├── .gitignore            (__pycache__/, *.pyc, converted/, report.json)
├── LICENSE               (GPL-3.0)
├── README.md
├── d2mz.py
└── test_synthetic.py
```

## Credits & licence

BAF reverse engineering derived from the public work of **H3Jones
(bruker-baf)**, itself derived from the R package **rtms** (N. Twarog), and
from the **baf2sql / kbaf2mzml** knowledge (bjonnh). Precursor re-phasing
inspired by lock-mass-based precursor correction strategies used in plant
metabolomics ([mzxml-precursor-corrector](https://github.com/elnurgar/mzxml-precursor-corrector) by E. Garayev).

**Authorship**: the d2mz code was written by **Vibe** (AI coding agent
powered by GLM, developed by Mistral AI) at the request of, and under the
scientific direction of, Stéphane GREFF (IMBE, Aix-Marseille Université).
The same statement is embedded in the `CITATION.cff` / `.bib` / `.txt`
files generated in every `reports/` folder.

**Transparency**: d2mz was developed through an interactive collaboration
between Stéphane GREFF (IMBE, Marseille) and Vibe, across successive versions
since November 2026, each version being validated by automated tests before
release. The explanatory notice (`docs/d2mz_explication.pdf`) was likewise
written by Vibe from the source code of each version.

**Funding**: This work was developed within the BeeMed project, funded by the French National Research Programme for Environmental and Occupational Health (PNR EST) of ANSES (project 22-EST-063).

GPL-3.0 licence — see `LICENSE`. Everyone may use, modify and redistribute the code freely; any redistributed product must in turn be published under GPL-3.0 with its sources, so no closed-source embedding in commercial products.

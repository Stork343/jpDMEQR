# Manuscript files

This directory contains the current SJS/Wiley working manuscript and Supporting Information.

## Principal files

- `sjs_profile_quantile_main.tex`
- `sjs_profile_quantile_supplement.tex`
- `sections/` (modular main-text and supplement sections included by the two drivers)
- `references.bib`
- compiled PDFs are generated locally from the two source files and are not required by the computational pipeline
- `REVISION_NOTES_zh.md`

The pre-reconstruction manuscript is preserved on the Git branch `legacy-pre-sjs-reconstruction`; large legacy PDFs are not duplicated in `main`.

## Compilation

The source uses the Wiley NJDv5 class with `APA,STIX1COL`. Font files are not stored in this repository. The small compatibility file `lettersp.sty` is included. Retrieve the non-font template support files with:

```bash
bash fetch_wiley_template_files.sh
```

For submission, compare those files with the current Wiley NJDv5 package distributed by Wiley, copy Wiley's `Fonts/` directory locally into this directory, and run:

```bash
bash compile.sh
```

## Status of numerical sections

Sections 5 and 6 are protocols with deliberate placeholders. Results must be inserted only from validated files produced under `results/` and must pass the output contract. The manuscript must never be edited to contain hand-entered simulation or application numbers.

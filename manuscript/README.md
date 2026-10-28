# Joint profile manuscript

`sjs_profile_quantile_main.tex` builds the article and `sjs_profile_quantile_supplement.tex` builds Supporting Information. Component sources are in `sections/`; `references.bib` supplies both bibliographies. The author and affiliation blocks remain in the two document entry files.

The current estimator is `joint_profile_design_20261026`. Its complete definition is in `docs/METHOD_SPECIFICATION.md`. The application uses the fixed-design coefficient interface and the separately tuned prediction learner in `research/joint_profile_manuscript_20261026/`.

## Rebuild numerical content

From the repository root, with the saved numerical records available, run:

```bash
TMPDIR=results/scratch Rscript research/joint_profile_manuscript_20261026/population_roots.R
TMPDIR=results/scratch Rscript research/joint_profile_manuscript_20261026/build_inference_assets.R
TMPDIR=results/scratch Rscript research/joint_profile_manuscript_20261026/build_prediction_assets.R
```

The builders require complete keys for the independent evaluation, its paired ablation, all module workflows, exploratory roles and final prediction fits. Numerical tables and text values are written to `generated/joint_profile/`. Figure PDFs are written to `generated/figures/`. The corresponding unrounded source tables are in `results/tables/joint_profile_manuscript_20261026/`.

The cohort and module-membership tables retain the completed metadata audit and public membership records. The conditional-scale population table uses the validated source in `results/tables/population_targets/manuscript_20261023_v1/`; its population parameter is the same unsmoothed complete profile target used in the article.

## Compile and inspect

From this directory, run `bash compile.sh` with XeLaTeX and BibTeX on `PATH`. The script compiles both documents and resolves their bibliography and cross references. Keep temporary files within the repository's ignored scratch directory by setting `TMPDIR=../results/scratch`.

From the repository root, run `research/joint_profile_manuscript_20261026/check_manuscript.py` with Python containing `pypdf`, `pdfplumber`, `pylatexenc` and `bibtexparser` version 2. The checker uses the installed `latexpand` utility to read the complete multifile source. It checks citation keys, labels, PDF text, embedded fonts and character bounds. Its report explicitly identifies programmatic inspection; no visual inspection is performed.

## Scope of the numerical records

The article reports the completed 24-setting independent evaluation with 60 datasets per setting and a paired removal of conditional covariate adjustment on the same datasets. It also reports the completed new-method application calculations. The separate 78-setting full method comparison remains pending. Its registry and population-target assets are retained in the repository.

The published-method application records contain one outer fold in one repetition. Their actual subject counts and convergence states appear in Supporting Information. The complete new-method prediction workload is reported separately. A returned comparator fit that reports nonconvergence has an explicit unavailable loss comparison.

## Reproduction sources

`research/joint_profile_manuscript_20261026/README.md` gives the execution entry points, output directories and numerical checks. Generated PDF files can be read directly; the TeX sources remain editable in this directory.

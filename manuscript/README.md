# Joint profile manuscript

`sjs_profile_quantile_main.tex` builds the article and `sjs_profile_quantile_supplement.tex` builds Supporting Information. Component sources are in `sections/`; `references.bib` supplies both bibliographies. The author and affiliation blocks remain in the two document entry files.

The current estimator is `joint_profile_design_20261026`. Its complete definition is in `docs/METHOD_SPECIFICATION.md`. The application uses the fixed-design coefficient interface and the separately tuned prediction learner in `research/joint_profile_manuscript_20261026/`.

## Rebuild numerical content

The completed registered new-method queue contains 45,500 datasets in
71 statistical settings. Its tables can be regenerated independently of
the comparison-method queues:

```bash
TMPDIR=results/scratch/runtime Rscript research/submission_completion_20261027/collect_simulation.R results/raw/simulation/submission_20261027_main_joint
TMPDIR=results/scratch/runtime Rscript research/submission_completion_20261027/summarise_simulation_collection.R results/tables/submission_20261028/collections/submission_20261027_main_joint/records.rds
TMPDIR=results/scratch/runtime Rscript research/submission_completion_20261027/audit_prescreen_targets.R
TMPDIR=results/scratch/runtime Rscript research/submission_completion_20261027/build_joint_study_assets.R
TMPDIR=results/scratch/runtime Rscript research/submission_completion_20261027/build_target_table.R
```

The collection step verifies every registered dataset and method identity.
It retains numerical failures in the coverage denominator. The generated
tables report the complete setting-level coverage and interval lengths,
with individual-coordinate binomial intervals and standard error ratios.
The target-displacement calculation uses the stored F05 and F06 supports
and their generating Gaussian covariance.

The 40 completed reference-processing module workflows have a separate
builder:

```bash
TMPDIR=results/scratch/runtime Rscript research/submission_completion_20261027/build_reference_module_assets.R
```

It verifies all 5,200 full-sample and subject-deletion fits, recomputes the
workflow standard errors from the deletion estimates, and checks the
simultaneous interval allocation before writing the reference module table.
The table builders use `knitr` and `kableExtra` to generate LaTeX tables
with repeated column headings across pages. Both packages are declared in
`DESCRIPTION` and `scripts/00_install_dependencies.R`.

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

The Python requirements are listed in
`research/joint_profile_manuscript_20261026/requirements.txt`. For an
installation within this project, run:

```bash
mkdir -p results/scratch/runtime
TMPDIR="$PWD/results/scratch/runtime" python3 -m pip install --pre --target results/scratch/python_dependencies -r research/joint_profile_manuscript_20261026/requirements.txt
PYTHONPATH="$PWD/results/scratch/python_dependencies" python3 research/joint_profile_manuscript_20261026/check_manuscript.py
```

## Scope of the numerical records

The registered statistical study retains all 71 settings and their
45,500 new-method datasets. Seven separate computation settings complete
the 78-setting registry. The simulation comparison queues have their own
completion records; their final combined summaries require every registered
method record.

The independent evaluation uses 24 settings with 60 datasets per setting
and a paired removal of conditional covariate adjustment. Its builder
requires all 1,440 matched dataset keys. Application tables combine 40
reference-processing module workflows with twelve alternative-processing
workflows, and retain every subject-deletion fit. Prediction and
exploratory results have separate complete-key checks.

The generated manuscript fragments retain their source manifests. A full
reproduction requires the corresponding raw records for each study.
The completed reference-processing module results and registered
new-method simulations can be analysed while the remaining calculation
queues continue.

## Reproduction sources

`research/joint_profile_manuscript_20261026/README.md` gives the execution entry points, output directories and numerical checks. Generated PDF files can be read directly; the TeX sources remain editable in this directory.

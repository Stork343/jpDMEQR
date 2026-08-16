# Empirical protocol: GSE65391 longitudinal paediatric SLE study

## 1. Scientific objective

The application asks whether whole-blood transcriptional features are associated with different parts of the contemporaneous complement C3 distribution after accounting for repeated observations, demographics, follow-up time, treatment and low-dimensional clinical covariates. The primary focus is the lower and central portions of the C3 distribution, where complement depletion may be most pronounced.

The analysis is observational. Estimated coefficients are regularised profile associations. They are not causal treatment effects and are not interpreted as biological mechanisms without external validation.

## 2. Why GSE65391 fits the design

The GEO series is a longitudinal paediatric systemic lupus erythematosus study with clinical parameters. The source record reports whole-blood transcriptomes for 158 SLE patients followed for up to four years, 924 SLE expression samples, and healthy samples/technical replicates used in the study. The platform is Illumina HumanHT-12 v4.0. Sample metadata include subject, visit, cumulative time, age, sex, race, treatment, SLEDAI, C3, C4 and numerous laboratory variables.

This structure is materially better aligned with the paper than the legacy dataset because it has many independent subjects, repeated visits, a bounded-cluster longitudinal structure and a continuous clinical outcome available at the same visit as high-dimensional expression measurements.

## 3. Official data sources

Accession page:

```text
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65391
```

Processed Series Matrix:

```text
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE65nnn/GSE65391/matrix/GSE65391_series_matrix.txt.gz
```

Supplementary directory:

```text
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE65nnn/GSE65391/suppl/
```

Raw/custom archive:

```text
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE65nnn/GSE65391/suppl/GSE65391_RAW.tar
```

Non-normalized expression files:

```text
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE65nnn/GSE65391/suppl/GSE65391_non-normalized_data_Illumina_HT12_V4_R1.txt.gz
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE65nnn/GSE65391/suppl/GSE65391_non-normalized_data_Illumina_HT12_V4_R2.txt.gz
```

The default pipeline downloads the processed Series Matrix first because it contains the batch-normalised values described in the GEO sample records and is sufficient for cohort auditing and a primary reproducible analysis. The raw/non-normalized files are optional sensitivity inputs and are substantially larger.

## 4. Download and provenance

Run:

```bash
Rscript scripts/application/00_download_GSE65391.R --processed=true --raw=false
```

Optional raw sensitivity download:

```bash
Rscript scripts/application/00_download_GSE65391.R --processed=true --raw=true
```

The script writes `data/raw/GSE65391/download_manifest.csv` with:

- accession;
- source URL;
- local path;
- download timestamp in UTC;
- file size;
- MD5 checksum;
- success/failure status;
- R/GEOquery version.

A partial or checksum-changing download is not used without an audit note.

## 5. Data construction

Run:

```bash
Rscript scripts/application/01_build_GSE65391.R
```

The build script performs the following steps.

### 5.1 ExpressionSet extraction

1. Load the Series Matrix through `GEOquery::getGEO` or the downloaded matrix file.
2. Confirm one platform and select `GPL10558`.
3. Extract the probe-by-sample expression matrix.
4. Extract GEO phenotype columns and parse every `characteristics_ch1*` field of the form `key: value` into a data dictionary.
5. Preserve original strings and create cleaned typed variables separately.

### 5.2 Identifier reconciliation

Construct unique visit identifiers from GEO accession, subject and visit. Verify:

- every expression column has one phenotype row;
- no two SLE biological visits share the same final identifier;
- technical replicates are explicitly labelled;
- subject and visit numbers agree with sample titles;
- rows are sorted only after identifiers have been checked.

Discrepancies stop the pipeline and are written to a reconciliation report.

### 5.3 Cohort groups

Create mutually exclusive flags:

- SLE biological sample;
- healthy biological sample;
- healthy technical replicate;
- SLE batch replicate, if present;
- unresolved.

The outcome model uses SLE biological visits only. Healthy samples and technical replicates are retained for QC and batch diagnostics.

## 6. Primary outcome gate

The provisional primary response is

```text
Y_ij = log(C3_ij)
```

for the contemporaneous visit. Log transformation is used only if all retained C3 values are positive. Before any expression-outcome association is fitted, run:

```bash
Rscript scripts/application/02_audit_GSE65391.R
```

The C3 analysis proceeds only if all frozen gates pass:

1. at least 100 SLE subjects have usable expression and C3;
2. at least 500 eligible visits remain;
3. at least 75 subjects have two or more eligible visits;
4. C3 missingness among SLE expression visits is at most 30%;
5. all retained C3 values are positive;
6. the outcome has at least 20 distinct values and no single value accounts for more than 20% of visits;
7. each primary quantile has enough local observations for a stable curvature estimate;
8. no source batch contains only one outcome region after cohort restrictions.

If any gate fails, the script stops with `OUTCOME_GATE_FAILED`. It does not automatically replace C3 with C4, SLEDAI or another outcome. A change of primary outcome requires a committed protocol amendment before fitting.

### 6.1 Target quantiles

Primary:

```text
tau = 0.25, 0.50, 0.75
```

Sensitivity:

```text
tau = 0.10
```

The 0.10 analysis is run only when effective sample size and local-density diagnostics meet the frozen threshold.

## 7. Clinical variables

### 7.1 Always-included fixed covariates

Main adjustment set:

- centred cumulative time since baseline;
- age at visit or baseline age plus time, with no duplicate parameterisation;
- sex;
- race/ancestry categories with small levels pooled by a prespecified rule;
- treatment indicators or a reproducible treatment-category encoding;
- source batch;
- training/test-set indicator if residual source-processing differences remain.

All have penalty factor zero.

### 7.2 SLEDAI adjustment

Because SLEDAI may lie on the pathway between transcriptional activity and complement depletion or may act as a broad severity confounder, it is not automatically included in the primary model. Two prespecified models are reported:

- Model A: demographic, time, treatment and batch adjustment;
- Model B: Model A plus contemporaneous SLEDAI.

Differences are interpreted as sensitivity to clinical-severity adjustment, not as evidence of causality.

### 7.3 Random-effect design

Primary:

```text
random intercept by subject
```

Sensitivity:

```text
random intercept + centred cumulative-time slope
```

The slope model is attempted only if the visit distribution supports it. The following are recorded: subjects with at least three visits, within-subject time variation, nuisance-block condition numbers and convergence rate. If the slope model is unstable, the intercept model remains primary and the failure is reported.

## 8. Expression preprocessing

### 8.1 Primary processed-data route

1. Use the processed/batch-normalised expression values in the Series Matrix.
2. Confirm their scale, finite values and sample alignment.
3. Use healthy technical replicates to quantify residual batch and replicate disagreement.
4. Obtain probe annotation from `GPL10558`, `illuminaHumanv4.db`, or a frozen platform-annotation file with version recorded.
5. Remove control probes, probes without a gene mapping, probes mapping to multiple genes, and probes failing detection/variation rules.
6. Collapse multiple probes mapping uniquely to the same gene using a rule learned on training subjects only. Primary rule: retain the probe with largest training-subject median absolute deviation. Sensitivity: median across mapped probes.
7. Apply an outcome-blind training-fold variance filter. The exploratory high-dimensional block retains the top 5000 genes by training-fold MAD, with 2000 and 10000 as sensitivity values.
8. Standardise genes using training-subject means and standard deviations; apply the same transformation to validation subjects.

### 8.2 Raw/non-normalized sensitivity route

The optional route uses the R1/R2 non-normalized files. It must document:

- separation of expression and detection-p-value columns;
- background and control-probe handling;
- `limma::neqc` or another prespecified Illumina-appropriate normalisation;
- cross-batch adjustment described in the GEO record;
- technical-replicate agreement before and after processing;
- concordance with the processed Series Matrix.

The raw route is a sensitivity analysis. It does not replace the processed route after observing which one gives more favourable associations.

### 8.3 Leakage prevention

No operation using outcome values or all-subject distributional information may be performed before subject splitting. In particular:

- no full-data gene ranking by C3 association;
- no full-data standardisation before cross-validation;
- no ComBat including validation subjects when evaluating prediction;
- no probe-to-gene selection based on full-data outcome association.

Outcome-blind platform annotation may use all samples. Any data-adaptive distributional filter is learned on training subjects.

## 9. Confirmatory and exploratory analyses

### 9.1 Confirmatory low-dimensional layer

The confirmatory layer uses a frozen list of immune-module scores chosen before examining C3 associations. At minimum, the list should include interferon-response modules with public, reproducible gene sets. Additional plasmablast or neutrophil modules require exact, versioned gene definitions.

Module scores are calculated within training folds using a frozen rule, such as the mean of standardized member genes when a minimum proportion of genes is observed. These scores are unpenalised and receive profile confidence intervals at each quantile.

The module definition file must contain:

```text
module_id, source, source_version, gene_symbol, inclusion_reason
```

An empty or changing module file blocks confirmatory claims.

### 9.2 Exploratory gene-level layer

The exploratory layer uses high-dimensional genes with the l1 penalty. It has two tasks:

1. subject-held-out prediction and stability assessment;
2. split-sample coordinate inference for a small number of genes selected on independent subjects.

The same subjects are never used both to select a gene and to form its ordinary pointwise Wald interval within a split. Roles are swapped across repeated splits and combined by a frozen aggregation rule. The result is labelled exploratory split-sample inference.

Ordinary BH-FDR applied to thousands of same-sample debiased p-values is not part of the primary analysis.

## 10. Subject splitting and validation

### 10.1 Frozen folds

Generate repeated five-fold subject-level cross-validation with five repetitions. Stratification may use source batch and coarse baseline C3 category, but never individual visits. Save the fold table before model fitting.

### 10.2 Tuning

Within each outer training set:

- select `lambda_beta` by inner subject-level CV on pinball loss;
- use frozen `h` and `Lambda` grids;
- choose Dantzig `mu` by the smallest feasible prespecified multiplier;
- choose screening dimension by a frozen grid evaluated within training subjects;
- do not examine outer-fold outcomes when tuning.

### 10.3 Predictive metrics

At each quantile report:

- held-out pinball loss;
- paired difference from pooled QR-Lasso;
- skill relative to intercept-only quantile prediction;
- empirical quantile calibration;
- subject-level distribution of loss differences;
- selected model size and selection stability;
- convergence/failure rate;
- runtime.

New-subject prediction does not use latent true random effects. A conditional follow-up prediction analysis may estimate random effects from earlier visits, but it must be separately labelled.

## 11. Comparators

Primary empirical comparators:

1. pooled high-dimensional quantile Lasso;
2. longitudinal QGEE-SCAD or the closest faithful implementation of Zu et al.;
3. low-dimensional LQMM on prespecified clinical/module covariates;
4. proposed profile method;
5. proposed profile method with random-intercept-only versus intercept-slope sensitivity.

Bayesian Gaussian mixed shrinkage is supplemental and reported as a mean-model prediction benchmark, not as quantile inference.

## 12. Required audit outputs

Before fitting:

- source sample counts by group and technical-replicate status;
- parsed metadata dictionary;
- C3/C4/SLEDAI missingness;
- eligible-subject and eligible-visit counts;
- cluster-size and follow-up-time distributions;
- treatment-category frequencies;
- race/sex distributions;
- batch-by-group and batch-by-outcome summaries;
- technical-replicate agreement;
- expression scale and probe-annotation flow;
- exact number of probes/genes removed at each step;
- frozen subject-fold table.

A cohort-flow table must reconcile the 996 source samples to the final analysis visits.

## 13. Required result outputs

Primary manuscript outputs:

1. cohort and preprocessing flow table;
2. profile estimates and 95% intervals for confirmatory modules across quantiles;
3. subject-held-out pinball-loss comparison table;
4. coefficient-by-quantile figure for confirmatory modules;
5. paired held-out loss-difference figure;
6. exploratory split-sample gene table with stability frequencies;
7. optimisation, Dantzig and profile-identity diagnostics;
8. sensitivity table for SLEDAI adjustment, random slope and preprocessing route.

All numbers are generated from scripts. No manual result transcription is allowed.

## 14. Stop conditions

The empirical analysis stops for review if:

- the C3 outcome gate fails;
- sample/metadata identifiers cannot be reconciled;
- the processed expression scale is inconsistent across source batches and cannot be diagnosed;
- technical-replicate disagreement indicates unresolved normalization failure;
- fewer than 90% of primary proposed-method fits converge;
- the profile identity test fails on empirical design subsets;
- Dantzig feasibility is below 90% for confirmatory coordinates;
- subject-level folds are substantially imbalanced and cannot be repaired without outcome leakage.

A stop condition produces a report; it does not trigger an undisclosed change of outcome, model or preprocessing.

# Empirical analysis freeze decisions

**Application:** GSE65391 longitudinal paediatric SLE transcriptomic study  
**Primary provisional outcome:** contemporaneous complement C3  
**Status:** data download, construction, and audit are authorised; confirmatory/full fitting is **not** authorised until every empirical gate in this document passes  
**Companion files:** `docs/EMPIRICAL_PROTOCOL_GSE65391.md`, `config/application.yml`, and `docs/METHOD_SPECIFICATION.md`

## 1. Purpose and authority

This document freezes the decisions that must be made before inspecting expression--outcome associations. It prevents outcome switching, preprocessing leakage, post hoc module selection, and empirical fitting with an implementation that has not passed the simulation geometry gates.

Authority order for the application is:

1. `docs/METHOD_SPECIFICATION.md` for the estimator and inferential target;
2. this document for frozen empirical decisions and stop rules;
3. `docs/EMPIRICAL_PROTOCOL_GSE65391.md` for the full processing protocol;
4. `config/application.yml` for machine-readable settings;
5. `docs/RESULTS_CONTRACT.md` for outputs;
6. application code;
7. manuscript prose;
8. legacy GSE121239 code/results.

The legacy GSE121239 analysis is not reused. In particular, AFFX/control-probe findings, old FDR tables, and old calibration results are not carried into the new application.

## 2. Frozen scientific question

The application asks whether whole-blood transcriptional features are associated with different parts of the **contemporaneous C3 distribution** after accounting for repeated measurements and prespecified demographic, time, treatment, batch, and clinical-severity covariates.

The estimand is a regularised profile association. It is not a causal treatment effect, a mechanistic effect, or an unregularised random-effect parameter.

## 3. Frozen data source

Use GEO accession `GSE65391`, platform `GPL10558`.

Primary source:

```text
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE65nnn/GSE65391/matrix/GSE65391_series_matrix.txt.gz
```

Additional provenance/sensitivity sources:

```text
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65391
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE65nnn/GSE65391/suppl/
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE65nnn/GSE65391/suppl/GSE65391_RAW.tar
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE65nnn/GSE65391/suppl/GSE65391_non-normalized_data_Illumina_HT12_V4_R1.txt.gz
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE65nnn/GSE65391/suppl/GSE65391_non-normalized_data_Illumina_HT12_V4_R2.txt.gz
```

The processed Series Matrix is the primary route. Raw/non-normalized processing is a prespecified sensitivity route; it may not replace the primary route after observing which produces more favourable associations.

Every downloaded file must have source URL, UTC timestamp, byte size, checksum, local path, and status in `data/raw/GSE65391/download_manifest.csv`.

## 4. Cohort freeze

### 4.1 Included observations

The outcome analysis includes SLE biological visits with:

- reconciled subject and visit identifiers;
- one expression profile assigned to the biological visit after a frozen replicate rule;
- a valid contemporaneous C3 value;
- required adjustment variables or an explicitly documented missing-data treatment;
- no unresolved group/batch/replicate status.

### 4.2 Excluded from the outcome model

- healthy biological samples;
- healthy technical replicates;
- technical/batch replicates not chosen by the frozen replicate rule;
- unresolved samples;
- visits without a valid contemporaneous C3 value;
- duplicate biological visits after identifier reconciliation.

Healthy and technical-replicate samples may be used for QC and normalization diagnostics only.

### 4.3 Equal subject-level split unit

All training, validation, screening, and inference roles are assigned at the subject level. No visit from the same subject may appear in two roles within a split.

## 5. Primary outcome and non-switching rule

The provisional primary response is

\[
Y_{ij}=\log(C3_{ij}),
\]

provided all retained C3 values are positive. Run the audit before any expression--C3 model.

The C3 analysis proceeds only if all gates pass:

1. at least 100 SLE subjects with usable expression and C3;
2. at least 500 eligible SLE visits;
3. at least 75 subjects with at least two eligible visits;
4. C3 missingness among SLE expression visits at most 30%;
5. all retained C3 values positive;
6. at least 20 distinct outcome values;
7. no single value accounts for more than 20% of eligible visits;
8. sufficient observations near each primary quantile for stable smoothed curvature;
9. no source batch is confined to one outcome region after restrictions;
10. subject/visit reconciliation has no unresolved critical conflict.

A failed gate produces `OUTCOME_GATE_FAILED` and stops confirmatory fitting. The pipeline must not automatically switch to C4, SLEDAI, flare status, or another outcome. Any new primary outcome requires a committed protocol amendment made before association fitting on that outcome.

## 6. Frozen quantiles

Primary:

```text
tau = 0.25, 0.50, 0.75
```

Sensitivity:

```text
tau = 0.10
```

The `tau=0.10` fit is run only when the audit confirms adequate local effective sample size and curvature. It is never substituted for a failed primary quantile.

## 7. Frozen clinical adjustment models

All listed low-dimensional covariates have penalty factor zero.

### Model A

- centred cumulative time since baseline;
- age represented once, either age at visit or baseline age plus time;
- sex;
- race/ancestry under a frozen small-cell pooling rule;
- reproducibly encoded treatment indicators/categories;
- source batch or processing batch.

### Model B

Model A plus contemporaneous SLEDAI.

Model B is a clinical-severity sensitivity analysis. Differences between Models A and B are not interpreted causally.

If a planned variable is unavailable or too sparse, its handling must be determined in the metadata audit report and committed before expression--C3 fitting. Do not decide covariates by p-values or predictive improvement on the full outcome data.

## 8. Frozen random-effect design

Primary:

```text
random intercept by subject
```

Sensitivity:

```text
random intercept + centered cumulative-time slope
```

The slope sensitivity is attempted only if:

- enough subjects have at least three eligible visits;
- within-subject time variation is nondegenerate;
- the nuisance block is numerically identified;
- the strict profile identity continues to pass on empirical-design subsets;
- the convergence rate is acceptable under the frozen criterion.

An unstable slope model does not replace or invalidate the primary random-intercept model. Its instability is reported as a sensitivity result.

## 9. Frozen expression-processing routes

### 9.1 Primary processed route

1. Load the processed Series Matrix.
2. reconcile expression columns with phenotype rows;
3. confirm scale, finite values, and platform;
4. remove control probes;
5. remove probes without a unique gene mapping;
6. remove probes mapped to multiple genes;
7. apply any detection/variation rule specified before model fitting;
8. collapse multiple uniquely mapped probes to one gene using the probe with the largest **training-subject** MAD;
9. apply an outcome-blind training-fold MAD/variance filter;
10. standardise using training-subject means and standard deviations;
11. apply the training transformation to validation/inference subjects.

Primary retained high-dimensional size:

```text
5000 genes
```

Sensitivity sizes:

```text
2000 and 10000 genes
```

### 9.2 Raw/non-normalized sensitivity route

The raw route must predeclare and record:

- expression and detection-p-value column separation;
- control/background handling;
- the exact normalization algorithm and package versions;
- batch adjustment;
- technical-replicate concordance before/after processing;
- concordance with the processed route.

The raw route is labelled sensitivity and does not become primary after results are observed.

## 10. Leakage rules

The following operations are fold-local and use outer-training subjects only:

- probe collapse when data-adaptive;
- MAD/variance filtering;
- standardisation;
- outcome-guided screening;
- selection of `lambda_beta`, screening dimension, or other predictive tuning;
- any batch correction that estimates distributional parameters from subjects;
- module score transformations that depend on empirical distributions.

Permitted global operations are limited to outcome-blind platform annotation and deterministic identifier reconciliation.

Prohibited:

- ranking genes by C3 association on all subjects;
- full-data standardisation before cross-validation;
- ComBat/normalization using validation subjects when estimating held-out performance;
- selecting a probe among mapped probes using full-data outcome association;
- using the same subjects to select an exploratory gene and form its ordinary pointwise Wald interval;
- ordinary BH-FDR over thousands of same-sample post-selection p-values as the primary analysis.

## 11. Confirmatory layer freeze

Confirmatory inference uses a low-dimensional, externally defined immune-module list frozen before C3 association fitting.

Required module file:

```text
config/module_signatures.csv
```

with columns

```text
module_id,source,source_version,gene_symbol,inclusion_reason
```

The file must be nonempty, deduplicated, versioned, and checksum-frozen. At minimum, interferon-response modules may be included only with a public, reproducible source and exact gene membership. Any additional plasmablast, neutrophil, or other module requires the same standard.

Module scores use a frozen rule, normally the mean of training-standardised member genes when a prespecified minimum member proportion is observed. Module definitions cannot be changed after inspecting C3 coefficients.

Confirmatory outputs are profile estimates and 95% intervals for prespecified module/covariate coordinates. Same-sample exploratory gene discoveries are not promoted to confirmatory modules.

## 12. Exploratory gene-level layer

The exploratory layer has two distinct tasks:

1. subject-held-out prediction and selection stability;
2. split-sample coordinate inference for a small number of genes selected on independent subjects.

For split-sample inference:

- one subject subset selects genes;
- a disjoint subset estimates and forms intervals;
- roles may be swapped across repeated predetermined splits;
- the aggregation rule is frozen before fitting;
- selection frequency and forced-inclusion status are reported;
- results are labelled exploratory split-sample inference.

No genome-wide confirmatory FDR claim is permitted without a separate simultaneous/post-selection theory and protocol.

## 13. Frozen validation design

Use repeated five-fold subject-level cross-validation with five repetitions. Save the fold table before model fitting.

Stratification may use source batch and a coarse baseline C3 category determined before fitting. It must not split visits independently or use future/outer-fold outcomes for tuning.

Within each outer training set:

- select `lambda_beta` by inner subject-level CV on pinball loss;
- use the frozen `h`, `Lambda`, and Dantzig multiplier grids;
- select the smallest feasible Dantzig tolerance on the frozen grid;
- select screening dimension within training subjects;
- fit all adaptive preprocessing only on training subjects.

New-subject prediction uses fixed effects and no unavailable latent subject effect. Conditional prediction from earlier visits is a separate, clearly labelled analysis.

## 14. Empirical comparators

Primary empirical comparisons:

1. pooled high-dimensional quantile Lasso;
2. faithful QGEE-SCAD where the adapter passes its benchmark gate;
3. low-dimensional LQMM for the clinical/module layer;
4. `PROFILE-DQR` random-intercept primary model;
5. `PROFILE-DQR` intercept-plus-time-slope sensitivity.

Comparator target differences must be stated. A Gaussian Bayesian mixed shrinkage model, if retained, belongs only in a separate mean-prediction supplement and is not a quantile coverage comparator.

## 15. Method-implementation gate

Data download, construction, and audit may proceed in parallel with simulation development.

No confirmatory or full high-dimensional fit may begin until the exact implementation commit has:

- passed strict profile gradient/Hessian/Schur identity tests;
- passed required Dantzig tests;
- passed the simulation micro-preflight;
- passed required unit/package checks;
- recorded its commit and configuration checksum.

The empirical pipeline must save the implementation commit. Changing the estimator after fitting requires rerunning all empirical fits under a new versioned result directory.

## 16. Required pre-fit audit outputs

Before association fitting, produce and commit only non-sensitive summaries/manifests:

- source counts by cohort and replicate status;
- parsed metadata dictionary;
- reconciliation report;
- C3/C4/SLEDAI missingness table;
- eligible-subject and eligible-visit counts;
- cluster-size and follow-up-time distributions;
- treatment, race, sex, and batch frequencies;
- batch-by-cohort and batch-by-outcome summaries;
- technical-replicate agreement;
- expression-scale audit;
- probe/gene annotation flow;
- counts removed at each processing step;
- C3 outcome-gate report;
- frozen subject-fold table;
- checksum of `config/application.yml` and `config/module_signatures.csv`;
- package/session manifest.

The cohort-flow table must reconcile all source samples to final eligible visits.

## 17. Required final outputs

1. cohort/preprocessing flow table;
2. confirmatory module/covariate estimates and 95% intervals across primary quantiles;
3. subject-held-out pinball-loss comparison with paired uncertainty;
4. coefficient-by-quantile figure for confirmatory coordinates;
5. paired held-out loss-difference figure;
6. exploratory split-sample gene table with stability/selection frequencies;
7. optimisation, KKT, nuisance-gradient, Dantzig, and profile-identity diagnostics;
8. sensitivity table for SLEDAI adjustment, random slope, `tau=0.10`, filter size, and processing route;
9. failure and warning table;
10. complete provenance/session/config manifest.

All manuscript values are generated programmatically. No manual transcription or hand-edited result table is permitted.

## 18. Stop conditions

Stop and write a review report when:

- the C3 gate fails;
- identifiers cannot be reconciled;
- processed expression scales/batches cannot be diagnosed;
- technical replicates indicate unresolved normalization failure;
- the confirmatory module file is missing, empty, or changed after association fitting;
- strict profile geometry fails on the current implementation;
- fewer than 90% of primary proposed fits converge;
- Dantzig feasibility is below 90% for confirmatory coordinates;
- subject folds are materially imbalanced and cannot be repaired without leakage;
- a requested covariate encoding cannot be frozen without inspecting association results.

A stop condition does not authorise undisclosed outcome, preprocessing, model, quantile, or module switching.

## 19. Mandatory execution order

1. download official files and create the checksum manifest;
2. build expression/metadata objects and reconcile identifiers;
3. run the C3 and cohort audit;
4. freeze cohort rules, covariate encodings, probe annotation/collapse rules, module definitions, and subject folds;
5. confirm the method implementation commit passes its numerical gates;
6. run smoke fits on frozen training subsets;
7. run the full repeated subject-level validation;
8. run confirmatory module inference;
9. run exploratory split-sample gene analysis;
10. generate all tables/figures from immutable outputs.

Until Steps 1--5 pass, only download/build/audit work is authorised.
# Data directory

No subject-level data are committed to GitHub.

Expected local layout after running the application scripts:

```text
data/
  raw/GSE65391/
    GSE65391_series_matrix.txt.gz
    GSE65391_RAW.tar                         optional
    GSE65391_non-normalized_data_*.txt.gz   optional
    download_manifest.csv
  interim/GSE65391/
    gse65391_eset.rds
    expression_processed_probe.rds
    phenotype_parsed.csv
    probe_annotation.csv
  derived/GSE65391/
    outcome_audit.csv
    eligible_visits.csv
    subject_folds.csv
    analysis_cohort.rds
```

The download script records source URL, file size, timestamp and MD5 checksum. Raw and derived files are gitignored. Only aggregate, disclosure-safe manuscript tables may be committed after review.

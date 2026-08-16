# Results directory

Generated results are not committed by default. The expected local layout is:

```text
results/
  raw/simulation/              one file per experiment/replication chunk
  raw/application/             fold-level empirical fits
  intermediate/                merged validated objects
  tables/                      final CSV/LaTeX tables
  figures/                     final PDF/PNG figures
  logs/                        run, warning, failure and session logs
```

Every run must also copy the exact configuration file and save `sessionInfo()`.

Do not edit summary tables manually. `scripts/simulation/03_summarise.R` and `scripts/application/04_summarise_GSE65391.R` are the only supported routes from raw outputs to manuscript-ready tables.

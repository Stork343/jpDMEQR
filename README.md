# jpDMEQR

`jpDMEQR` provides screening, jointly penalised estimation, and debiased inference
for high-dimensional mixed-effects quantile regression.

## Install (GitHub)

```r
install.packages("remotes")
remotes::install_github("Stork343/jpDMEQR")
```

## Install (local)

```r
install.packages("/path/to/code/jpDMEQR", repos = NULL, type = "source")
```

## Quick start

```r
library(jpDMEQR)

dat <- generate_mixed_qr_data(n = 100, p = 200, scenario = "S1", tau = 0.5)
screen <- CA_IQR_SIS(dat$y, dat$X, dat$cluster_id, d_target = 40)
Xr <- dat$X[, screen, drop = FALSE]

fit <- JP_DME_QR(
  y = dat$y,
  X = Xr,
  Z = dat$Z,
  cluster_id = dat$cluster_id,
  tau = 0.5,
  h = dat$N^(-1 / 3),
  lambda1 = sqrt(log(ncol(Xr)) / dat$N),
  lambda2 = 1
)

inf <- Debias_Inference(
  y = dat$y,
  X = Xr,
  Z = dat$Z,
  cluster_id = dat$cluster_id,
  beta_hat = fit$beta,
  gamma_hat = fit$gamma,
  tau = 0.5,
  h = dat$N^(-1 / 3),
  lambda2 = 1,
  method = "RIDGE"
)
```

## Reproducibility scripts

- `inst/scripts/run_simulation_study.R`
- `inst/scripts/run_real_data_analysis.R`
- `inst/REPRODUCIBILITY.md`

These scripts are intended for end-to-end reproduction of manuscript outputs.

## Vignette

```r
vignette("jpDMEQR-workflow", package = "jpDMEQR")
```

## Citation

```r
citation("jpDMEQR")
```

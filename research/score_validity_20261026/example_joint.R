source("research/score_validity_20261026/load_joint.R")
dat <- generate_profile_qr_data_v2(n_clusters = 200L, p = 500L, s = 5L, tau = 0.5,
                                   q = 2L, seed = 20261026L)
fit <- ds_fit_profile(dat$y, dat$X, dat$Z, dat$cluster_id, tau = 0.5,
                       coordinates = c("x00001", "x00002", "x00006", "x00007"))
print(fit$table, row.names = FALSE)
print(fit$timing)

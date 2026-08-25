#!/usr/bin/env Rscript
# Round-6 ladder aggregation: variance identity, four-level ratios with 5,000
# paired bootstrap (seed 20260825), V1-V4 classification.
# Usage: Rscript scripts/simulation/07_ladder_aggregate_v6.R <run_id> <n_shards>
args <- commandArgs(trailingOnly = TRUE)
root <- "."
run_id <- args[1] %||% "pilot_v6_ladder_P06"
n_shards <- as.integer(args[2] %||% "8")
out_dir <- file.path(root, "results", "raw", "simulation", run_id)

rows <- lapply(seq_len(n_shards), function(s) {
  f <- file.path(out_dir, sprintf("ladder_%s_s%d.csv", run_id, s))
  if (file.exists(f)) read.csv(f, stringsAsFactors = FALSE) else NULL
})
rows <- rows[!vapply(rows, is.null, logical(1))]
tab <- do.call(rbind, rows)
tab <- tab[!is.na(tab$T) & is.finite(tab$T), ]
cat(sprintf("ladder rows: %d (reps %s)\n", nrow(tab),
            paste(sort(unique(tab$replicate)), collapse = ",")))

set.seed(20260825)
nboot <- 5000L
out <- lapply(sort(unique(tab$coordinate)), function(coord) {
  d <- tab[tab$coordinate == coord, ]
  B <- nrow(d)
  v <- function(x) stats::var(x) * (B - 1) / B   # MC variance (n denominator)
  VarT <- v(d$T); VarL <- v(d$L); VarQ <- v(d$Q); VarR <- v(d$R)
  CovLQ <- stats::cov(d$L, d$Q) * (B - 1) / B
  CovLR <- stats::cov(d$L, d$R) * (B - 1) / B
  CovQR <- stats::cov(d$Q, d$R) * (B - 1) / B
  ident <- VarL + VarQ + VarR + 2 * CovLQ + 2 * CovLR + 2 * CovQR
  ident_err <- abs(ident - VarT)
  boot_ratio <- function(scores) {
    # scores: list of per-rep vectors; ratio = mean(m1)/var(L+Q)
    m1 <- scores[[1]]; LQ <- d$L + d$Q
    replicate(nboot, {
      idx <- sample.int(B, B, replace = TRUE)
      mean(m1[idx]) / v(LQ[idx])
    })
  }
  bM_fit <- boot_ratio(list(d$V_fit_fit))
  bM_popfit <- boot_ratio(list(d$V_pop_fit))
  bM_poptgt <- boot_ratio(list(d$V_pop_target))
  bM_T <- replicate(nboot, {
    idx <- sample.int(B, B, replace = TRUE)
    v(d$T[idx]) / v((d$L + d$Q)[idx])
  })
  ci <- function(x) quantile(x, c(0.025, 0.975), names = FALSE)
  data.frame(
    coordinate = coord, B = B,
    VarT = VarT, VarL = VarL, VarQ = VarQ, VarR = VarR,
    CovLQ = CovLQ, CovLR = CovLR, CovQR = CovQR,
    identity_sum = ident, identity_error = ident_err,
    E_V_fit_fit = mean(d$V_fit_fit), E_V_fit_target = mean(d$V_fit_target),
    E_V_pop_fit = mean(d$V_pop_fit), E_V_pop_target = mean(d$V_pop_target),
    M_fit = mean(d$V_fit_fit) / VarT,   # placeholder, replaced below
    M_fit_ci_lo = ci(bM_fit)[1], M_fit_ci_hi = ci(bM_fit)[2],
    M_pop_fit = mean(d$V_pop_fit) / VarL,
    M_pop_fit_ci_lo = ci(bM_popfit)[1], M_pop_fit_ci_hi = ci(bM_popfit)[2],
    M_pop_target = mean(d$V_pop_target) / VarL,
    M_pop_target_ci_lo = ci(bM_poptgt)[1], M_pop_target_ci_hi = ci(bM_poptgt)[2],
    M_T = VarT / v(d$L + d$Q),
    M_T_ci_lo = ci(bM_T)[1], M_T_ci_hi = ci(bM_T)[2],
    stringsAsFactors = FALSE
  )
})
agg <- do.call(rbind, out)
# Fix M_fit: ratio of E(V_fit_fit) to Var(L+Q), bootstrap CI from bM_fit drawn
# with the LQ denominator (consistent with the amendment).
for (i in seq_len(nrow(agg))) {
  d <- tab[tab$coordinate == agg$coordinate[i], ]
  B <- nrow(d); v <- function(x) stats::var(x) * (B - 1) / B
  LQ <- d$L + d$Q
  bM <- replicate(nboot, { idx <- sample.int(B, B, replace = TRUE); mean(d$V_fit_fit[idx]) / v(LQ[idx]) })
  agg$M_fit[i] <- mean(d$V_fit_fit) / v(LQ)
  qq <- quantile(bM, c(0.025, 0.975), names = FALSE)
  agg$M_fit_ci_lo[i] <- qq[1]; agg$M_fit_ci_hi[i] <- qq[2]
}
print(agg, digits = 4)

# V1-V4 classification
cat("\n=== CLASSIFICATION (per amendment section 4) ===\n")
for (i in seq_len(nrow(agg))) {
  a <- agg[i, ]
  mfit_ok <- a$M_fit_ci_lo <= 1 && a$M_fit_ci_hi >= 1
  mT_gt <- a$M_T_ci_lo > 1
  mpf_ok <- a$M_pop_fit_ci_lo <= 1 && a$M_pop_fit_ci_hi >= 1
  mpt_ok <- a$M_pop_target_ci_lo <= 1 && a$M_pop_target_ci_hi >= 1
  cat(sprintf("%s: M_fit=%.3f [%.3f,%.3f] (CI contains 1: %s) | M_T=%.3f [%.3f,%.3f] (>1: %s) | identity err=%.2e\n",
              a$coordinate, a$M_fit, a$M_fit_ci_lo, a$M_fit_ci_hi, mfit_ok,
              a$M_T, a$M_T_ci_lo, a$M_T_ci_hi, mT_gt, a$identity_error))
}
cat("\nBranch: if M_fit CI contains 1 for all coords and M_T>1 -> V1; if M_fit<1 with pop ratios OK -> V2; if M_pop_fit<1 with M_pop_target OK -> V3; else V4.\n")
write.csv(agg, file.path(out_dir, "ladder_aggregate.csv"), row.names = FALSE)
cat("aggregate written to:", file.path(out_dir, "ladder_aggregate.csv"), "\n")
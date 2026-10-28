source("research/score_validity_20261026/load_joint.R")
args <- commandArgs(trailingOnly = TRUE)
worker <- as.integer(args[1L])
workers <- as.integer(args[2L])
registry <- read.csv("research/score_validity_20261026/evaluation.csv", stringsAsFactors = FALSE)
output <- "results/scratch/joint_profile_manuscript_20261026/ablation"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
request_file <- file.path(output, paste0("worker_", worker, "_request.rds"))
if (file.exists(request_file)) stop("This ablation worker already exists.")
saveRDS(list(registry = registry, worker = worker, workers = workers, seed_base = 1800000000L,
             variant = "outcome_and_profile_selection", started = Sys.time()), request_file)
current <- NULL
globalCallingHandlers(error = function(condition) {
  if (!is.null(current)) saveRDS(condition, file.path(current, "failure.rds"))
})
completed <- 0L
for (case in seq_len(nrow(registry))) {
  cfg <- registry[case, ]
  for (replicate in seq_len(cfg$replications)) {
    if ((replicate - 1L) %% workers != worker - 1L) next
    seed <- 1800000000L + 10000L * case + replicate
    current <- file.path(output, cfg$case_id, sprintf("replicate_%03d", replicate))
    if (dir.exists(current)) stop("An ablation unit already exists.")
    dir.create(current, recursive = TRUE)
    saveRDS(list(config = cfg, replicate = replicate, seed = seed), file.path(current, "request.rds"))
    dat <- generate_profile_qr_data_v2(n_clusters = cfg$n_clusters, p = cfg$p, s = cfg$s,
                                       tau = cfg$tau, q = cfg$q, error_dist = cfg$error_dist,
                                       signal = cfg$signal, rho_x = cfg$rho_x, seed = seed)
    coordinates <- c(1L, 2L, 6L, 7L)
    h <- cfg$n_clusters^(-.3)
    prepared <- ds_prepare(dat$y, dat$X, dat$Z, dat$cluster_id, cfg$tau, coordinates,
                            diag(cfg$q), rep(1, cfg$p), h, include_design = FALSE)
    proportions <- vapply(prepared$folds, function(fold) length(fold$validation_subjects), integer(1)) / prepared$n
    details <- lapply(seq_along(coordinates), function(j) {
      fits <- lapply(prepared$folds, function(fold) {
        model <- fold$models[[j]]
        ds_pair_small(model$data, model$initial, cfg$tau, h, diag(cfg$q), gram = model$gram)
      })
      estimate <- sum(proportions * vapply(fits, `[[`, numeric(1), "estimate"))
      se <- sqrt(sum(vapply(fits, function(fit) sum((fit$influence - mean(fit$influence))^2), numeric(1)))) / prepared$n
      shifts <- vapply(prepared$folds, function(fold) {
        support <- fold$graphs[[j]]$support
        omitted <- setdiff(dat$active, support)
        if (!length(omitted)) return(0)
        H <- cfg$rho_x^abs(outer(support, support, "-"))
        C <- cfg$rho_x^abs(outer(support, omitted, "-"))
        drop(solve(H, C %*% dat$beta0[omitted]))[1L]
      }, numeric(1))
      data.frame(coordinate = colnames(dat$X)[coordinates[j]], estimate = estimate, se = se,
                   lower = estimate - qnorm(.975) * se, upper = estimate + qnorm(.975) * se,
                   target_shift = sum(proportions * shifts),
                   max_fold_target_shift = max(abs(shifts)))
    })
    table <- do.call(rbind, details)
    table$truth <- dat$beta0[coordinates]
    table$covered <- table$lower <= table$truth & table$upper >= table$truth
    table$case_id <- cfg$case_id
    table$replicate <- replicate
    saveRDS(list(request = list(config = cfg, replicate = replicate, seed = seed), table = table,
                 selected_sets = lapply(prepared$folds, function(fold) lapply(fold$graphs, `[[`, "support"))),
             file.path(current, "result.rds"))
    completed <- completed + 1L
  }
  cat("Ablation worker", worker, "completed", cfg$case_id, "\n")
}
saveRDS(list(completed = completed, finished = Sys.time()), file.path(output, paste0("worker_", worker, "_complete.rds")))

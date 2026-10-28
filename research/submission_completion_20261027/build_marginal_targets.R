source("research/submission_completion_20261027/load.R")
library(statmod)
output <- "results/intermediate/joint_profile_submission_20261027/marginal_targets"
dir.create(output, recursive = TRUE, showWarnings = FALSE)

marginal_quadrature <- function(config, order) {
  rule <- statmod::gauss.quad(order, kind = "legendre")
  boundaries <- c(-10, -3, 3, 10)
  nodes <- weights <- numeric()
  for (j in 1:3) {
    half <- diff(boundaries[j + 0:1]) / 2
    centre <- mean(boundaries[j + 0:1])
    x <- centre + half * rule$nodes
    nodes <- c(nodes, x)
    weights <- c(weights, half * rule$weights * dnorm(x))
  }
  gaussian <- statmod::gauss.quad.prob(order * 2L, dist = "normal")
  list(x = nodes, weights = weights, b = config$sigma_b0 * gaussian$nodes,
    b_weights = gaussian$weights, sigma = exp(.35 * pmax(-3, pmin(3, nodes))),
    centring = error_quantile_function_v2(config$tau, config$error_dist, config$tau))
}

marginal_score <- function(theta, config, quadrature) {
  x <- quadrature$x
  sigma <- quadrature$sigma
  threshold <- theta[1L] + theta[2L] * x
  if (config$error_dist == "gaussian") {
    total_sd <- sqrt(sigma^2 + config$sigma_b0^2)
    argument <- (threshold + sigma * quadrature$centring) / total_sd
    probability <- pnorm(argument)
    density <- dnorm(argument) / total_sd
  } else if (config$error_dist == "t3") {
    argument <- sqrt(3) * (outer(threshold, quadrature$b, "-") / sigma + quadrature$centring)
    probability <- as.numeric(pt(argument, df = 3) %*% quadrature$b_weights)
    density <- as.numeric((sqrt(3) * dt(argument, df = 3) / sigma) %*% quadrature$b_weights)
  } else stop("This marginal integration requires Gaussian or t3 errors.")
  design <- cbind(1, x)
  score <- as.numeric(crossprod(design, quadrature$weights * (probability - config$tau)))
  curvature <- crossprod(design, design * (quadrature$weights * density))
  list(score = score, curvature = curvature)
}

marginal_root <- function(config, order) {
  quadrature <- marginal_quadrature(config, order)
  theta <- c(0, 0)
  for (iteration in 1:50) {
    evaluated <- marginal_score(theta, config, quadrature)
    if (max(abs(evaluated$score)) <= 1e-12) break
    step <- solve(evaluated$curvature, evaluated$score)
    accepted <- FALSE
    for (power in 0:30) {
      candidate <- theta - 2^(-power) * step
      candidate_score <- marginal_score(candidate, config, quadrature)$score
      if (sum(candidate_score^2) < sum(evaluated$score^2)) {
        theta <- candidate
        accepted <- TRUE
        break
      }
    }
    if (!accepted) stop("Marginal population Newton iteration failed its decrease condition.")
  }
  evaluated <- marginal_score(theta, config, quadrature)
  if (max(abs(evaluated$score)) > 1e-10) stop("Marginal population score tolerance was not attained.")
  list(theta = theta, score = evaluated$score, curvature = evaluated$curvature,
    order = order, iterations = iteration)
}

registry <- submission_registry("main")
records <- list()
for (id in c("B16", "B19")) {
  config <- config_row_to_list_v2(registry[registry$experiment_id == id, ])
  destination <- file.path(output, paste0(id, ".rds"))
  if (file.exists(destination)) {
    asset <- readRDS(destination)
    if (!identical(asset$config, config) || !asset$deterministic_target || asset$max_kkt_residual > 1e-8)
      stop("The saved marginal target does not match its registered design.")
    records[[id]] <- data.frame(experiment_id = id, intercept = asset$intercept_marginal,
      slope_shift = asset$slope_shift, maximum_score = asset$max_kkt_residual,
      refinement_difference = asset$quadrature_refinement_difference,
      maximum_standardised_validation_score = max(abs(unlist(lapply(asset$validation, `[[`, "standardised_score")))))
    next
  }
  if (config$design_type != "ar1" || config$q != 1L || config$x_b_corr != 0 ||
      config$informative_size || config$random_effect_dist != "normal" ||
      config$nonlinear_re_strength != 0 || !config$heteroskedastic)
    stop("The conditional Gaussian integration does not match this population design.")
  roots <- lapply(c(128L, 256L, 512L), function(order) marginal_root(config, order))
  difference <- max(abs(roots[[2L]]$theta - roots[[3L]]$theta))
  saveRDS(roots, file.path(output, paste0(id, "_quadrature.rds")))
  if (difference > 1e-8) stop("Marginal population quadrature refinement did not agree.")
  fitted <- roots[[3L]]
  beta <- make_beta_v2(config$p, config$s, config$signal)
  beta[1L] <- beta[1L] + fitted$theta[2L]
  # Gaussian AR1的条件均值使其余坐标的score由X1的score确定。
  validation_config <- config
  validation_config$p <- max(config$s + 6L, 14L)
  validation_config$n_clusters <- 100000L
  validation <- vector("list", 4L)
  for (replicate in 1:4) {
    seed <- seed_from_id_v2(20261027L, id, replicate, "marginal_population_validation")
    data <- simulate_from_config_v2(validation_config, seed)
    prediction <- fitted$theta[1L] + as.numeric(data$X %*% beta[seq_len(ncol(data$X))])
    score_rows <- cbind(intercept = 1, data$X) * (as.numeric(data$y <= prediction) - config$tau)
    score_clusters <- rowsum(score_rows, data$cluster_id, reorder = FALSE) /
      as.numeric(table(factor(data$cluster_id, levels = unique(data$cluster_id))))
    estimate <- colMeans(score_clusters)
    mcse <- apply(score_clusters, 2L, sd) / sqrt(nrow(score_clusters))
    validation[[replicate]] <- list(seed = seed, clusters = nrow(score_clusters), score = estimate,
      score_mcse = mcse, standardised_score = estimate / mcse,
      dimension = validation_config$p, full_dimension = config$p,
      reduction = "Gaussian_AR1_prefix_has_the_same_joint_distribution_with_all_active_coordinates")
    cat("Marginal score validation", id, replicate, "completed\n")
  }
  asset <- list(experiment_id = id, config = config, deterministic_target = TRUE,
    target_type = "marginal_linear_quantile_projection_with_unpenalised_intercept",
    beta_marginal = beta, intercept_marginal = fitted$theta[1L], slope_shift = fitted$theta[2L],
    max_kkt_residual = max(abs(fitted$score)), quadrature_refinement_difference = difference,
    target_mc_se_by_coordinate = setNames(rep(0, length(beta)), names(beta)),
    quadrature = roots, integration_interval = c(-10, 10),
    omitted_normal_probability = 2 * pnorm(-10), validation = validation,
    implementation_commit = current_commit_v2(), completed = Sys.time())
  saveRDS(asset, destination)
  records[[id]] <- data.frame(experiment_id = id, intercept = asset$intercept_marginal,
    slope_shift = asset$slope_shift, maximum_score = asset$max_kkt_residual,
    refinement_difference = difference,
    maximum_standardised_validation_score = max(abs(unlist(lapply(validation, `[[`, "standardised_score")))))
}
write.csv(do.call(rbind, records), file.path(output, "checks.csv"), row.names = FALSE)
print(do.call(rbind, records))

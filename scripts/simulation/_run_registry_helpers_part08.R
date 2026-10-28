validate_micro_registry_execution_v2 <- function(root, result, registry) {
  expected <- do.call(rbind, lapply(seq_len(nrow(registry)), function(i) {
    expand.grid(experiment_id = registry$experiment_id[i],
                replicate = seq_len(as.integer(registry$replications[i])),
                method_id = parse_pipe_character_v2(registry$methods[i]),
                stringsAsFactors = FALSE)
  }))
  keys <- c("experiment_id", "replicate", "method_id")
  make_key <- function(table) do.call(paste, c(table[keys], sep = "|"))
  method_complete <- nrow(result$metrics) == nrow(expected) &&
    !anyDuplicated(make_key(result$metrics)) &&
    setequal(make_key(result$metrics), make_key(expected))
  theory_complete <- nrow(result$theory) == nrow(expected) &&
    !anyDuplicated(make_key(result$theory)) &&
    setequal(make_key(result$theory), make_key(expected))
  coordinate_methods <- c(
    "PROFILE-DQR", "PROFILE-DQR-TRUE-SUPPORT", "PROFILE-DQR-TRUE-NUISANCE",
    "PROFILE-DQR-POP-H", "PROFILE-DQR-SPLIT", "POSTREFIT-EXACT-H",
    "SQR-DEBIASED-IID", "BIAS-ADJ-LQMM"
  )
  expected_coordinates <- expected[expected$method_id %in% coordinate_methods, , drop = FALSE]
  coordinate_counts <- as.integer(registry$target_coordinate_count[
    match(expected_coordinates$experiment_id, registry$experiment_id)])
  actual_coordinate_counts <- table(make_key(result$coordinates))
  required_coordinate_keys <- make_key(expected_coordinates)
  coordinate_complete <- nrow(result$coordinates) == sum(coordinate_counts) &&
    all(required_coordinate_keys %in% names(actual_coordinate_counts)) &&
    all(as.integer(actual_coordinate_counts[required_coordinate_keys]) == coordinate_counts) &&
    !anyDuplicated(result$coordinates[c(keys, "coordinate")])
  schema_complete <- identical(names(result$metrics), contract_schema_columns_v2(root, 2L)) &&
    identical(names(result$coordinates), contract_schema_columns_v2(root, 3L)) &&
    identical(names(result$theory), contract_schema_columns_v2(root, 4L)) &&
    identical(names(result$screening), contract_schema_columns_v2(root, 5L))
  failures <- result$metrics$status == "failed"
  known_failures <- all(!is.na(result$metrics$failure_stage[failures])) &&
    all(nzchar(result$metrics$failure_stage[failures])) &&
    !any(result$metrics$failure_stage[failures] == "unknown") &&
    all(!is.na(result$metrics$failure_message[failures])) &&
    all(nzchar(result$metrics$failure_message[failures]))
  failure_accounting <- validate_paired_failed_coordinate_accounting_v2(
    result$metrics, result$coordinates
  )
  expected_seeds <- vapply(seq_len(nrow(result$metrics)), function(i) {
    seed_from_id_v2(20260817L, result$metrics$experiment_id[i],
                    result$metrics$replicate[i], "training")
  }, integer(1))
  seeds_complete <- identical(as.integer(result$metrics$seed), expected_seeds)
  statuses_complete <- all(result$metrics$status %in% c("ok", "warning", "failed"))
  list(
    pass = length(result$errors) == 0L && method_complete && theory_complete &&
      coordinate_complete && schema_complete && known_failures &&
      isTRUE(failure_accounting) && seeds_complete && statuses_complete,
    method_complete = method_complete, theory_complete = theory_complete,
    coordinate_complete = coordinate_complete, schema_complete = schema_complete,
    known_failure_stages = known_failures, failed_coordinate_accounting = failure_accounting,
    seeds_complete = seeds_complete, expected_method_rows = nrow(expected),
    expected_coordinate_rows = sum(coordinate_counts),
    all_methods_successful = all(result$metrics$status %in% c("ok", "warning")),
    method_status_counts = as.list(table(result$metrics$status)),
    runner_failure_count = length(result$errors)
  )
}

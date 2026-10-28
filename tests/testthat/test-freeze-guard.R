root <- v2_require_root()
source(file.path(root, "scripts", "00_source_v2.R"), local = FALSE)
source_v2_module(root, "profile_v2", envir = .GlobalEnv)
source_v2_module(root, "simulation_v2", envir = .GlobalEnv)
sys.source(file.path(root, "scripts", "simulation", "_run_registry_helpers.R"),
           envir = environment())

test_that("the micro-preflight registry runs before a freeze manifest exists", {
  expect_false(requires_freeze_manifest_v2(
    file.path(root, "config", "simulation_preflight.csv")
  ))
  expect_false(requires_freeze_manifest_v2("simulation_preflight.csv"))
  expect_false(requires_freeze_manifest_v2(
    file.path("config", "simulation_preflight.csv")
  ))
})

test_that("main registries require a matching freeze manifest after the pilot", {
  expect_true(requires_freeze_manifest_v2(file.path(root, "config", "simulation_main.csv")))
  expect_true(requires_freeze_manifest_v2(file.path(root, "config", "simulation_main_b200.csv")))
  expect_false(requires_freeze_manifest_v2(file.path(root, "config", "simulation_pilot_v5.csv")))
  expect_true(requires_freeze_manifest_v2("simulation_main.csv"))
})

test_that("the registry file name decides the validation strength", {
  expect_identical(registry_type_for_path_v2("simulation_main.csv"), "main")
  expect_identical(registry_type_for_path_v2("config/simulation_main_b200.csv"), "main")
  expect_identical(registry_type_for_path_v2("config/simulation_pilot.csv"), "pilot")
  expect_identical(registry_type_for_path_v2("config/simulation_preflight.csv"), "preflight")
  expect_error(registry_type_for_path_v2("config/simulation_quick_repro.csv"),
               "Unknown registry file name")
})

test_that("the freeze manifest fingerprints every piece of gate evidence", {
  config_path <- file.path(root, "config", "simulation_preflight.csv")
  cfg <- read_simulation_registry_v2(config_path)
  evidence <- freeze_evidence_v2(root, cfg, config_path)
  expect_true(all(c("label", "path", "sha256") %in% names(evidence)))
  labels <- as.character(evidence$label)
  expect_true("config" %in% labels)
  expect_true("geometry_manifest" %in% labels)
  expect_true("micro_preflight_manifest" %in% labels)
  expect_true("pilot_gate_manifest" %in% labels)
  expect_setequal(labels[grepl("^benchmark:", labels)],
                  paste0("benchmark:", registry_methods_v2(cfg)))
  expect_setequal(labels[grepl("^profile_target:", labels)],
                  paste0("profile_target:", cfg$experiment_id[cfg$target_mode == "profile_mc"]))
  expect_setequal(labels[grepl("^population_direction:", labels)],
                  paste0("population_direction:", cfg$experiment_id[
                    vapply(cfg$methods, function(x) "PROFILE-DQR-POP-H" %in%
                             parse_pipe_character_v2(x), logical(1))]))
  expect_false(anyDuplicated(labels) > 0)
  # 记录中当前存在的依据文件必须带校验和，缺失的以 NA 记录，校验时同样视为不一致。
  expect_true(all(nzchar(evidence$sha256[file.exists(evidence$path)])))
  expect_true(all(is.na(evidence$sha256[!file.exists(evidence$path)])))
})

test_that("the freeze evidence check refuses changed or removed entries", {
  config_path <- file.path(root, "config", "simulation_preflight.csv")
  cfg <- read_simulation_registry_v2(config_path)
  evidence <- freeze_evidence_v2(root, cfg, config_path)
  expect_true(verify_freeze_evidence_v2(root, config_path, evidence))

  altered <- evidence
  altered$sha256[altered$label == "geometry_manifest"] <- paste0(rep("0", 64), collapse = "")
  expect_error(verify_freeze_evidence_v2(root, config_path, altered),
               "Freeze evidence changed since the manifest was written: geometry_manifest")

  trimmed <- evidence[evidence$label != "config", , drop = FALSE]
  expect_error(verify_freeze_evidence_v2(root, config_path, trimmed),
               "Freeze evidence changed since the manifest was written: config")
})

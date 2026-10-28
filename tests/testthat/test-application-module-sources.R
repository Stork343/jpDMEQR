root <- v2_require_root()
source(file.path(root, "scripts", "00_source_v2.R"), local = FALSE)
source_v2_module(root, "profile_v2", envir = .GlobalEnv)
source_v2_module(root, "application_v2", envir = .GlobalEnv)

testthat::test_that("公开module成员及官方Gene ID恢复保持完整来源", {
  rule <- yaml::read_yaml(file.path(root, "config", "module_scoring_20261022.yml"))
  signatures <- read_frozen_module_definitions_v2(file.path(root, rule$membership_file))
  members <- read_frozen_module_score_members_v2(file.path(root, rule$score_membership_file))
  testthat::expect_equal(nrow(signatures), 180L)
  testthat::expect_equal(nrow(members), 181L)
  testthat::expect_equal(as.integer(table(factor(members$module_id, levels = rule$module_order))),
                         c(27L, 53L, 18L, 59L, 24L))
  recovered <- signatures[signatures$source_gene_symbol %in% c("LOC129607", "LOC26010", "LOC728358"), ]
  testthat::expect_equal(recovered$gene_symbol[match(c("LOC129607", "LOC26010", "LOC728358"),
                                                    recovered$source_gene_symbol)],
                         c("CMPK2", "SPATS2L", "DEFA1B"))
  testthat::expect_true(all(recovered$mapping_status == "unique_official_gene_id"))
  testthat::expect_true(all(nzchar(recovered$source_mapping_url)))
  testthat::expect_true(all(members$member_type[members$module_id == "M5.15"] == "ILMN_probe"))
})

testthat::test_that("真实表达评分保持原始计量单位并仅使用training分布", {
  expression_path <- file.path(root, "data", "interim", "GSE65391", "expression_processed_probe.rds")
  annotation_path <- file.path(root, "data", "interim", "GSE65391", "probe_annotation_20261022.csv")
  testthat::skip_if_not(file.exists(expression_path) && file.exists(annotation_path),
                        "本地GSE65391表达和官方注释未提供。")
  rule <- yaml::read_yaml(file.path(root, "config", "module_scoring_20261022.yml"))
  signatures <- read_frozen_module_definitions_v2(file.path(root, rule$membership_file))
  members <- read_frozen_module_score_members_v2(file.path(root, rule$score_membership_file))
  annotation <- read.csv(annotation_path, stringsAsFactors = FALSE)
  expression <- readRDS(expression_path)
  testthat::expect_identical(rownames(expression), annotation$probe_id)
  source_dir <- file.path(root, "data", "raw", "public_modules", "banchereau2016_20261022")
  original <- read.csv(file.path(source_dir, "selected_probe_membership.csv"), stringsAsFactors = FALSE)
  counts <- vapply(rule$module_order, function(id) {
    sum(original$source_probe_id[original$module_id == id] %in% rownames(expression))
  }, integer(1))
  testthat::expect_equal(unname(counts), c(30L, 54L, 17L, 50L, 24L))
  testthat::expect_setequal(members$member_id[members$module_id == "M5.15"],
                            original$source_probe_id[original$module_id == "M5.15"])
  wanted <- annotation$probe_id[annotation$gene_symbol %in%
                                   members$member_id[members$member_type == "gene_symbol"]]
  wanted <- unique(c(wanted, members$member_id[members$member_type == "ILMN_probe"]))
  expression <- expression[rownames(expression) %in% wanted, seq_len(60L), drop = FALSE]
  train <- seq_len(40L)
  with_validation <- score_frozen_module_representations_v2(
    expression, annotation, signatures, members, train, rule
  )
  training_only <- score_frozen_module_representations_v2(
    expression[, train, drop = FALSE], annotation, signatures, members, train, rule
  )
  testthat::expect_identical(with_validation$center, training_only$center)
  testthat::expect_identical(with_validation$scale, training_only$scale)
  testthat::expect_identical(with_validation$gene_probe_map, training_only$gene_probe_map)
  for (id in rule$module_order) {
    rows <- members[members$module_id == id, ]
    probe_ids <- if (rows$member_type[1] == "ILMN_probe") rows$member_id else {
      with_validation$gene_probe_map$probe_id[match(rows$member_id,
                                                  with_validation$gene_probe_map$gene_symbol)]
    }
    direct <- colMeans(expression[probe_ids, , drop = FALSE])
    testthat::expect_equal(unname(with_validation$raw_scores[, id]), unname(direct), tolerance = 1e-12)
  }
  testthat::expect_true(all(with_validation$audit$eligible))
  testthat::expect_equal(colMeans(with_validation$scores[train, , drop = FALSE]),
                         setNames(rep(0, 5L), rule$module_order), tolerance = 1e-12)
})

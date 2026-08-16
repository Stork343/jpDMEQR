# GSE65391 ingestion, metadata parsing, audit and subject-split helpers.

clean_geo_name_v2 <- function(x) {
  x <- tolower(trimws(x))
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

standardize_geo_missing_v2 <- function(x) {
  x <- as.character(x)
  miss <- c("", "NA", "N/A", "Data Not Available", "Not Applicable",
            "not available", "not applicable", "null", "NULL")
  x[trimws(x) %in% miss] <- NA_character_
  x
}

split_characteristic_v2 <- function(value) {
  value <- standardize_geo_missing_v2(value)
  if (is.na(value)) return(c(key = NA_character_, value = NA_character_))
  loc <- regexpr(":", value, fixed = TRUE)
  if (loc[1] < 0) return(c(key = clean_geo_name_v2(value), value = NA_character_))
  key <- substr(value, 1, loc[1] - 1)
  val <- substr(value, loc[1] + 1, nchar(value))
  c(key = clean_geo_name_v2(key), value = trimws(val))
}

parse_geo_characteristics_v2 <- function(pdata) {
  pdata <- as.data.frame(pdata, stringsAsFactors = FALSE)
  char_cols <- grep("characteristics.*ch1", names(pdata), ignore.case = TRUE, value = TRUE)
  if (!length(char_cols)) stop("No GEO characteristics_ch1 columns found.")

  parsed_rows <- vector("list", nrow(pdata))
  duplicate_rows <- vector("list", nrow(pdata))
  all_keys <- character()

  for (i in seq_len(nrow(pdata))) {
    values <- unlist(pdata[i, char_cols, drop = FALSE], use.names = FALSE)
    pairs <- lapply(values, split_characteristic_v2)
    keys <- vapply(pairs, `[[`, character(1), "key")
    vals <- vapply(pairs, `[[`, character(1), "value")
    keep <- !is.na(keys) & nzchar(keys)
    keys <- keys[keep]; vals <- vals[keep]
    all_keys <- union(all_keys, keys)

    duplicated_keys <- unique(keys[duplicated(keys)])
    duplicate_rows[[i]] <- if (length(duplicated_keys)) {
      data.frame(sample_row = i, key = duplicated_keys, stringsAsFactors = FALSE)
    } else NULL

    row <- setNames(as.list(rep(NA_character_, length(unique(keys)))), unique(keys))
    for (key in unique(keys)) {
      candidates <- standardize_geo_missing_v2(vals[keys == key])
      candidates <- candidates[!is.na(candidates)]
      if (length(candidates)) row[[key]] <- candidates[length(candidates)]
    }
    parsed_rows[[i]] <- row
  }

  out <- as.data.frame(matrix(NA_character_, nrow(pdata), length(all_keys)),
                       stringsAsFactors = FALSE)
  names(out) <- all_keys
  for (i in seq_len(nrow(out))) {
    if (length(parsed_rows[[i]])) {
      out[i, names(parsed_rows[[i]])] <- parsed_rows[[i]]
    }
  }
  rownames(out) <- rownames(pdata)
  duplicates <- do.call(rbind, duplicate_rows[!vapply(duplicate_rows, is.null, logical(1))])
  list(parsed = out, duplicate_keys = duplicates, source_columns = char_cols)
}

coerce_numeric_fields_v2 <- function(metadata, fields) {
  for (field in intersect(fields, names(metadata))) {
    metadata[[field]] <- suppressWarnings(as.numeric(standardize_geo_missing_v2(metadata[[field]])))
  }
  metadata
}

build_gse65391_metadata_v2 <- function(eset) {
  if (!requireNamespace("Biobase", quietly = TRUE)) {
    stop("Package 'Biobase' is required.")
  }
  pdata <- Biobase::pData(eset)
  parsed <- parse_geo_characteristics_v2(pdata)
  meta <- cbind(
    data.frame(
      geo_accession = rownames(pdata),
      title = as.character(pdata$title %||% rownames(pdata)),
      source_name = as.character(pdata$source_name_ch1 %||% NA_character_),
      stringsAsFactors = FALSE,
      row.names = rownames(pdata)
    ),
    parsed$parsed
  )

  numeric_fields <- c(
    "visit", "visit_count", "cumulative_time",
    "days_since_diagnosis", "days_since_last_visit", "age", "wbc",
    "neutrophil_count", "lymphocyte_count", "monocyte_count",
    "platelet_count", "esr", "hgb", "hct", "cr", "alb", "ds_dna",
    "c3", "c4", "ast", "alt", "ldh", "sledai", "disease_activity",
    "batch"
  )
  meta <- coerce_numeric_fields_v2(meta, numeric_fields)

  disease_source <- meta$disease_state %||% meta$disease %||% meta$title
  disease_source <- toupper(as.character(disease_source))
  meta$group_type <- ifelse(
    grepl("SLE", disease_source) | grepl("WHOLE BLOOD-SLE", toupper(meta$title)),
    "SLE",
    ifelse(grepl("HEALTH", disease_source) | grepl("HEALTH", toupper(meta$title)),
           "HEALTHY", "UNRESOLVED")
  )

  # GEO uses subject strings such as SLE-40. Preserve the source string and create
  # a stable character subject identifier rather than coercing it to an integer.
  subject_source <- if ("subject" %in% names(meta)) {
    as.character(meta$subject)
  } else {
    rep(NA_character_, nrow(meta))
  }
  subject_source <- standardize_geo_missing_v2(subject_source)
  missing_subject <- is.na(subject_source) | !nzchar(subject_source)
  title_subject <- sub(".*-(SLE-[0-9]+)-.*", "\\1", meta$title)
  subject_source[missing_subject & grepl("SLE-[0-9]+", meta$title)] <-
    title_subject[missing_subject & grepl("SLE-[0-9]+", meta$title)]
  meta$subject_id <- subject_source

  batch_rep <- if ("batch_replicate" %in% names(meta)) {
    toupper(as.character(meta$batch_replicate))
  } else {
    rep("FALSE", nrow(meta))
  }
  meta$batch_replicate_flag <- batch_rep %in% c("TRUE", "1", "YES")
  meta$technical_replicate_flag <- meta$group_type == "HEALTHY" &
    (meta$batch_replicate_flag | grepl("HEALTHY-2", toupper(meta$title)))

  visit_number <- meta$visit %||% suppressWarnings(as.numeric(sub(".*-V([0-9]+)-.*", "\\1", meta$title)))
  meta$visit_number <- suppressWarnings(as.numeric(visit_number))
  meta$visit_id <- paste(meta$subject_id, sprintf("V%03d", meta$visit_number), sep = "_")
  meta$sample_id <- meta$geo_accession

  list(
    metadata = meta,
    duplicate_characteristic_keys = parsed$duplicate_keys,
    source_pdata = pdata
  )
}

build_gse65391_objects_v2 <- function(eset) {
  if (!requireNamespace("Biobase", quietly = TRUE)) stop("Biobase is required.")
  expression <- Biobase::exprs(eset)
  storage.mode(expression) <- "double"
  if (any(!is.finite(expression))) stop("Expression matrix contains non-finite values.")
  meta_obj <- build_gse65391_metadata_v2(eset)
  metadata <- meta_obj$metadata

  missing_expr <- setdiff(colnames(expression), rownames(metadata))
  missing_meta <- setdiff(rownames(metadata), colnames(expression))
  reconciliation <- data.frame(
    type = c(rep("expression_without_metadata", length(missing_expr)),
             rep("metadata_without_expression", length(missing_meta))),
    id = c(missing_expr, missing_meta),
    stringsAsFactors = FALSE
  )
  if (nrow(reconciliation)) {
    stop("Expression/metadata identifiers do not reconcile; inspect the saved report.")
  }
  metadata <- metadata[colnames(expression), , drop = FALSE]

  list(
    expression = expression,
    metadata = metadata,
    duplicate_characteristic_keys = meta_obj$duplicate_characteristic_keys,
    reconciliation = reconciliation,
    source_pdata = meta_obj$source_pdata
  )
}

safe_fraction_v2 <- function(num, den) if (den > 0) num / den else NA_real_

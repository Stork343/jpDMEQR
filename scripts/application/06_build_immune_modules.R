library(GEOquery)

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1] else "inspect-platform"
root <- normalizePath(".", mustWork = TRUE)
source_dir <- file.path(root, "data", "raw", "public_modules", "banchereau2016_20261022")
archive <- file.path(source_dir, "TcGSA_s007.zip")
source_objects <- new.env(parent = emptyenv())
load(unz(archive, "ReproducibleRFiles/GMTs_PLOScb.RData"), envir = source_objects)
modules <- source_objects$gmt_modulesV2
if (length(modules$genesets) != 260L || length(modules$geneset.names) != 260L) {
  stop("The public source does not contain the expected G2 repertoire.")
}
all_members <- do.call(rbind, lapply(seq_along(modules$genesets), function(i) {
  data.frame(module_id = modules$geneset.names[i],
             module_label = modules$geneset.descriptions[i],
             source_probe_id = as.character(modules$genesets[[i]]),
             stringsAsFactors = FALSE)
}))
selected_ids <- c("M1.2", "M3.4", "M4.11", "M5.12", "M5.15")
selected <- all_members[all_members$module_id %in% selected_ids, , drop = FALSE]
if (mode == "export-probes") {
  write.csv(all_members, file.path(source_dir, "g2_full_probe_membership.csv"), row.names = FALSE)
  write.csv(selected, file.path(source_dir, "selected_probe_membership.csv"), row.names = FALSE)
  print(as.data.frame(table(selected$module_id)))
} else if (mode == "freeze-scores") {
  genes <- read.csv(file.path(root, "config", "module_signatures.csv"), stringsAsFactors = FALSE,
                    colClasses = c(source_gene_id = "character", source_mapping_url = "character"))
  gene_coverage <- read.csv(file.path(source_dir, "module_mapping_coverage.csv"), stringsAsFactors = FALSE)
  probe_coverage <- read.csv(file.path(source_dir, "probe_level_expression_coverage.csv"), stringsAsFactors = FALSE)
  gene_ids <- c("M1.2", "M3.4", "M4.11", "M5.12")
  if (!all(gene_coverage$full_membership_available[match(gene_ids, gene_coverage$module_id)]) ||
      !probe_coverage$full_membership_available[match("M5.15", probe_coverage$module_id)]) {
    stop("The prespecified module representation lacks full measured membership.")
  }
  score_members <- do.call(rbind, lapply(c(gene_ids, "M5.15"), function(id) {
    if (id %in% gene_ids) {
      rows <- genes[genes$module_id == id, ]
      member_id <- rows$gene_symbol
      source_member_id <- rows$source_gene_symbol
      source_gene_id <- rows$source_gene_id
      source_mapping_url <- rows$source_mapping_url
      member_type <- "gene_symbol"
      source <- rows$source
      source_version <- rows$source_version
    } else {
      rows <- selected[selected$module_id == id, ]
      member_id <- source_member_id <- rows$source_probe_id
      source_gene_id <- source_mapping_url <- ""
      member_type <- "ILMN_probe"
      source <- "https://journals.plos.org/ploscompbiol/article/file?type=supplementary&id=10.1371/journal.pcbi.1004310.s007"
      source_version <- "TcGSA_s007_GMTs_PLOScb_RData_20141020_G2"
    }
    data.frame(module_id = id, member_type = member_type, member_id = member_id,
               source_member_id = source_member_id, weight = 1 / length(member_id),
               source_gene_id = source_gene_id, source_mapping_url = source_mapping_url,
               source = source, source_version = source_version,
               outcome_association_inspected = FALSE, stringsAsFactors = FALSE)
  }))
  target <- file.path(root, "config", "module_score_members_20261022.csv")
  if (file.exists(target)) stop("The versioned module-score membership already exists.")
  write.csv(score_members, target, row.names = FALSE)
  write.csv(data.frame(md5 = unname(tools::md5sum(target)),
                       source = "G2_public_gene_and_original_probe_memberships",
                       frozen_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
                       outcome_association_inspected = FALSE,
                       score = "equal_weight_mean_processed_normalised_expression"),
            paste0(target, ".freeze.csv"), row.names = FALSE)
  print(as.data.frame(table(score_members$module_id, score_members$member_type)))
} else if (mode == "inspect-platform") {
  integrity <- system2("gzip", c("-t", shQuote(file.path(source_dir, "GPL10558.annot.gz"))))
  if (integrity != 0L) stop("Official GPL10558 annotation download is incomplete.")
  gpl <- GEOquery::getGEO(filename = file.path(source_dir, "GPL10558.annot.gz"),
                          AnnotGPL = TRUE, destdir = source_dir)
  platform <- GEOquery::Table(gpl)
  print(dim(platform))
  print(names(platform))
  print(utils::head(platform, 2))
  print(GEOquery::Meta(gpl))
  write.csv(platform, file.path(source_dir, "GPL10558_official_annotation.csv"), row.names = FALSE)
} else if (mode == "inspect-submitter") {
  gpl <- GEOquery::getGEO(filename = file.path(source_dir, "GPL10558_submitter_complete.soft"),
                          destdir = source_dir)
  platform <- GEOquery::Table(gpl)
  print(dim(platform))
  print(names(platform))
  print(utils::head(platform, 2))
  if (nrow(platform) != 48107L) stop("The original GPL10558 table is incomplete.")
  write.csv(platform, file.path(source_dir, "GPL10558_submitter_complete_annotation.csv"), row.names = FALSE)
} else if (mode == "audit-expression-probes") {
  expression <- readRDS(file.path(root, "data", "interim", "GSE65391",
                                  "expression_processed_probe.rds"))
  if (is.null(rownames(expression)) || anyDuplicated(rownames(expression))) {
    stop("Expression probe identifiers must be unique.")
  }
  selected$present_in_expression <- selected$source_probe_id %in% rownames(expression)
  selected$source_weight <- 1 / as.numeric(table(selected$module_id)[selected$module_id])
  coverage <- do.call(rbind, lapply(selected_ids, function(id) {
    rows <- selected[selected$module_id == id, ]
    data.frame(module_id = id, module_label = rows$module_label[1],
               source_probe_count = nrow(rows), observed_probe_count = sum(rows$present_in_expression),
               missing_probe_count = sum(!rows$present_in_expression),
               full_membership_available = all(rows$present_in_expression),
               score_unit = "ILMN_probe", source_weight = 1 / nrow(rows),
               expression_probe_count = nrow(expression), expression_sample_count = ncol(expression),
               stringsAsFactors = FALSE)
  }))
  write.csv(selected, file.path(source_dir, "probe_level_expression_members.csv"), row.names = FALSE)
  write.csv(coverage, file.path(source_dir, "probe_level_expression_coverage.csv"), row.names = FALSE)
  print(coverage)
} else if (mode == "inspect-tmod") {
  tmod_archive <- file.path(source_dir, "tmod_0.50.13.tar.gz")
  if (system2("gzip", c("-t", shQuote(tmod_archive))) != 0L) {
    stop("Public tmod archive is incomplete.")
  }
  utils::untar(tmod_archive, files = c("tmod/data/tmod.rda", "tmod/DESCRIPTION"), exdir = source_dir)
  load(file.path(source_dir, "tmod", "data", "tmod.rda"), envir = source_objects)
  print(ls(source_objects))
  str(source_objects$tmod, max.level = 1)
  print(utils::head(source_objects$tmod$gs))
  print(utils::head(source_objects$tmod$gv))
} else if (mode %in% c("build", "check", "refresh-geneid")) {
  source(file.path(root, "scripts", "00_source_v2.R"))
  source_v2_module(root, "application_v2", envir = environment())
  tmod_archive <- file.path(source_dir, "tmod_0.50.13.tar.gz")
  if (system2("gzip", c("-t", shQuote(tmod_archive))) != 0L ||
      system2("gzip", c("-t", shQuote(file.path(source_dir, "GPL10558.annot.gz")))) != 0L) {
    stop("A public module or platform source is incomplete.")
  }
  utils::untar(tmod_archive, files = "tmod/data/tmod.rda", exdir = source_dir)
  load(file.path(source_dir, "tmod", "data", "tmod.rda"), envir = source_objects)
  tmod <- source_objects$tmod
  gpl <- GEOquery::getGEO(filename = file.path(source_dir, "GPL10558.annot.gz"),
                          AnnotGPL = TRUE, destdir = source_dir)
  platform <- GEOquery::Table(gpl)
  submitter <- GEOquery::getGEO(filename = file.path(source_dir, "GPL10558_submitter_complete.soft"),
                                destdir = source_dir)
  vendor <- GEOquery::Table(submitter)
  if (nrow(vendor) != nrow(platform) || !setequal(vendor$ID, platform$ID)) {
    stop("Original submitter and curated GPL10558 probe inventories do not reconcile.")
  }
  expression <- readRDS(file.path(root, "data", "interim", "GSE65391",
                                  "expression_processed_probe.rds"))
  current <- data.frame(probe_id = rownames(expression), stringsAsFactors = FALSE)
  if (anyDuplicated(current$probe_id) || !nrow(current)) stop("Invalid expression probe inventory.")
  annotation <- assemble_gpl10558_probe_annotation_v2(
    current$probe_id, platform, vendor,
    "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL10nnn/GPL10558/annot/GPL10558.annot.gz",
    as.character(GEOquery::Meta(gpl)$date)
  )
  public_members <- do.call(rbind, lapply(selected_ids, function(id) {
    index <- which(tmod$gs$SourceID == "DC" & tmod$gs$original.ID == id)
    if (length(index) != 1L) stop("The public G2 gene-level module identity is not unique.")
    data.frame(module_id = id, source_module_id = tmod$gs$ID[index],
               module_label = tmod$gs$Title[index],
               source_gene_symbol = as.character(tmod$gv[tmod$gs2gv[[index]]]),
               module_original_url = tmod$gs$URL[index], stringsAsFactors = FALSE)
  }))
  aliases <- strsplit(ifelse(is.na(annotation$submitter_synonyms), "", annotation$submitter_synonyms),
                      ";", fixed = TRUE)
  aliases <- lapply(aliases, trimws)
  keep <- annotation$mapping_status == "unique_gene"
  gene_summary <- jsonlite::read_json(file.path(source_dir, "ncbi_gene_summary.json"))$result
  public_members$gene_symbol <- public_members$source_gene_symbol
  public_members$mapping_status <- "unmapped_member"
  public_members$mapped_probe_count <- 0L
  public_members$source_gene_id <- ""
  public_members$source_mapping_url <- ""
  for (i in seq_len(nrow(public_members))) {
    symbol <- public_members$source_gene_symbol[i]
    match_exact <- keep & !is.na(annotation$gene_symbol) & annotation$gene_symbol == symbol
    match_original <- keep & !is.na(annotation$submitter_symbol) & annotation$submitter_symbol == symbol
    match_vendor <- keep & !is.na(annotation$submitter_ilmn_gene) & annotation$submitter_ilmn_gene == symbol
    match_alias <- keep & vapply(aliases, function(a) symbol %in% a, logical(1))
    candidates <- if (any(match_exact)) which(match_exact) else if (any(match_original))
      which(match_original) else if (any(match_vendor)) which(match_vendor) else which(match_alias)
    if (!length(candidates) && grepl("^LOC[0-9]+$", symbol)) {
      gene_id <- sub("^LOC", "", symbol)
      official <- gene_summary[[gene_id]]
      if (!is.null(official)) {
        if (official$uid != gene_id || official$organism$taxid != 9606L ||
            official$nomenclaturestatus != "Official" || nzchar(official$status) ||
            nzchar(official$currentid)) stop("An official LOC Gene ID mapping is not current and unique.")
        candidates <- which(keep & !is.na(annotation$entrez_id) &
                              annotation$entrez_id == gene_id &
                              annotation$gene_symbol == official$nomenclaturesymbol)
        public_members$source_gene_id[i] <- gene_id
        public_members$source_mapping_url[i] <- paste0("https://www.ncbi.nlm.nih.gov/gene/", gene_id)
      }
    }
    canonical <- unique(annotation$gene_symbol[candidates])
    if (length(canonical) == 1L) {
      public_members$gene_symbol[i] <- canonical
      public_members$mapping_status[i] <- if (nzchar(public_members$source_gene_id[i])) {
        "unique_official_gene_id"
      } else if (canonical == symbol) "exact_symbol" else "unique_submitter_alias"
      public_members$mapped_probe_count[i] <- sum(keep & annotation$gene_symbol == canonical, na.rm = TRUE)
    } else if (length(canonical) > 1L) {
      public_members$mapping_status[i] <- "ambiguous_symbol"
    }
  }
  membership <- do.call(rbind, lapply(split(public_members,
                                            paste(public_members$module_id, public_members$gene_symbol, sep = "::")),
                                     function(rows) {
    data.frame(
      module_id = rows$module_id[1],
      source = "https://cran.r-project.org/src/contrib/tmod_0.50.13.tar.gz",
      source_version = "tmod_0.50.13_DC_G2_GPL10558_20160809_NCBI_Gene_20261022",
      gene_symbol = rows$gene_symbol[1],
      inclusion_reason = paste0("G2既定", rows$module_label[1], "功能分类"),
      source_gene_symbol = paste(sort(unique(rows$source_gene_symbol)), collapse = "|"),
      source_gene_id = paste(sort(unique(rows$source_gene_id[nzchar(rows$source_gene_id)])), collapse = "|"),
      source_mapping_url = paste(sort(unique(rows$source_mapping_url[nzchar(rows$source_mapping_url)])), collapse = "|"),
      source_module_id = rows$source_module_id[1], module_label = rows$module_label[1],
      mapping_status = paste(sort(unique(rows$mapping_status)), collapse = "|"),
      stringsAsFactors = FALSE
    )
  }))
  membership <- membership[order(match(membership$module_id, selected_ids), membership$gene_symbol), ]
  rownames(membership) <- NULL
  coverage <- do.call(rbind, lapply(selected_ids, function(id) {
    source_rows <- public_members[public_members$module_id == id, ]
    fitted_rows <- membership[membership$module_id == id, ]
    data.frame(
      module_id = id, module_label = source_rows$module_label[1],
      public_gene_members = nrow(source_rows), frozen_gene_members = nrow(fitted_rows),
      mapped_public_members = sum(source_rows$mapped_probe_count > 0L),
      unmapped_public_members = sum(source_rows$mapping_status == "unmapped_member"),
      ambiguous_public_members = sum(source_rows$mapping_status == "ambiguous_symbol"),
      full_membership_available = all(source_rows$mapped_probe_count > 0L),
      original_probe_members = sum(selected$module_id == id),
      original_probes_in_expression = sum(selected$source_probe_id[selected$module_id == id] %in% current$probe_id),
      stringsAsFactors = FALSE
    )
  }))
  module_path <- file.path(root, "config", "module_signatures.csv")
  annotation_path <- file.path(root, "data", "interim", "GSE65391", "probe_annotation_20261022.csv")
  if (mode == "check") {
    actual <- read.csv(module_path, stringsAsFactors = FALSE,
                       colClasses = c(source_gene_id = "character", source_mapping_url = "character"))
    if (!identical(actual, membership)) stop("Frozen module CSV differs from the public-source rebuild.")
    read_frozen_module_definitions_v2(module_path)
  } else {
    if (mode == "refresh-geneid") {
      previous <- file.path(source_dir, "module_signatures_GPL10558_version.csv")
      previous_freeze <- paste0(previous, ".freeze.csv")
      if (file.exists(previous) || file.exists(previous_freeze)) {
        prior_freeze <- read.csv(previous_freeze, stringsAsFactors = FALSE)
        if (!file.exists(previous) ||
            unname(tools::md5sum(module_path)) != prior_freeze$md5 ||
            unname(tools::md5sum(previous)) != prior_freeze$md5) {
          stop("The preceding module version does not match the preserved freeze.")
        }
      } else {
        if (!file.copy(module_path, previous) ||
            !file.copy(paste0(module_path, ".freeze.csv"), previous_freeze)) {
          stop("The preceding module membership could not be preserved.")
        }
      }
      if (!identical(lapply(read.csv(annotation_path, stringsAsFactors = FALSE), as.character),
                     lapply(annotation, as.character))) {
        stop("The versioned probe annotation differs from the official-source rebuild.")
      }
    } else if (file.exists(module_path) || file.exists(annotation_path)) {
      stop("Frozen module or versioned annotation output already exists.")
    }
    if (mode == "build") write.csv(annotation, annotation_path, row.names = FALSE)
    write.csv(membership, module_path, row.names = FALSE)
    write.csv(data.frame(md5 = unname(tools::md5sum(module_path)),
                         source = "tmod_0.50.13_DC_G2_and_GPL10558_official_sources",
                         frozen_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
                         source_version = "G2_NCBI_Gene_20261022", outcome_association_inspected = FALSE,
                         full_membership_available = all(coverage$full_membership_available)),
              paste0(module_path, ".freeze.csv"), row.names = FALSE)
  }
  write.csv(public_members, file.path(source_dir, "public_member_mapping.csv"), row.names = FALSE)
  write.csv(coverage, file.path(source_dir, "module_mapping_coverage.csv"), row.names = FALSE)
  write.csv(as.data.frame(table(annotation$mapping_status)),
            file.path(source_dir, "annotation_mapping_flow.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(source_dir, "builder_sessionInfo.txt"))
  print(coverage)
} else {
  stop("Unknown module-builder mode.")
}

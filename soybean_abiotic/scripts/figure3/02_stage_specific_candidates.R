# Source: 图3.docx, extracted text lines 1143-1513
# 独立运行；阶段DEG表和全细胞类型DEG表
# See README.md and docs/KNOWN_LIMITATIONS.md before running.

###
# P2/P3 stage DEG + whole-cell-type DEG -> Endodermis-specific differential response candidates.
# Based only on significant DEG lists; absence from other cell types does not prove no expression.
# Standalone run; only openxlsx is required. No ComplexHeatmap, OrgDb, or RDS needed.
# install.packages("openxlsx")

PAPER_DIR <- file.path(Sys.getenv("SOYBEAN_PROJECT_ROOT", unset = "data"), "muscat_results_ck_ref_paper")
STAGE_DEG_FILE <- Sys.getenv("SOYBEAN_STAGE_DEG_FILE", unset = "")
if (!nzchar(STAGE_DEG_FILE)) stop("Set SOYBEAN_STAGE_DEG_FILE to the actual standardized stage DEG CSV.")
DEG_FILE_CANDIDATES <- c(
  file.path(PAPER_DIR, "All_DEG_filtered_paper_thresholds.csv"),
  file.path(PAPER_DIR, "All_DEG_filtered_combined.csv"),
  file.path(PAPER_DIR, "All_DEG_filtered_paper_thresholds.csv.gz"),
  file.path(PAPER_DIR, "All_DEG_filtered_combined.csv.gz"))
OUT_DIR <- file.path(PAPER_DIR, "Gene_heatmaps_ALKBH9B_network")
OUTPUT_XLSX <- file.path(OUT_DIR, "Endodermis_P2_P3_specific_DEG_candidates.xlsx")
PRINT_FILE <- file.path(OUT_DIR, "Endodermis_P2_P3_specific_DEG_print.txt")
TREATMENTS <- c("High_Temp", "Low_Temp", "PEG")
STAGES <- c("P2", "P3")
ENDO_LABELS <- c("Endodermis")  # If labels differ, add them here explicitly
FDR_CUT <- 0.05
LOG2FC_CUT <- log2(1.5)
DETECTION_CUT <- 0.10
USE_SAVED_DEG_FLAG <- TRUE
STRIP_GENE_VERSION <- TRUE # Only unify format/suffix, not coordinate conversion between reference genomes

read_csv <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  con <- if (grepl("[.]gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  utils::read.csv(con, check.names = FALSE, stringsAsFactors = FALSE)
}

col_pick <- function(d, candidates, required = TRUE) {
  hit <- intersect(candidates, names(d))
  if (length(hit)) return(hit[1])
  if (required) stop("Missing column(s): ", paste(candidates, collapse = "/"),
                     "; existing columns: ", paste(names(d), collapse = ", "))
  NULL
}

as_flag <- function(x) {
  z <- tolower(trimws(as.character(x)))
  if (any(!is.na(z) & !z %in% c("true", "false", "1", "0", ""))) {
    stop("Logical flag column contains unknown values.")
  }
  z %in% c("true", "1")
}

gene_key <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- sub("^GLYMA[._-]?", "GLYMA-", x)
  if (STRIP_GENE_VERSION) x <- sub("(V[0-9]+)?([.][0-9]+)*$", "", x)
  x
}

ct_key <- function(x) tolower(gsub("[[:space:]_.-]+", "", trimws(as.character(x))))

parse_group <- function(x) {
  z <- tolower(gsub("[[:space:]_.-]+", "", as.character(x)))
  ctl <- "(ckref|reference|control|ctrl|mock|ref|ck)"
  reverse <- grepl(paste0("^", ctl, "vs"), z)
  z <- sub(paste0("^", ctl, "vs"), "", z)
  z <- sub(paste0("vs", ctl, "$"), "", z)
  ans <- rep(NA_character_, length(z))
  ans[z %in% c("hightemp", "hightemperature", "heat", "ht", "hs")] <- "High_Temp"
  ans[z %in% c("lowtemp", "lowtemperature", "cold", "lt", "cs")] <- "Low_Temp"
  ans[z %in% c("peg", "pegtreatment", "drought", "osmotic")] <- "PEG"
  data.frame(treatment = ans, multiplier = ifelse(reverse, -1, 1))
}

parse_stage <- function(x) {
  z <- toupper(trimws(as.character(x)))
  z[z %in% c("S2", "MIDDLE")] <- "P2"
  z[z %in% c("S3", "LATE")] <- "P3"
  z
}

get_fdr <- function(d) {
  for (co in intersect(c("FDR", "p_adj.loc", "p_val_adj", "padj", "FDR_global", "p_adj.glb"), names(d))) {
    z <- suppressWarnings(as.numeric(d[[co]]))
    if (any(is.finite(z))) {
      if (any(is.finite(z) & (z < 0 | z > 1))) stop("FDR values not in 0~1: ", co)
      return(list(value = z, column = co))
    }
  }
  list(value = rep(NA_real_, nrow(d)), column = "Unavailable")
}

standardize <- function(raw, source, is_stage) {
  gc <- col_pick(raw, c("gene", "Gene", "gene_id", "GENE", "gene_raw"))
  tc <- col_pick(raw, c("treatment", "comparison", "contrast", "Comparison", "group"))
  fc <- col_pick(raw, c("logFC", "log2FC", "avg_log2FC", "log2FoldChange"), is_stage)
  cc <- col_pick(raw, c("cell_type", "celltype", "CellType", "cluster_id", "cluster"), FALSE)
  tr <- parse_group(raw[[tc]])
  cmp <- col_pick(raw, c("contrast", "comparison", "Comparison"), FALSE)
  if (!is.null(cmp) && cmp != tc) {
    cp <- parse_group(raw[[cmp]])
    both <- !is.na(cp$treatment) & !is.na(tr$treatment)
    if (any(cp$treatment[both] != tr$treatment[both])) stop(source, ": treatment/contrast inconsistent.")
    tr$multiplier[both] <- cp$multiplier[both]
  }
  bad <- grepl("temp|heat|cold|peg", as.character(raw[[tc]]), ignore.case = TRUE) & is.na(tr$treatment)
  if (any(bad)) stop(source, ": cannot recognize treatment labels: ",
                     paste(unique(raw[[tc]][bad]), collapse = ", "),
                     "; supported: highTemp_vs_REF / lowTemp_vs_REF / PEG_vs_REF, no automatic timepoint merging.")
  fdr <- get_fdr(raw)
  if (is_stage) {
    sc <- col_pick(raw, c("state", "stage", "period", "CK_reference_state"))
    stage <- parse_stage(raw[[sc]])
    # First file comes from known Endodermis analysis; cluster_id may be P1/P2/P3, not other cell types.
    celltype <- rep("Endodermis", nrow(raw))
    keep_ct <- rep(TRUE, nrow(raw))
    explicit_ct <- col_pick(raw, c("cell_type", "celltype", "CellType"), FALSE)
    if (!is.null(explicit_ct)) keep_ct <- ct_key(raw[[explicit_ct]]) %in% ct_key(ENDO_LABELS)
  } else {
    if (is.null(cc)) stop("Whole-cell DEG table missing cell_type/cluster_id column.")
    stage <- rep(NA_character_, nrow(raw))
    celltype <- as.character(raw[[cc]])
    keep_ct <- rep(TRUE, nrow(raw))
  }
  d <- data.frame(
    Gene = as.character(raw[[gc]]),
    gene_key = gene_key(raw[[gc]]),
    treatment = tr$treatment,
    stage = stage,
    celltype = celltype,
    log2FC = if (is.null(fc)) rep(NA_real_, nrow(raw)) else suppressWarnings(as.numeric(raw[[fc]])) * tr$multiplier,
    FDR = fdr$value,
    FDR_column = rep(fdr$column, nrow(raw)),
    source_file = rep(source, nrow(raw))
  )
  valid <- is.finite(d$log2FC) & is.finite(d$FDR)
  if ("tested" %in% names(raw)) valid <- valid & as_flag(raw$tested)
  if ("method" %in% names(raw)) valid <- valid & !grepl("descriptive", raw$method, ignore.case = TRUE)
  if (is_stage) {
    if (USE_SAVED_DEG_FLAG && "is_DEG" %in% names(raw)) {
      selected <- valid & as_flag(raw$is_DEG)
    } else {
      detection <- if ("max_detection" %in% names(raw)) suppressWarnings(as.numeric(raw$max_detection)) else rep(NA_real_, nrow(raw))
      if (any(is.finite(detection) & (detection < 0 | detection > 1))) stop("Detection rate must be 0~1.")
      selected <- valid & d$FDR <= FDR_CUT & abs(d$log2FC) >= LOG2FC_CUT &
        (!is.finite(detection) | detection >= DETECTION_CUT)
    }
  } else {
    # Second file is already filtered DEG: keep existing filtering criteria, no new threshold.
    selected <- if ("is_DEG" %in% names(raw)) as_flag(raw$is_DEG) else rep(TRUE, nrow(raw))
    if ("tested" %in% names(raw)) selected <- selected & as_flag(raw$tested)
    if ("method" %in% names(raw)) selected <- selected & !grepl("descriptive", raw$method, ignore.case = TRUE)
  }
  d$valid_test <- valid
  d$selected_DEG <- selected
  keep <- d$treatment %in% TREATMENTS & keep_ct
  if (is_stage) keep <- keep & d$stage %in% STAGES
  d <- d[keep, , drop = FALSE]
  if (!nrow(d)) stop(source, ": no data for selected treatments/stages.")
  if (anyNA(d$gene_key) || any(!nzchar(d$gene_key)) || anyNA(d$celltype) || any(!nzchar(trimws(d$celltype)))) {
    stop(source, ": empty gene or cell type.")
  }
  d$celltype_key <- ct_key(d$celltype)
  keys <- if (is_stage) c("gene_key", "treatment", "stage") else c("gene_key", "treatment", "celltype_key")
  if (anyDuplicated(d[keys])) {
    stop(source, ": duplicate records for same gene/comparison or version conflict; please check, no automatic min FDR selection.")
  }
  d$direction <- ifelse(!is.finite(d$log2FC), "Unknown",
                        ifelse(d$log2FC > 0, "Up", ifelse(d$log2FC < 0, "Down", "Zero")))
  if (is_stage) {
    for (s in STAGES) {
      for (tr in TREATMENTS) {
        if (!any(d$valid_test[d$stage == s & d$treatment == tr])) {
          stop(s, " / ", tr, ": no valid statistical test; cannot treat as no DEG.")
        }
      }
    }
  }
  d
}

collapse <- function(x) paste(sort(unique(as.character(x))), collapse = ";")

check_candidates <- function(stage, allcells) {
  cell_deg <- allcells[allcells$selected_DEG, , drop = FALSE]
  if (!any(cell_deg$celltype_key %in% ct_key(ENDO_LABELS))) {
    stop("Whole-cell DEG table has no Endodermis label, please check ENDO_LABELS.")
  }
  cand <- stage[stage$selected_DEG, , drop = FALSE]
  cand$Endodermis_DEG_in_celltype_table <- rep(FALSE, nrow(cand))
  cand$Endodermis_log2FC <- cand$Endodermis_FDR <- cand$Endodermis_direction <- rep("", nrow(cand))
  cand$Other_celltypes_with_DEG <- rep("", nrow(cand))
  cand$N_other_celltypes_with_DEG <- integer(nrow(cand))
  cand$Other_celltype_records_available <- rep(FALSE, nrow(cand))
  for (tr in TREATMENTS) {
    e <- cell_deg[cell_deg$treatment == tr & cell_deg$celltype_key %in% ct_key(ENDO_LABELS), , drop = FALSE]
    o <- cell_deg[cell_deg$treatment == tr & !cell_deg$celltype_key %in% ct_key(ENDO_LABELS), , drop = FALSE]
    e_groups <- split(e, e$gene_key)
    o_groups <- split(o$celltype, o$gene_key)
    for (i in which(cand$treatment == tr)) {
      ee <- e_groups[[cand$gene_key[i]]]
      oo <- o_groups[[cand$gene_key[i]]]
      cand$Endodermis_DEG_in_celltype_table[i] <- !is.null(ee) && nrow(ee) > 0
      cand$Endodermis_log2FC[i] <- collapse(ee$log2FC)
      cand$Endodermis_FDR[i] <- collapse(ee$FDR)
      cand$Endodermis_direction[i] <- collapse(ee$direction)
      cand$Other_celltypes_with_DEG[i] <- collapse(oo)
      cand$N_other_celltypes_with_DEG[i] <- length(unique(ct_key(oo)))
      cand$Other_celltype_records_available[i] <- nrow(o) > 0
    }
    if (!nrow(o)) {
      warning(tr, ": no DEG records in other cell types; specificity remains NA, cannot declare specificity.")
    }
  }
  cand$Endodermis_only_DEG_same_stress <- ifelse(
    cand$Other_celltype_records_available,
    cand$Endodermis_DEG_in_celltype_table & cand$N_other_celltypes_with_DEG == 0,
    NA
  )
  cand$Direction_matches_whole_Endodermis <- ifelse(
    cand$Endodermis_direction %in% c("Up", "Down") & cand$direction %in% c("Up", "Down"),
    cand$Endodermis_direction == cand$direction,
    NA
  )
  others <- cell_deg[!cell_deg$celltype_key %in% ct_key(ENDO_LABELS), , drop = FALSE]
  other_any <- unique(others$gene_key)
  covered <- all(TREATMENTS %in% others$treatment)
  cand$No_other_celltype_DEG_across_all3 <- if (covered) !cand$gene_key %in% other_any else rep(NA, nrow(cand))
  cand$Check_result <- ifelse(
    !cand$Other_celltype_records_available, "Cannot_assess_other_celltypes",
    ifelse(!cand$Endodermis_DEG_in_celltype_table, "Absent_from_whole_Endodermis_DEG_list",
           ifelse(cand$N_other_celltypes_with_DEG > 0, "Also_DEG_in_other_celltypes", "Endodermis_only_DEG_candidate"))
  )
  cand
}

make_gene_table <- function(rows) {
  keys <- unique(rows$gene_key)
  data.frame(
    Gene = vapply(keys, function(k) rows$Gene[match(k, rows$gene_key)], character(1)),
    gene_key = keys,
    Supported_stresses = vapply(keys, function(k) collapse(rows$treatment[rows$gene_key == k]), character(1)),
    N_stresses = vapply(keys, function(k) length(unique(rows$treatment[rows$gene_key == k])), integer(1)),
    stringsAsFactors = FALSE
  )
}

build_results <- function(audit) {
  sheets <- list()
  count <- list()
  for (s in STAGES) {
    candidates <- audit[audit$stage == s, , drop = FALSE]
    chosen <- candidates[candidates$Endodermis_only_DEG_same_stress %in% TRUE, , drop = FALSE]
    for (tr in TREATMENTS) {
      nm <- paste0(s, "_", c(High_Temp = "Heat", Low_Temp = "Cold", PEG = "PEG")[[tr]])
      sheets[[nm]] <- chosen[chosen$treatment == tr, , drop = FALSE]
      count[[nm]] <- data.frame(
        stage = s,
        treatment = tr,
        N_stage_DEG = sum(candidates$treatment == tr),
        N_selected = nrow(sheets[[nm]]),
        N_unassessable = sum(candidates$treatment == tr & is.na(candidates$Endodermis_only_DEG_same_stress))
      )
    }
    any_genes <- make_gene_table(chosen)
    common <- any_genes[any_genes$N_stresses == 3, , drop = FALSE]
    sheets[[paste0(s, "_any_stress_genes")]] <- any_genes
    sheets[[paste0(s, "_common3_genes")]] <- common
    sheets[[paste0(s, "_common3_details")]] <- chosen[chosen$gene_key %in% common$gene_key, , drop = FALSE]
    # Stricter candidate: at least one treatment passes, and gene absent from any other cell type DEG list across all three stresses.
    strict <- chosen[chosen$No_other_celltype_DEG_across_all3 %in% TRUE, , drop = FALSE]
    sheets[[paste0(s, "_strict_any_genes")]] <- make_gene_table(strict)
    sheets[[paste0(s, "_all_candidate_audit")]] <- candidates
  }
  shared <- intersect(sheets$P2_common3_genes$gene_key, sheets$P3_common3_genes$gene_key)
  sheets$P2_P3_shared_common3 <- sheets$P2_common3_genes[sheets$P2_common3_genes$gene_key %in% shared, , drop = FALSE]
  sheets <- c(list(Counts = do.call(rbind, count)), sheets)
  sheets
}

write_excel <- function(sheets) {
  wb <- openxlsx::createWorkbook()
  header <- openxlsx::createStyle(fgFill = "#21918C", fontColour = "white", textDecoration = "bold", wrapText = TRUE)
  for (nm in names(sheets)) {
    d <- sheets[[nm]]
    openxlsx::addWorksheet(wb, nm)
    openxlsx::writeData(wb, nm, d, headerStyle = header, withFilter = nrow(d) > 0,
                        keepNA = TRUE, na.string = "NA")
    openxlsx::freezePane(wb, nm, firstRow = TRUE, firstCol = TRUE)
    openxlsx::setColWidths(wb, nm, cols = seq_len(ncol(d)), widths = 23)
    openxlsx::setRowHeights(wb, nm, rows = 1, heights = 45)
    if (nm == "README") {
      openxlsx::setColWidths(wb, nm, cols = 2, widths = 110)
      openxlsx::addStyle(wb, nm, openxlsx::createStyle(wrapText = TRUE), rows = 2:(nrow(d)+1), cols = 2)
      openxlsx::setRowHeights(wb, nm, rows = 2:(nrow(d)+1), heights = 50)
    }
  }
  openxlsx::saveWorkbook(wb, OUTPUT_XLSX, overwrite = TRUE)
}

print_results <- function(sheets) {
  # Block print to avoid max.print truncation; also used to save full text report.
  show_all <- function(d) {
    if (!nrow(d)) {
      print(d, row.names = FALSE)
      return(invisible(NULL))
    }
    for (start in seq.int(1L, nrow(d), by = 100L)) {
      print(d[start:min(start + 99L, nrow(d)), , drop = FALSE], row.names = FALSE, max = 10000)
    }
  }
  show_all(sheets$Counts)
  for (s in STAGES) {
    for (tr in c("Heat", "Cold", "PEG")) {
      nm <- paste0(s, "_", tr)
      cat("\n", nm, ": Endodermis specific DEG candidates\n", sep = "")
      show_all(sheets[[nm]][c("Gene", "log2FC", "FDR", "direction", "Direction_matches_whole_Endodermis")])
    }
    for (suffix in c("common3_genes", "strict_any_genes")) {
      nm <- paste0(s, "_", suffix)
      cat("\n", nm, "\n", sep = "")
      show_all(sheets[[nm]])
    }
  }
}

main <- function() {
  if (!requireNamespace("openxlsx", quietly = TRUE)) stop("Please install: install.packages('openxlsx')")
  existing <- DEG_FILE_CANDIDATES[file.exists(DEG_FILE_CANDIDATES)]
  if (!length(existing)) stop("Cannot find whole-cell DEG file.")
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
  message("Stage DEG: ", STAGE_DEG_FILE, "\nWhole-cell DEG: ", existing[1])
  stage <- standardize(read_csv(STAGE_DEG_FILE), STAGE_DEG_FILE, TRUE)
  cells <- standardize(read_csv(existing[1]), existing[1], FALSE)
  audit <- check_candidates(stage, cells)
  sheets <- build_results(audit)
  covered <- cells[cells$selected_DEG, , drop = FALSE]
  sheets$Celltype_DEG_records <- as.data.frame(table(
    Treatment = factor(covered$treatment, TREATMENTS),
    Celltype = covered$celltype
  ))
  sheets$Gene_ID_mapping <- unique(rbind(
    stage[c("Gene", "gene_key", "source_file")],
    cells[c("Gene", "gene_key", "source_file")]
  ))
  notes <- data.frame(
    Item = c(
      "Stage DEG file", "Whole-cell DEG file", "Selection rule", "any_stress_genes",
      "common3_genes", "strict_any_genes", "Direction", "Scope", "NA",
      "Stage interpretation", "Background scope"
    ),
    Explanation = c(
      STAGE_DEG_FILE,
      existing[1],
      "DEG in a given stage/treatment, also Endodermis DEG in whole-cell table, and absent from other cell type DEG lists for the same treatment.",
      "At least one treatment in this stage passes the above filter; each gene listed once.",
      "Genes passing the filter in all three treatments (High_Temp, Low_Temp, PEG) for this stage.",
      "At least one treatment passes, and gene absent from any other cell type DEG list across all three stresses.",
      "Default does not require matching logFC direction between stage and whole Endodermis; Direction_matches_whole_Endodermis can be used for further filtering.",
      "Can only be called Endodermis differential response candidates. Absence from significant DEG list does not prove no expression in other cells, nor is it a direct statistical test between cell types.",
      "When other cell types have no DEG records, specificity remains NA; lack of records may also come from low-expression filtering, untested, or non-significant.",
      "P2/P3 are filtered separately; P2 genes are not required to be absent in P3, and vice versa.",
      "Other cell types exclusion covers High_Temp, Low_Temp, PEG; not ABA, IAA, or other treatments."
    )
  )
  sheets <- c(list(README = notes), sheets)
  write_excel(sheets)
  printed <- capture.output(print_results(sheets))
  writeLines(printed, PRINT_FILE, useBytes = TRUE)
  cat(paste(printed, collapse = "\n"), "\n")
  message("Excel: ", OUTPUT_XLSX, "\nFull print report: ", PRINT_FILE)
  invisible(sheets)
}

if (!isTRUE(getOption("p23_specificity.skip_main"))) main()

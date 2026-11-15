#!/usr/bin/env Rscript
# Run on the server that can read MERGED_RDS, or source this file in server RStudio.
# This exports real values only. It does not fit UMAP or compute DE from single cells.

MERGED_RDS <- "/data01/liuxy/Project/CK2_added_v3_Results/seurat_obj_with_CK2_v3.rds"
OUT_DIR <- "/data01/liuxy/Project/soybean_abiotic_exports"
GENES <- c("Glyma.05G103600", "Glyma.15G040200", "Glyma.20G176100", "Glyma.08G163500", "Glyma.11G147700")
EXPORT_GROUPS <- c("CK", "Heat", "Cold", "PEG") # Set c("CK") for CK-only marker plots.
ASSAY <- "RNA"
GROUP_COL <- NULL    # e.g. "Group1" or "condition"; set if inference is ambiguous.
CELLTYPE_COL <- NULL # e.g. "celltype_final" or "cell_type".
UMAP_REDUCTION <- NULL # Prefers umap, then ref.umap, then refumap.
DF_COLUMN <- NULL    # Optional existing classification; only confirmed singlets retained.
# Add exact sample aliases only after checking the experimental sample sheet.
# SAMPLE_GROUP_MAP <- c("your_control_sample"="CK", "your_heat_sample"="Heat")
SAMPLE_GROUP_MAP <- setNames(character(), character())
# FEATURE_MAP <- c("Glyma.05G103600"="Cortex-1") # Only if object rownames use aliases.
FEATURE_MAP <- setNames(character(), character())
# All genes are searchable by cell-type/treatment summaries. Files are loaded
# in small gene batches by the website rather than as a whole transcriptome.
EXPORT_ALL_GENE_SUMMARIES <- TRUE
SUMMARY_GENES_PER_FILE <- 200L
# Optional full single-cell UMAP vectors, stored sparsely. Large transcriptomes
# may exceed static hosting budgets; Shiny can query them without JSON export.
EXPORT_ALL_GENE_UMAP <- TRUE
MAX_SINGLE_CELL_EXPORT_MB <- 750

required <- c("Seurat", "SeuratObject", "Matrix", "jsonlite", "digest")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install missing R packages first: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages(library(Seurat))
if (!file.exists(MERGED_RDS)) stop("RDS not accessible from this R session: ", MERGED_RDS)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(OUT_DIR)) stop("Cannot create output directory: ", OUT_DIR)
write_csv <- function(x, filename) utils::write.csv(x, file.path(OUT_DIR, filename), row.names = FALSE, fileEncoding = "UTF-8")
key <- function(x) gsub("[^a-z0-9]", "", tolower(trimws(as.character(x))))
aliases <- list(
  CK = c("CK", "CK1", "CK2", "CK_1", "CK_2", "control", "control1", "control2", "untreated"),
  Heat = c("Heat", "Heat1", "Heat2", "highTemp", "highTemp1", "highTemp2", "37C", "8.37temp"),
  Cold = c("Cold", "Cold1", "Cold2", "lowTemp", "lowTemp1", "lowTemp2", "6C", "9.6temp"),
  PEG = c("PEG", "PEG1", "PEG2")
)
lookup <- setNames(rep(names(aliases), lengths(aliases)), key(unlist(aliases, use.names = FALSE)))
if (length(SAMPLE_GROUP_MAP)) {
  if (is.null(names(SAMPLE_GROUP_MAP)) || any(!nzchar(names(SAMPLE_GROUP_MAP))) || any(!SAMPLE_GROUP_MAP %in% names(aliases))) stop("SAMPLE_GROUP_MAP must have named sample IDs and values CK/Heat/Cold/PEG.")
  lookup[key(names(SAMPLE_GROUP_MAP))] <- unname(SAMPLE_GROUP_MAP)
}
map_groups <- function(x) unname(lookup[key(x)])

obj <- readRDS(MERGED_RDS)
if (!inherits(obj, "Seurat")) stop("RDS must contain a Seurat object.")
if (!ASSAY %in% names(obj@assays)) stop("Assay not found: ", ASSAY)
md <- obj[[]]
write_csv(data.frame(column = names(md)), "metadata_columns.csv")
candidates <- intersect(c("PublicationGroup", "group_id", "Group1", "Group", "group", "condition", "treatment", "Treatment", "orig.ident", "sample", "sample_id", "sample_name"), names(md))
audit <- do.call(rbind, lapply(candidates, function(nm) {
  tt <- table(as.character(md[[nm]]), useNA = "ifany")
  data.frame(column = nm, value = names(tt), cells = as.integer(tt), mapped_group = map_groups(names(tt)))
}))
if (!is.null(audit)) write_csv(audit, "group_metadata_audit.csv")
if (!is.null(GROUP_COL)) {
  if (!GROUP_COL %in% names(md)) stop("GROUP_COL is not present. Inspect metadata_columns.csv.")
  groups <- map_groups(md[[GROUP_COL]]); group_source <- GROUP_COL
} else {
  groups <- rep(NA_character_, nrow(md)); used <- character()
  for (nm in candidates) {
    z <- map_groups(md[[nm]])
    conflict <- !is.na(groups) & !is.na(z) & groups != z
    if (any(conflict)) stop("Conflicting treatment labels in metadata. Inspect group_metadata_audit.csv and set GROUP_COL / SAMPLE_GROUP_MAP explicitly.")
    take <- is.na(groups) & !is.na(z)
    if (any(take)) {groups[take] <- z[take]; used <- c(used, nm)}
  }
  group_source <- paste(used, collapse = ", ")
}
if (!any(groups == "CK", na.rm = TRUE)) stop("No CK labels identified. Inspect group_metadata_audit.csv; set GROUP_COL or exact SAMPLE_GROUP_MAP. REF is not assumed to be CK.")
keep <- !is.na(groups) & groups %in% EXPORT_GROUPS
if (!is.null(DF_COLUMN)) {
  if (!DF_COLUMN %in% names(md)) stop("DF_COLUMN is missing.")
  df <- tolower(trimws(as.character(md[[DF_COLUMN]])))
  if (any(!is.na(df) & !df %in% c("singlet", "doublet", ""))) stop("Unexpected DoubletFinder labels.")
  keep <- keep & !is.na(df) & df == "singlet"
}
if (!any(keep)) stop("No cells remain after treatment / optional singlet selection.")
selected_groups <- setNames(groups[keep], rownames(md)[keep])
write_csv(data.frame(group = c(names(aliases), "Unmapped_or_other"), cells = c(vapply(names(aliases), function(g) sum(groups == g, na.rm = TRUE), integer(1)), sum(is.na(groups)))), "input_group_counts.csv")
obj <- subset(obj, cells = rownames(md)[keep]); DefaultAssay(obj) <- ASSAY
groups <- unname(selected_groups[colnames(obj)])
if (!any(groups == "CK")) stop("No CK cells remain after filtering.")
if (ncol(obj) > 500000L) stop("More than 500,000 selected cells. Export fewer treatment groups for the browser viewer.")

if (is.null(CELLTYPE_COL)) {
  cols <- intersect(c("PublicationCellType", "celltype_final", "cell_type_final", "cell_type", "celltype", "CellType", "predicted.cell_type", "predicted.celltype", "annotation", "predicted.id"), names(obj[[]]))
  cols <- cols[vapply(cols, function(nm) any(!is.na(obj[[]][[nm]]) & nzchar(trimws(as.character(obj[[]][[nm]])))), logical(1))]
  if (!length(cols)) stop("Cell-type annotations not found. Set CELLTYPE_COL.")
  CELLTYPE_COL <- cols[1]
}
if (!CELLTYPE_COL %in% names(obj[[]])) stop("CELLTYPE_COL not found.")
ct <- trimws(as.character(obj[[]][[CELLTYPE_COL]])); ct[is.na(ct) | !nzchar(ct)] <- "Unannotated"
if (is.null(UMAP_REDUCTION)) {
  reductions <- intersect(c("umap", "ref.umap", "refumap"), names(obj@reductions))
  if (!length(reductions)) stop("Existing UMAP not found. Set UMAP_REDUCTION.")
  UMAP_REDUCTION <- reductions[1]
}
if (!UMAP_REDUCTION %in% names(obj@reductions)) stop("Configured UMAP reduction is missing.")
emb <- Embeddings(obj, reduction = UMAP_REDUCTION)
if (ncol(emb) < 2L || !all(colnames(obj) %in% rownames(emb))) stop("UMAP/cell alignment is invalid.")
emb <- emb[colnames(obj), 1:2, drop = FALSE]
if (any(!is.finite(emb))) stop("UMAP has non-finite coordinates.")
if (inherits(obj[[ASSAY]], "Assay5")) {
  obj <- JoinLayers(obj, assay = ASSAY)
  normalize_here <- !"data" %in% Layers(obj[[ASSAY]])
  if (normalize_here) obj <- NormalizeData(obj, assay = ASSAY, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  mat <- GetAssayData(obj, assay = ASSAY, layer = "data")
} else {
  mat <- GetAssayData(obj, assay = ASSAY, slot = "data")
  normalize_here <- !nrow(mat) || !ncol(mat)
  if (normalize_here) {
    obj <- NormalizeData(obj, assay = ASSAY, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
    mat <- GetAssayData(obj, assay = ASSAY, slot = "data")
  }
}
expression_label <- if (normalize_here) "RNA LogNormalize: log1p(CP10K)" else "Existing RNA normalized data; original normalization retained"
feature_key <- function(x) tolower(sub("\\.[0-9]+$", "", x))
features <- setNames(vapply(GENES, function(g) {
  if (g %in% names(FEATURE_MAP)) {
    f <- unname(FEATURE_MAP[g]); if (!f %in% rownames(mat)) stop("Mapped feature missing: ", f); return(f)
  }
  if (g %in% rownames(mat)) return(g)
  hit <- rownames(mat)[feature_key(rownames(mat)) == feature_key(g)]
  if (length(hit) > 1L) stop("Ambiguous feature versions for ", g, "; set FEATURE_MAP.")
  if (!length(hit)) return(NA_character_)
  hit
}, character(1)), GENES)
write_csv(data.frame(gene = GENES, object_feature = unname(features)), "gene_match_report.csv")
if (anyNA(features)) warning("Missing genes: ", paste(names(features)[is.na(features)], collapse = ", "), ". See gene_match_report.csv.")
GENES <- names(features)[!is.na(features)]
features <- features[GENES]
if (!length(GENES)) stop("None of the requested marker genes exists in the RNA assay. Inspect gene_match_report.csv.")
expr <- as.matrix(mat[unname(features), colnames(obj), drop = FALSE]); rownames(expr) <- GENES
if (any(!is.finite(expr)) || any(expr < 0)) stop("RNA expression must be finite and nonnegative; scaled residuals are not supported.")
cells <- data.frame(id = colnames(obj), x = unname(emb[, 1]), y = unname(emb[, 2]), celltype = ct, group = groups, stringsAsFactors = FALSE)
expression <- setNames(lapply(GENES, function(g) I(as.numeric(expr[g, ]))), GENES)
payload <- list(schema_version = 1L, source_label = basename(MERGED_RDS), created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), expression_label = expression_label, group_source = group_source, celltype_column = CELLTYPE_COL, reduction = UMAP_REDUCTION, cells = cells, expression = expression)
json_path <- file.path(OUT_DIR, "atlas.json")
jsonlite::write_json(payload, json_path, dataframe = "rows", auto_unbox = TRUE, digits = 7, pretty = FALSE)
if (isTRUE(EXPORT_ALL_GENE_SUMMARIES)) {
  message("Computing all-gene expression summaries by treatment and cell type...")
  mat <- mat[, cells$id, drop = FALSE]
  # Sparse aggregation retains every selected cell; no sampling or imputation.
  strata <- unique(cells[, c("group", "celltype"), drop = FALSE])
  strata <- strata[order(match(strata$group, EXPORT_GROUPS), strata$celltype), , drop = FALSE]
  strata_key <- paste(strata$group, strata$celltype, sep = "\t")
  membership <- match(paste(cells$group, cells$celltype, sep = "\t"), strata_key)
  strata$n_cells <- tabulate(membership, nbins = nrow(strata))
  design <- Matrix::sparseMatrix(i = seq_len(nrow(cells)), j = membership, x = 1, dims = c(nrow(cells), nrow(strata)))
  means <- sweep(as.matrix(mat %*% design), 2, strata$n_cells, "/")
  percents <- 100 * sweep(as.matrix((mat > 0) %*% design), 2, strata$n_cells, "/")
  summary_ids <- rownames(mat)
  for (g in names(features)) summary_ids[match(features[g], rownames(mat))] <- g
  if (anyDuplicated(summary_ids)) stop("Summary gene IDs are ambiguous after alias mapping.")
  cell_order_hash <- digest::digest(cells$id, algo = "xxhash64")
  # Both artifacts receive the same signature to prevent mismatched cell order.
  payload$cell_order_hash <- cell_order_hash
  jsonlite::write_json(payload, json_path, dataframe = "rows", auto_unbox = TRUE, digits = 7, pretty = FALSE)
  nonzero_count <- Matrix::nnzero(mat)
  estimated_mb <- nonzero_count * 26 / 1024^2
  include_umap <- isTRUE(EXPORT_ALL_GENE_UMAP) && estimated_mb <= MAX_SINGLE_CELL_EXPORT_MB
  if (isTRUE(EXPORT_ALL_GENE_UMAP) && !include_umap) warning("Estimated sparse single-cell export exceeds the configured budget (", round(estimated_mb), " MB). Exporting all-gene summaries and the five selected UMAP vectors; use Shiny for arbitrary-gene UMAP.")
  batch_size <- if (include_umap) min(SUMMARY_GENES_PER_FILE, 20L) else SUMMARY_GENES_PER_FILE
  batch_dir <- file.path(OUT_DIR, "gene_summaries")
  dir.create(batch_dir, showWarnings = FALSE)
  batches <- split(seq_along(summary_ids), ceiling(seq_along(summary_ids) / batch_size))
  gene_files <- setNames(vector("list", length(summary_ids)), summary_ids)
  for (b in seq_along(batches)) {
    ids <- batches[[b]]
    relative <- sprintf("gene_summaries/%05d.json", b)
    records <- setNames(lapply(ids, function(i) {
      record <- list(mean = I(as.numeric(means[i, ])), pct = I(as.numeric(percents[i, ])))
      if (include_umap) {
        sparse_row <- summary(methods::as(mat[i, , drop = FALSE], "TsparseMatrix"))
        positive <- sparse_row$x > 0
        record$umap <- list(i = I(as.integer(sparse_row$j[positive] - 1L)), v = I(as.numeric(sparse_row$x[positive])))
      }
      record
    }), summary_ids[ids])
    jsonlite::write_json(list(schema_version = 1L, genes = records), file.path(OUT_DIR, relative), auto_unbox = TRUE, digits = 7, pretty = FALSE)
    gene_files[ids] <- rep(list(relative), length(ids))
  }
  index <- list(schema_version = 1L, type = "gene_summary_index", source_label = basename(MERGED_RDS), expression_label = expression_label, cell_order_hash = cell_order_hash,
                strata = strata, gene_files = gene_files, umap_exported = include_umap, n_cells = nrow(cells), n_genes = length(summary_ids))
  jsonlite::write_json(index, file.path(OUT_DIR, "expression_index.json"), dataframe = "rows", auto_unbox = TRUE, pretty = FALSE)
  message("Wrote ", length(summary_ids), " gene summaries in ", length(batches), " files. Arbitrary-gene UMAP vectors exported: ", include_umap)
  message("For all-gene local browsing: download the entire export folder, then select Import expression folder on the website. For publication: put its contents under docs/.")
}
write_csv(data.frame(group = names(table(groups)), cells = as.integer(table(groups))), "export_group_counts.csv")
capture.output(sessionInfo(), file = file.path(OUT_DIR, "sessionInfo.txt"))
message("Exported ", nrow(cells), " cells, ", length(GENES), " genes. Group metadata: ", group_source, "; annotations: ", CELLTYPE_COL, "; UMAP: ", UMAP_REDUCTION)
message("JSON: ", json_path, " (", round(file.info(json_path)$size / 1024^2, 1), " MB)")
if (file.info(json_path)$size > 100 * 1024^2) warning("JSON exceeds the browser's 100 MB import limit. Export fewer groups/genes or use the Shiny app.")
message("Download atlas.json to your computer, open the website, then select Import atlas JSON. For shared automatic loading, publish it as docs/atlas.json in the GitHub repository.")

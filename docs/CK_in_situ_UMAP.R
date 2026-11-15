#!/usr/bin/env Rscript
# CK five-gene UMAP + browser atlas export.
# Run on the Linux server: Rscript CK_in_situ_UMAP.R

MERGED_RDS <- "/data01/liuxy/Project/CK2_added_v3_Results/seurat_obj_with_CK2_v3.rds"
OUT_DIR <- "/data01/liuxy/Project/soybean_abiotic_exports/CK_in_situ_UMAP"
GENES <- c("Glyma.05G103600", "Glyma.15G040200", "Glyma.12G234200",
           "Glyma.08G163500", "Glyma.10G022500")
ASSAY <- "RNA"
UMAP_REDUCTION <- NULL  # Auto-select umap, then ref.umap; set explicitly if necessary.
GROUP_COL <- NULL       # Set exact metadata column if automatic selection is ambiguous.
CELLTYPE_COL <- NULL    # e.g. cell_type, celltype_final, predicted.cell_type.
CK_VALUES <- c("CK", "CK1", "CK2", "CK_1", "CK_2", "Control", "Control1", "Control2")
# This is a multi-treatment merged object: CK must be identified in metadata.
# Never treat all merged cells as CK. REF is not automatically mapped to CK.
CK_ONLY_OBJECT <- FALSE
DF_COLUMN <- NULL       # Optional explicit DF.classifications column; no new QC is run.
POINT_SIZE <- 0.30
PLOT_ALPHA <- 1
RASTER_DPI <- 600L
SEED <- 42L
# Optional exact alias map if rownames use symbols instead of Glyma IDs:
# FEATURE_MAP <- c("Glyma.05G103600" = "Cortex-1")
FEATURE_MAP <- setNames(character(), character())

required <- c("Seurat", "SeuratObject", "Matrix", "ggplot2", "patchwork", "ggrastr", "jsonlite")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install required packages: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({library(Seurat); library(ggplot2); library(patchwork)})
if (!file.exists(MERGED_RDS)) stop("RDS not found: ", MERGED_RDS)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(OUT_DIR)) stop("Cannot create output directory: ", OUT_DIR)
set.seed(SEED)
options(ggrastr.default.dpi = RASTER_DPI)
write_csv <- function(x, name) utils::write.csv(x, file.path(OUT_DIR, name), row.names = FALSE, fileEncoding = "UTF-8")

merged <- readRDS(MERGED_RDS)
if (!inherits(merged, "Seurat")) stop("RDS must contain a Seurat object, not a list.")
if (!ASSAY %in% names(merged@assays)) stop("Assay not found: ", ASSAY)
md <- merged[[]]
# Inspect these files if automatic column detection stops.
write_csv(data.frame(column = colnames(md)), "metadata_columns.csv")
candidate_groups <- intersect(c("PublicationGroup", "Group1", "Group", "group", "condition", "treatment", "Treatment", "orig.ident", "sample", "sample_id"), colnames(md))
audit <- do.call(rbind, lapply(candidate_groups, function(nm) {
  tt <- table(as.character(md[[nm]]), useNA = "ifany")
  data.frame(column = nm, value = names(tt), cells = as.integer(tt))
}))
if (!is.null(audit)) write_csv(audit, "group_metadata_audit.csv")
norm <- function(x) tolower(trimws(as.character(x)))
ck_match <- function(x) !is.na(x) & norm(x) %in% norm(CK_VALUES)
if (!is.null(GROUP_COL)) {
  if (!GROUP_COL %in% colnames(md)) stop("GROUP_COL does not exist: ", GROUP_COL)
  keep <- ck_match(md[[GROUP_COL]])
  if (!any(keep)) stop("No CK values matched GROUP_COL. Inspect group_metadata_audit.csv and set CK_VALUES.")
  group_source <- GROUP_COL
} else {
  hits <- candidate_groups[vapply(candidate_groups, function(nm) any(ck_match(md[[nm]])), logical(1))]
  if (length(hits)) {
    masks <- lapply(hits, function(nm) ck_match(md[[nm]]))
    if (!all(vapply(masks, identical, logical(1), masks[[1]]))) {
      stop("CK selection differs across metadata columns: ", paste(hits, collapse = ", "), ". Set GROUP_COL and CK_VALUES explicitly.")
    }
    keep <- masks[[1]]; group_source <- hits[1]
  } else {
    if (!isTRUE(CK_ONLY_OBJECT)) stop("No CK metadata matched. Set GROUP_COL / CK_VALUES, or confirm CK_ONLY_OBJECT for a dedicated CK object.")
    # Refuse clearly labelled non-CK conditions even in CK-only fallback mode.
    nonck <- c("heat", "cold", "peg", "hightemp", "lowtemp", "aba", "iaa", "sa", "meja", "ebl")
    detected <- unlist(lapply(candidate_groups, function(nm) norm(md[[nm]])), use.names = FALSE)
    if (any(grepl("^(heat|cold|peg|hightemp|lowtemp|aba|iaa|sa|meja|ebl)([_.-]?[0-9]+)?$", detected), na.rm = TRUE)) stop("Non-CK treatments detected; disable CK_ONLY_OBJECT and set GROUP_COL / CK_VALUES.")
    keep <- rep(TRUE, nrow(md)); group_source <- "explicit CK_ONLY_OBJECT assumption"
    warning("No CK label matched; treating this dedicated CK reference object as CK because CK_ONLY_OBJECT=TRUE. Check group_metadata_audit.csv.")
  }
}
message("CK selection source: ", group_source, "; selected cells: ", sum(keep))

if (!is.null(DF_COLUMN)) {
  if (!DF_COLUMN %in% colnames(md)) stop("DF_COLUMN not found.")
  status <- norm(md[[DF_COLUMN]])
  if (any(!is.na(status) & !status %in% c("singlet", "doublet", ""))) stop("Unexpected DoubletFinder classification values.")
  # When explicitly requested, only retain confirmed singlets.
  keep <- keep & !is.na(status) & status == "singlet"
}
if (!any(keep)) stop("No CK cells remain.")
obj <- subset(merged, cells = rownames(md)[keep])
DefaultAssay(obj) <- ASSAY
rm(merged); invisible(gc())

if (is.null(CELLTYPE_COL)) {
  ct_candidates <- intersect(c("PublicationCellType", "celltype_final", "cell_type_final", "cell_type", "celltype", "CellType", "predicted.cell_type", "predicted.celltype", "annotation", "predicted.id"), colnames(obj[[]]))
  ct_candidates <- ct_candidates[vapply(ct_candidates, function(nm) any(!is.na(obj[[]][[nm]]) & nzchar(trimws(as.character(obj[[]][[nm]])))), logical(1))]
  if (!length(ct_candidates)) stop("No cell type column found; set CELLTYPE_COL. See metadata_columns.csv.")
  # Prefer the final annotation; log the chosen column for review.
  CELLTYPE_COL <- ct_candidates[1]
}
if (!CELLTYPE_COL %in% colnames(obj[[]])) stop("CELLTYPE_COL not found.")
celltypes <- trimws(as.character(obj[[]][[CELLTYPE_COL]]))
celltypes[is.na(celltypes) | !nzchar(celltypes)] <- "Unannotated"
message("Cell type column: ", CELLTYPE_COL)

if (is.null(UMAP_REDUCTION)) {
  reductions <- intersect(c("umap", "ref.umap", "refumap"), names(obj@reductions))
  if (!length(reductions)) stop("No existing UMAP found. Set UMAP_REDUCTION; this script does not concatenate or refit sample UMAPs.")
  UMAP_REDUCTION <- reductions[1]
}
if (!UMAP_REDUCTION %in% names(obj@reductions)) stop("UMAP reduction not found: ", UMAP_REDUCTION)
emb <- Embeddings(obj, reduction = UMAP_REDUCTION)
if (ncol(emb) < 2 || !all(colnames(obj) %in% rownames(emb))) stop("UMAP dimensions/cell names are invalid.")
emb <- emb[colnames(obj), 1:2, drop = FALSE]
if (any(!is.finite(emb))) stop("UMAP contains missing or infinite coordinates.")

# Join RNA v5 layers before reading data; do not use integrated or scale.data values.
if (inherits(obj[[ASSAY]], "Assay5")) {
  obj <- JoinLayers(obj, assay = ASSAY)
  if (!"data" %in% Layers(obj[[ASSAY]])) {
    if (!"counts" %in% Layers(obj[[ASSAY]])) stop("RNA counts/data layers missing.")
    obj <- NormalizeData(obj, assay = ASSAY, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
    expression_label <- "RNA log1p(counts / total counts * 10000); normalized by export script"
  } else expression_label <- "RNA data layer (existing normalized values; verify original normalization)"
  mat <- GetAssayData(obj, assay = ASSAY, layer = "data")
} else {
  mat <- GetAssayData(obj, assay = ASSAY, slot = "data")
  if (!nrow(mat) || !ncol(mat)) {
    obj <- NormalizeData(obj, assay = ASSAY, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
    mat <- GetAssayData(obj, assay = ASSAY, slot = "data")
    expression_label <- "RNA log1p(counts / total counts * 10000); normalized by export script"
  } else expression_label <- "RNA data slot (existing normalized values; verify original normalization)"
}
feature_key <- function(x) tolower(sub("\\.[0-9]+$", "", x))
features <- vapply(GENES, function(g) {
  if (g %in% names(FEATURE_MAP)) {
    f <- unname(FEATURE_MAP[g]); if (!f %in% rownames(mat)) stop("Mapped feature not found: ", f); return(f)
  }
  if (g %in% rownames(mat)) return(g)
  candidates <- rownames(mat)[feature_key(rownames(mat)) == feature_key(g)]
  if (length(candidates) == 1) return(candidates)
  if (length(candidates) > 1) stop("Ambiguous versioned feature: ", g, ". Set FEATURE_MAP.")
  NA_character_
}, character(1))
write_csv(data.frame(gene = GENES, object_feature = unname(features)), "gene_match_report.csv")
if (anyNA(features)) warning("Missing genes: ", paste(GENES[is.na(features)], collapse = ", "), ". See gene_match_report.csv; set exact FEATURE_MAP if necessary.")
available_genes <- GENES[!is.na(features)]
if (!length(available_genes)) stop("No requested marker genes found in the RNA assay.")
expr <- as.matrix(mat[unname(features[!is.na(features)]), colnames(obj), drop = FALSE])
rownames(expr) <- available_genes
if (any(!is.finite(expr)) || any(expr < 0)) stop("RNA data must be finite nonnegative normalized expression; do not use scaled residuals.")

plot_data <- data.frame(cell = colnames(obj), UMAP1 = emb[, 1], UMAP2 = emb[, 2], CellType = celltypes, Group = "CK", stringsAsFactors = FALSE)
base_theme <- theme_classic(base_size = 13, base_family = "Helvetica") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14), plot.subtitle = element_text(hjust = 0.5, size = 9), legend.title = element_text(size = 9), legend.text = element_text(size = 8), plot.margin = margin(5, 5, 5, 5))
save_plot <- function(p, name, width = 5, height = 4.5) {
  ggsave(file.path(OUT_DIR, paste0(name, ".pdf")), p, device = grDevices::pdf, width = width, height = height, units = "in", useDingbats = FALSE, family = "Helvetica", bg = "white")
  ggsave(file.path(OUT_DIR, paste0(name, ".png")), p, width = width, height = height, units = "in", dpi = RASTER_DPI, bg = "white")
}
ct_colors <- c("Cortex"="#2C80C8", "Dividing Cell"="#F58A17", "Endodermis"="#BD70A2", "Epidermis"="#91746A", "Procambium"="#098577", "Meristem"="#79AB39", "Pericycle"="#2B7E37", "Phloem"="#DE3035", "QC"="#BF2327", "Root Cap"="#8491B2", "Xylem"="#7F3298")
extra <- setdiff(unique(celltypes), names(ct_colors))
if (length(extra)) ct_colors <- c(ct_colors, setNames(grDevices::hcl.colors(length(extra), "Dark 3"), extra))
reference_plot <- ggplot(plot_data, aes(UMAP1, UMAP2, color = CellType)) +
  ggrastr::rasterise(geom_point(size = POINT_SIZE, alpha = PLOT_ALPHA), dpi = RASTER_DPI, dev = "cairo") +
  scale_color_manual(values = ct_colors) + coord_fixed() + labs(title = "CK cell types", x = "UMAP_1", y = "UMAP_2") + base_theme +
  guides(color = guide_legend(override.aes = list(size = 3)))
save_plot(reference_plot, "CK_celltypes", 7, 5)
plots <- list(); summaries <- list()
for (g in GENES) {
  if (!g %in% available_genes) {
    p <- ggplot() + annotate("text", x = 0, y = 0, label = "Absent from RNA assay") + labs(title = g) + theme_void()
    plots[[g]] <- p; next
  }
  d <- plot_data; d$Expression <- as.numeric(expr[g, ])
  positive <- d$Expression[d$Expression > 0]
  cap <- if (length(positive)) as.numeric(quantile(positive, .99, names = FALSE)) else 1
  d <- d[order(d$Expression), , drop = FALSE]
  # Zero-expression background is gray; positive cells are overlaid low to high.
  p <- ggplot(d, aes(UMAP1, UMAP2)) +
    ggrastr::rasterise(geom_point(color = "#D9DDDA", size = POINT_SIZE, alpha = PLOT_ALPHA), dpi = RASTER_DPI, dev = "cairo") +
    ggrastr::rasterise(geom_point(data = d[d$Expression > 0, , drop = FALSE], aes(color = Expression), size = POINT_SIZE, alpha = PLOT_ALPHA), dpi = RASTER_DPI, dev = "cairo") +
    scale_color_gradientn(colours = c("#ECE9A6", "#E7AF66", "#B75A6B", "#45316D"), limits = c(0, cap), oob = function(x, range) pmin(pmax(x, range[1]), range[2]), name = "RNA expression") +
    coord_fixed(xlim = range(plot_data$UMAP1), ylim = range(plot_data$UMAP2)) +
    labs(title = g, subtitle = "CK | positive-expression q99 color cap", x = "UMAP_1", y = "UMAP_2") + base_theme
  plots[[g]] <- p
  save_plot(p, paste0("CK_UMAP_", g))
  summaries[[g]] <- do.call(rbind, lapply(split(as.numeric(expr[g, ]), celltypes), function(v) data.frame(n_cells = length(v), mean_expression = mean(v), pct_expressing = 100 * mean(v > 0))))
  summaries[[g]]$celltype <- rownames(summaries[[g]]); summaries[[g]]$gene <- g; rownames(summaries[[g]]) <- NULL
}
# Keep an independent expression scale/legend for each gene.
combined <- wrap_plots(c(list(reference_plot + theme(legend.position = "none")), plots), ncol = 3)
save_plot(combined, "CK_in_situ_5genes_with_celltypes_2x3", 15, 9)
write_csv(do.call(rbind, summaries), "CK_gene_expression_summary.csv")

cells <- data.frame(id = plot_data$cell, x = unname(plot_data$UMAP1), y = unname(plot_data$UMAP2), celltype = plot_data$CellType, group = "CK", stringsAsFactors = FALSE)
expression <- setNames(lapply(available_genes, function(g) I(as.numeric(expr[g, ]))), available_genes)
atlas <- list(schema_version = 1L, source_label = basename(MERGED_RDS), created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), expression_label = expression_label, reduction = UMAP_REDUCTION, celltype_column = CELLTYPE_COL, group_selection = group_source, cells = cells, expression = expression)
jsonlite::write_json(atlas, file.path(OUT_DIR, "atlas.json"), dataframe = "rows", auto_unbox = TRUE, digits = 7, pretty = FALSE)
write_csv(data.frame(cell = plot_data$cell, group = "CK", celltype = celltypes), "selected_CK_cells.csv")
capture.output(sessionInfo(), file = file.path(OUT_DIR, "sessionInfo.txt"))
message("Done: ", OUT_DIR, "\nImport atlas.json into the website to view the real CK UMAP.")

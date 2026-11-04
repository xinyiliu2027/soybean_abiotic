# Source: Figure 3 document, extracted text lines 285-1142
# Standalone script requiring raw RNA counts for CK and three treatments.
# See README.md and docs/KNOWN_LIMITATIONS.md before running.

# Retain the original UMAP and add P1/P2/P3 labels without changing the analysis workflow.
# Use explicit dplyr/SeuratObject namespaces to resolve Assays() type conflicts.
# =============================================================================
# CK Endodermis reference -> High_Temp / Low_Temp / PEG response
# Input must contain raw RNA counts for CK and all three treatments.
# Three pseudotime intervals do not establish three discrete differentiation states; validate using markers and replicates.
# =============================================================================

# 0. User configuration.
MERGED_RDS <- file.path(Sys.getenv("SOYBEAN_PROJECT_ROOT", unset = "data"), "muscat_results_ck_ref_paper/CK_reference/CKref_merged_seurat_annotated.rds")
OUTPUT_ROOT <- file.path(Sys.getenv("SOYBEAN_PROJECT_ROOT", unset = "data"), "muscat_results_ck_ref_paper/CK_reference/Endodermis_CK_states_originalUMAP")
CELLTYPE_COL <- NULL                 # NULL detects the column automatically; alternatively specify cell_type.
TREATMENT_COL <- NULL                # NULL detects the column automatically; alternatively specify treatment.
SAMPLE_COL <- NULL                   # Prefer an explicit independent-sample column, such as sample_id.
BIOLOGICAL_REPLICATES_CONFIRMED <- FALSE
# Set TRUE only after confirming SAMPLE_COL identifies independent biological replicates; libraries/cell barcodes are not replicates.
# FALSE still exports state analyses, proportions and descriptive log2FC, without pseudoreplicated P values.
BATCH_COL <- NULL                    # Optional sample-level batch covariate.
BLOCK_COL <- NULL                    # Optional paired-individual column; use different SAMPLE_ID values for different treatments of the same individual.
ENDODERMIS_LABELS <- NULL            # Optionally specify c("Endodermis"); otherwise use pattern matching.
TREATMENT_LABELS <- list(
  CK = c("CK", "Control", "Mock"),
  High_Temp = c("HT", "HS", "Heat", "High_Temp", "HighTemp", "High_temperature", "\u9ad8\u6e29"),
  Low_Temp = c("LT", "CS", "Cold", "Low_Temp", "LowTemp", "Low_temperature", "\u4f4e\u6e29"),
  PEG = c("PEG", "Drought", "Osmotic", "\u5e72\u65f1")
)
# Default matching ignores case, whitespace, underscores and hyphens but otherwise requires exact matches; add CK1/CK2 explicitly if needed.

UMAP_REDUCTION <- "umap"             # Use the UMAP embedding already present in the input RDS.
UMAP_DIMS <- c(1L, 2L)
PERIOD_IDS <- c("P1", "P2", "P3")    # Order along CK pseudotime.
# To label the first interval P2, use c("P2", "P1", "P3").
UMAP_POINT_SIZE <- 0.35
UMAP_LABEL_SIZE <- 5
UMAP_LABEL_REPEL <- TRUE             # When available, ggrepel moves only text labels.
SAVE_UMAP_PNG <- TRUE
SEED <- 123L
N_HVG <- 2500L
N_PCS <- 20L
RESOLUTION <- 0.4
RESOLUTION_GRID <- c(0.1, 0.2, 0.4, 0.6, 0.8)
ROOT_CLUSTER <- NULL                 # After inspecting CK_cluster_markers.csv, specify a root cluster such as "0".
BIOLOGICAL_DIRECTION_CONFIRMED <- FALSE
LINEAGE_ID <- NULL                   # For branched trajectories, select and record the lineage with the most valid cells by default.
STATE_MODE <- "pseudotime_bins"      # Alternatively use curated_clusters.
CLUSTER_TO_STATE <- NULL             # For example, c("0" = "P1", "2" = "P2", "1" = "P3").
# Use curated_clusters only for marker-supported manual staging; unlisted clusters become Other_lineage.
# Default pseudotime_bins splits the main CK lineage at one-third and two-thirds; this does not establish exactly three biological states.
EARLY_MARKER_IDS <- character(0)     # Provide validated early-marker IDs for this species; do not invent ortholog assignments.
LATE_MARKER_IDS <- character(0)      # Provide validated mature-marker IDs for this species.
EXCLUDE_FROM_PCA <- character(0)     # Optionally exclude established technical/cell-cycle confounders for sensitivity analysis.
K_NEIGHBORS <- 30L
MIN_MAP_SCORE <- 0.60
MAP_DISTANCE_QUANTILE <- 0.99
MIN_CK_CELLS <- 90L
MIN_CELLS_PER_PB <- 20L
MIN_REPLICATES <- 2L                 # Minimum requirement; at least three replicates are preferable.
FDR_CUTOFF <- 0.05
LOG2FC_CUTOFF <- 0.5
RUN_GO <- TRUE
ORGDB_PACKAGE <- "org.G.max.eg.db"   # Soybean annotation package from the source; change it for other species.
ORGDB_KEYTYPE <- "GID"
GENE_ID_MAP_CSV <- NULL              # Optional CSV with gene and go_id columns; one-to-many mappings are supported.
MIN_GO_MAPPING_FRACTION <- 0.20
SAVE_SEURAT <- TRUE

# CRAN: install.packages(c("Seurat","dplyr","ggplot2","irlba","FNN","cluster"))
# Bioconductor: BiocManager::install(c("slingshot","edgeR","limma",
#   "ComplexHeatmap","circlize","clusterProfiler","AnnotationDbi"))
# The installed org.G.max.eg.db must match the annotation version used for the data.
required <- c("Seurat", "SeuratObject", "Matrix", "dplyr", "ggplot2",
              "irlba", "FNN", "cluster", "slingshot", "ComplexHeatmap", "circlize")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(dplyr); library(ggplot2)
})
if (BIOLOGICAL_REPLICATES_CONFIRMED &&
    !all(vapply(c("edgeR", "limma"), requireNamespace, logical(1), quietly = TRUE))) {
  stop("Formal pseudobulk testing requires edgeR and limma.")
}
set.seed(SEED)
OUT <- file.path(OUTPUT_ROOT, format(Sys.time(), "%Y%m%d_%H%M%S"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
for (s in c("CK", "mapping", "composition", "DE", "response", "GO", "UMAP")) {
  dir.create(file.path(OUT, s), showWarnings = FALSE)
}
notes <- character(0)
note <- function(x) {
  message(x)
  notes <<- c(notes, x)
  writeLines(notes, file.path(OUT, "analysis_notes.txt"), useBytes = TRUE)
}
csv <- function(x, path) write.csv(x, file.path(OUT, path), row.names = FALSE,
                                  na = "", fileEncoding = "UTF-8")
save_plot <- function(p, path, w = 9, h = 6) {
  ggsave(file.path(OUT, path), p, width = w, height = h, limitsize = FALSE)
}
heat <- function(m, path, title, zscore = FALSE, split = NULL) {
  if (nrow(m) < 2 || ncol(m) < 2) return(invisible(NULL))
  if (zscore) m <- t(scale(t(m)))
  m[!is.finite(m)] <- 0
  lim <- if (zscore) 2 else max(1, as.numeric(quantile(abs(m), .98)))
  grDevices::pdf(file.path(OUT, path), width = 10, height = 9)
  on.exit(grDevices::dev.off(), add = TRUE)
  ht <- ComplexHeatmap::Heatmap(m, name = if (zscore) "Row Z" else "log2FC",
    col = circlize::colorRamp2(c(-lim, 0, lim), c("#4575B4", "#FFFFBF", "#D73027")),
    cluster_rows = is.null(split), cluster_columns = FALSE, row_split = split,
    cluster_row_slices = FALSE, show_row_names = nrow(m) <= 40,
    column_title = title, use_raster = nrow(m) > 100)
  ComplexHeatmap::draw(ht)
}
pick_col <- function(md, specified, candidates, what, mandatory = TRUE) {
  if (!is.null(specified)) {
    if (!specified %in% colnames(md)) stop(what, "Column does not exist: ", specified)
    return(specified)
  }
  hits <- intersect(candidates, colnames(md))
  if (length(hits)) {
    note(paste0(what, "Using column: ", hits[1], if (length(hits) > 1)
      paste0("; other candidates: ", paste(hits[-1], collapse = ",")) else ""))
    return(hits[1])
  }
  if (mandatory) stop("Cannot find ", what, " column; configure explicitly. Metadata columns: ", paste(colnames(md), collapse = ","))
  NULL
}
norm_label <- function(x) tolower(gsub("[[:space:]_-]+", "", trimws(as.character(x))))
validate_rna_input <- function(obj) {
  if (!inherits(obj, "Seurat")) stop("Input must be a Seurat object; current class: ", paste(class(obj), collapse = ", "))
  # Read assay names directly to avoid Assays() returning an S4 container or a masked function.
  assay_names <- names(methods::slot(obj, "assays"))
  if (!"RNA" %in% assay_names) {
    stop("Input lacks RNA assay; available assays: ", paste(assay_names, collapse = ", "),
         ". Differential expression requires raw RNA counts.")
  }
  invisible(assay_names)
}

get_layer <- function(obj, layer) {
  if (inherits(obj[["RNA"]], "Assay5")) {
    SeuratObject::LayerData(obj, assay = "RNA", layer = layer)
  } else SeuratObject::GetAssayData(obj, assay = "RNA", slot = layer)
}
aggregate_counts <- function(counts, groups) {
  levels_g <- unique(as.character(groups))
  indicator <- Matrix::sparseMatrix(i = seq_along(groups), j = match(groups, levels_g),
    x = 1, dims = c(length(groups), length(levels_g)),
    dimnames = list(colnames(counts), levels_g))
  ans <- counts %*% indicator
  colnames(ans) <- levels_g
  ans
}

# Match original coordinates by cell name to avoid embedding/metadata order mismatches.
check_umap_coordinates <- function(reference, current, cells) {
  if (anyDuplicated(rownames(reference)) || anyDuplicated(rownames(current)) ||
      !all(cells %in% rownames(reference)) || !all(cells %in% rownames(current)))
    stop("UMAP cell names are duplicated or missing; coordinates cannot be matched reliably.")
  if (!isTRUE(all.equal(reference[cells, , drop = FALSE], current[cells, , drop = FALSE],
                        tolerance = 0, check.attributes = TRUE)))
    stop("UMAP coordinates changed; refusing to save.")
  invisible(TRUE)
}
assert_original_umap <- function(obj) {
  current <- SeuratObject::Embeddings(obj[[UMAP_REDUCTION]])
  check_umap_coordinates(original_umap, current, colnames(obj))
}
attach_original_xy <- function(dat) {
  if (!all(dat$cell %in% rownames(original_umap))) stop("Some cells lack original UMAP coordinates.")
  xy <- original_umap[dat$cell, UMAP_DIMS, drop = FALSE]
  dat$UMAP_1 <- xy[, 1]; dat$UMAP_2 <- xy[, 2]
  dat
}
period_label_positions <- function(dat, states) {
  # Anchor labels at the actual cell nearest the stage's two-dimensional median, avoiding empty regions.
  z <- dat[!is.na(dat$state) & dat$state %in% states &
             is.finite(dat$UMAP_1) & is.finite(dat$UMAP_2), , drop = FALSE]
  if (!nrow(z)) return(data.frame(group = character(), state = character(),
                                 UMAP_1 = numeric(), UMAP_2 = numeric(), n_cells = integer()))
  keys <- interaction(as.character(z$group), as.character(z$state), drop = TRUE)
  ans <- lapply(split(z, keys), function(d) {
    score <- (d$UMAP_1 - median(d$UMAP_1))^2 + (d$UMAP_2 - median(d$UMAP_2))^2
    a <- d[which.min(score), c("group", "state", "UMAP_1", "UMAP_2"), drop = FALSE]
    a$n_cells <- nrow(d); a
  })
  out <- do.call(rbind, ans); rownames(out) <- NULL
  out
}
original_umap_plot <- function(dat, title, color_col = "state", kind = "state",
                               facets = FALSE, full_context = FALSE, label_periods = TRUE) {
  if (!all(c("UMAP_1", "UMAP_2") %in% names(dat))) dat <- attach_original_xy(dat)
  dat$group <- factor(as.character(dat$group), levels = group_levels)
  # Use the same viewing limits and aspect ratio across groups without stretching or re-embedding individual groups.
  bounds <- if (full_context) original_umap[, UMAP_DIMS, drop = FALSE] else
    original_umap[colnames(endo), UMAP_DIMS, drop = FALSE]
  padded_range <- function(v) {
    r <- range(v); padding <- max(diff(r) * .04, 1e-6); r + c(-padding, padding)
  }
  p <- ggplot()
  if (full_context) {
    bg <- original_umap_context[original_umap_context$group %in% as.character(unique(dat$group)), , drop = FALSE]
    bg$group <- factor(bg$group, levels = group_levels)
    p <- p + geom_point(data = bg, aes(UMAP_1, UMAP_2), color = "grey87", size = .18, alpha = .45)
  }
  p <- p + geom_point(data = dat, aes(UMAP_1, UMAP_2, color = .data[[color_col]]),
                       size = UMAP_POINT_SIZE, alpha = .8)
  if (kind == "state") {
    colors <- c(setNames(c("#3B528B", "#21918C", "#FDE725"), state_names),
                Other_lineage = "#AAAAAA", Unassigned = "#D9D9D9")
    p <- p + scale_color_manual(values = colors, limits = names(colors), drop = FALSE,
                                 labels = c(unname(period_display), "Other lineage", "Unassigned"), name = "Period")
  } else if (kind == "pseudotime") {
    p <- p + scale_color_gradientn(colors = c("#440154", "#31688E", "#35B779", "#FDE725"),
      limits = range(pt[is.finite(pt)]), na.value = "grey80", name = "CK reference\npseudotime")
  } else if (kind == "score") {
    p <- p + scale_color_gradient(low = "grey85", high = "#B2182B", name = color_col)
  } else p <- p + labs(color = color_col)
  if (label_periods && "state" %in% names(dat)) {
    labs_df <- period_label_positions(dat, state_names)
    if (nrow(labs_df)) {
      # Use one label per stage in the combined plot and separate labels within each treatment panel.
      if (!facets) {
        temp <- dat; temp$group <- "all"
        labs_df <- period_label_positions(temp, state_names)
      }
      labs_df$label <- unname(period_display[as.character(labs_df$state)])
      if (facets) labs_df$group <- factor(labs_df$group, levels = group_levels)
      label_layer <- if (UMAP_LABEL_REPEL && requireNamespace("ggrepel", quietly = TRUE)) {
        ggrepel::geom_label_repel(data = labs_df, aes(UMAP_1, UMAP_2, label = label),
          inherit.aes = FALSE, seed = SEED, size = UMAP_LABEL_SIZE, fontface = "bold",
          fill = "white", color = "black", min.segment.length = 0, max.overlaps = Inf,
          box.padding = .4, show.legend = FALSE)
      } else geom_label(data = labs_df, aes(UMAP_1, UMAP_2, label = label),
        inherit.aes = FALSE, size = UMAP_LABEL_SIZE, fontface = "bold",
        fill = "white", color = "black", show.legend = FALSE)
      p <- p + label_layer
    }
  }
  if (facets) p <- p + facet_wrap(~group, nrow = 1, scales = "fixed")
  p + coord_fixed(ratio = 1, xlim = padded_range(bounds[, 1]), ylim = padded_range(bounds[, 2]), expand = FALSE) +
    labs(title = title, x = colnames(original_umap)[UMAP_DIMS[1]], y = colnames(original_umap)[UMAP_DIMS[2]]) +
    theme_classic(base_size = 12) + theme(panel.border = element_rect(color = "grey40", fill = NA),
      plot.title = element_text(face = "bold"), strip.background = element_rect(fill = "grey95"))
}
save_original_umap <- function(p, filename, w = 8, h = 7) {
  save_plot(p, paste0("UMAP/", filename, ".pdf"), w, h)
  if (SAVE_UMAP_PNG) ggsave(file.path(OUT, "UMAP", paste0(filename, ".png")),
                           plot = p, width = w, height = h, dpi = 300, limitsize = FALSE)
}

# 1. Read and validate inputs.
if (!file.exists(MERGED_RDS)) stop("File does not exist: ", MERGED_RDS)
raw <- readRDS(MERGED_RDS)
validate_rna_input(raw)
existing_reductions <- names(methods::slot(raw, "reductions"))
if (!UMAP_REDUCTION %in% existing_reductions) stop(
  "Specified original UMAP missing: ", UMAP_REDUCTION, "; available reductions: ",
  paste(existing_reductions, collapse = ", "), ". Set UMAP_REDUCTION appropriately; this script does not recompute UMAP.")
original_umap <- SeuratObject::Embeddings(raw[[UMAP_REDUCTION]])
if (length(UMAP_DIMS) != 2L || anyNA(UMAP_DIMS) || any(UMAP_DIMS != as.integer(UMAP_DIMS)) ||
    any(UMAP_DIMS < 1 | UMAP_DIMS > ncol(original_umap)) || anyDuplicated(UMAP_DIMS))
  stop("UMAP_DIMS must select two distinct valid dimensions of the original UMAP.")
if (is.null(rownames(original_umap)) || anyDuplicated(rownames(original_umap)) ||
    any(!is.finite(original_umap))) stop("Original UMAP contains invalid coordinates or cell names.")
if (length(PERIOD_IDS) != 3 || anyNA(PERIOD_IDS) || anyDuplicated(PERIOD_IDS) ||
    any(!grepl("^P[0-9]+$", PERIOD_IDS))) stop("PERIOD_IDS must contain three distinct labels consisting of P followed by digits.")
saveRDS(original_umap, file.path(OUT, "UMAP/original_UMAP_coordinates_snapshot.rds"))
note(paste("Preserving original UMAP: ", UMAP_REDUCTION, "; no UMAP recomputation, with coordinates verified for every cell."))
md <- raw[[]]
csv(data.frame(column = colnames(md)), "metadata_columns.csv")
ct_col <- pick_col(md, CELLTYPE_COL,
  c("cell_type", "merged_cell_type", "new_cell_type", "CellType", "celltype"), "cell type")
tr_col <- pick_col(md, TREATMENT_COL, c("treatment", "group", "condition", "Treatment"), "treatment")
sample_col <- pick_col(md, SAMPLE_COL,
  c("sample_id", "sample", "Sample", "sample_name", "biological_sample", "orig.ident"), "sample", FALSE)
csv(as.data.frame(table(as.character(md[[ct_col]]))), "celltype_labels.csv")
csv(as.data.frame(table(as.character(md[[tr_col]]))), "treatment_labels.csv")
ct_values <- if (is.null(ENDODERMIS_LABELS)) {
  unique(as.character(md[[ct_col]])[grepl("endodermis|\u5185\u76ae\u5c42", md[[ct_col]], ignore.case = TRUE)])
} else ENDODERMIS_LABELS
ct_values <- ct_values[!is.na(ct_values)]
note(paste("Endodermis labels: ", paste(ct_values, collapse = ", ")))
if (!length(ct_values)) stop("No Endodermis labels matched.")
tr_map <- do.call(rbind, lapply(names(TREATMENT_LABELS), function(g)
  data.frame(key = norm_label(TREATMENT_LABELS[[g]]), group = g)))
if (any(vapply(split(tr_map$group, tr_map$key), function(x) length(unique(x)) > 1, logical(1))))
  stop("TREATMENT_LABELS contains labels shared across groups.")
group_all <- tr_map$group[match(norm_label(md[[tr_col]]), tr_map$key)]
original_umap_context <- data.frame(cell = rownames(original_umap),
  group = group_all[match(rownames(original_umap), rownames(md))],
  UMAP_1 = original_umap[, UMAP_DIMS[1]], UMAP_2 = original_umap[, UMAP_DIMS[2]])
label_audit <- unique(data.frame(original_label = as.character(md[[tr_col]]), group = group_all))
csv(label_audit, "treatment_label_mapping.csv")
keep <- as.character(md[[ct_col]]) %in% ct_values & !is.na(group_all)
if (!any(keep)) stop("No Endodermis cells matched; inspect the labels and adjust configuration.")
endo <- subset(raw, cells = rownames(md)[keep])
endo$analysis_group <- group_all[match(colnames(endo), rownames(md))]
assert_original_umap(endo)
rm(raw); invisible(gc())
DefaultAssay(endo) <- "RNA"
if (inherits(endo[["RNA"]], "Assay5")) endo <- SeuratObject::JoinLayers(endo, assay = "RNA")
counts <- get_layer(endo, "counts")
counts <- counts[, colnames(endo), drop = FALSE]
if (!nrow(counts) || any(Matrix::colSums(counts) <= 0)) stop("RNA counts are empty or contain zero-library cells.")
vals <- if (inherits(counts, "sparseMatrix")) counts@x else as.numeric(counts)
if (any(!is.finite(vals)) || any(vals < 0) || any(abs(vals - round(vals)) > 1e-6))
  stop("RNA counts must be raw nonnegative integers; integrated or log-transformed values cannot be used for pseudobulk.")
rm(vals)
group_levels <- c("CK", "High_Temp", "Low_Temp", "PEG")
if (!all(group_levels %in% endo$analysis_group)) stop("Endodermis lacks these groups: ",
  paste(setdiff(group_levels, endo$analysis_group), collapse = ", "), "; check inputs and label mappings.")
if (is.null(sample_col)) {
  endo$analysis_sample <- paste0("pooled_", endo$analysis_group)
  if (BIOLOGICAL_REPLICATES_CONFIRMED) stop("No independent-sample column found.")
  note("No sample column: using within-group pooled IDs for descriptive comparisons only.")
} else {
  endo$analysis_sample <- as.character(endo[[]][[sample_col]])
  if (anyNA(endo$analysis_sample) || any(!nzchar(endo$analysis_sample))) stop("Sample IDs cannot be missing.")
}
sample_info <- unique(endo[[]][, c("analysis_sample", "analysis_group"), drop = FALSE])
if (anyDuplicated(sample_info$analysis_sample)) stop("A sample ID spans multiple treatments. Use unique sample IDs and place paired-individual IDs in BLOCK_COL.")
for (co in unique(c(BATCH_COL, BLOCK_COL))) {
  if (!co %in% colnames(endo[[]])) stop("Covariate column missing: ", co)
  sample_values <- split(as.character(endo[[]][[co]]), endo$analysis_sample)
  if (any(vapply(sample_values, function(x) anyNA(x) || any(!nzchar(x)) || length(unique(x)) != 1,
                 logical(1)))) stop("Covariate is missing or varies within a sample: ", co)
}
csv(as.data.frame(table(Sample = endo$analysis_sample, Group = endo$analysis_group)), "sample_cell_counts.csv")
if (!BIOLOGICAL_REPLICATES_CONFIRMED) note("Independent biological replicates not confirmed: report descriptive treatment effects only, skipping inferential tests and treatment DEG GO.")
endo <- NormalizeData(endo, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
expr <- get_layer(endo, "data")[, colnames(endo), drop = FALSE]
ck_cells <- colnames(endo)[endo$analysis_group == "CK"]
if (length(ck_cells) < MIN_CK_CELLS) stop("Insufficient CK Endodermis cells; required: ", MIN_CK_CELLS, "。")
ck <- subset(endo, cells = ck_cells)

# 2. Build a reference using CK alone.
ck <- FindVariableFeatures(ck, selection.method = "vst", nfeatures = N_HVG, verbose = FALSE)
hvg <- setdiff(VariableFeatures(ck), EXCLUDE_FROM_PCA)
x <- t(expr[hvg, ck_cells, drop = FALSE])
mu <- Matrix::colMeans(x)
sdv <- sqrt(pmax(0, (Matrix::colSums(x^2) - nrow(x) * mu^2) / (nrow(x) - 1)))
ok <- is.finite(sdv) & sdv > 1e-8
x <- x[, ok, drop = FALSE]; mu <- mu[ok]; sdv <- sdv[ok]; hvg <- hvg[ok]
npcs <- min(N_PCS, nrow(x) - 2L, ncol(x) - 2L)
if (npcs < 3) stop("Too few cells/genes for PCA.")
pca_fit <- irlba::prcomp_irlba(x, n = npcs, center = mu, scale. = sdv)
# Use the same projection formula for CK and query to avoid spurious mismatches from approximate SVD U*D residuals.
pca_weights <- sweep(pca_fit$rotation, 1, sdv, "/")
pc_ck <- sweep(as.matrix(x %*% pca_weights), 2, as.numeric(mu %*% pca_weights), "-")
rownames(pc_ck) <- ck_cells; colnames(pc_ck) <- paste0("CKPC_", seq_len(npcs))
pca_fit$x <- pc_ck
ck[["ckpca"]] <- CreateDimReducObject(embeddings = pc_ck, key = "CKPC_", assay = "RNA")
ck <- FindNeighbors(ck, reduction = "ckpca", dims = seq_len(npcs),
                    k.param = min(20, ncol(ck) - 1L), verbose = FALSE)
res_grid <- sort(unique(c(RESOLUTION_GRID, RESOLUTION)))
res_labels <- list(); res_summary <- list()
sil_idx <- sample(seq_len(nrow(pc_ck)), min(2000L, nrow(pc_ck)))
sil_dist <- dist(pc_ck[sil_idx, , drop = FALSE])
for (rr in res_grid) {
  ck <- FindClusters(ck, resolution = rr, random.seed = SEED, verbose = FALSE)
  labs <- as.character(Idents(ck)); res_labels[[as.character(rr)]] <- labs
  sizes <- table(labs); labs_sub <- as.integer(factor(labs[sil_idx]))
  sil <- if (length(unique(labs_sub)) > 1 && length(unique(labs_sub)) < length(labs_sub))
    mean(cluster::silhouette(labs_sub, sil_dist)[, "sil_width"]) else NA_real_
  res_summary[[as.character(rr)]] <- data.frame(resolution = rr, n_clusters = length(sizes),
    min_cluster_size = min(sizes), mean_silhouette = sil)
}
csv(dplyr::bind_rows(res_summary), "CK/CK_resolution_diagnostics.csv")
csv(data.frame(cell = ck_cells, as.data.frame(res_labels, check.names = FALSE)), "CK/CK_resolution_labels.csv")
ck$CK_cluster <- res_labels[[as.character(RESOLUTION)]]
Idents(ck) <- "CK_cluster"
csv(as.data.frame(table(Cluster = ck$CK_cluster, Sample = ck$analysis_sample)), "CK/CK_cluster_by_sample.csv")
note(paste("CK resolution =", RESOLUTION, "Detected ", length(unique(ck$CK_cluster)), " transcriptional subclusters; cluster count does not directly establish differentiation-stage count."))
ck_markers <- FindAllMarkers(ck, only.pos = TRUE, min.pct = .10, logfc.threshold = .25,
                             test.use = "wilcox", verbose = FALSE)
csv(ck_markers, "CK/CK_cluster_markers_exploratory.csv")
note("CK cluster/state marker P values are exploratory because clustering and marker selection use the same data; they are not independently validated differentiation evidence.")
assert_original_umap(ck)
ck_umap <- attach_original_xy(data.frame(cell = colnames(ck), group = "CK",
  cluster = ck$CK_cluster, sample = ck$analysis_sample))
save_original_umap(original_umap_plot(ck_umap, "CK clusters - original UMAP",
  color_col = "cluster", kind = "category", label_periods = FALSE), "CK_clusters_originalUMAP")
save_original_umap(original_umap_plot(ck_umap, "CK samples - original UMAP",
  color_col = "sample", kind = "category", label_periods = FALSE), "CK_samples_originalUMAP")
if (length(unique(ck$CK_cluster)) < 2) {
  saveRDS(ck, file.path(OUT, "CK/CK_reference.rds"))
  stop("Only one cluster at the current resolution. CK diagnostics saved; three stages will not be forced. Check markers, cell counts and resolution.")
}

# 3. Infer the CK trajectory and define three states.
if (!is.null(ROOT_CLUSTER) && !ROOT_CLUSTER %in% ck$CK_cluster) stop("ROOT_CLUSTER is not present in CK_cluster.")
if (BIOLOGICAL_DIRECTION_CONFIRMED && is.null(ROOT_CLUSTER)) stop("A marker-validated ROOT_CLUSTER is required before confirming developmental direction.")
args_ss <- list(data = pc_ck, clusterLabels = ck$CK_cluster, approx_points = 150)
if (!is.null(ROOT_CLUSTER)) args_ss$start.clus <- ROOT_CLUSTER
ss <- do.call(slingshot::slingshot, args_ss)
ptmat <- slingshot::slingPseudotime(ss)
rownames(ptmat) <- ck_cells
paths <- slingshot::slingLineages(ss)
csv(data.frame(cell = ck_cells, ptmat, check.names = FALSE), "CK/CK_all_lineage_pseudotimes.csv")
lin_summary <- data.frame(lineage = seq_len(ncol(ptmat)),
  n_valid = colSums(is.finite(ptmat)), cluster_path = vapply(paths, paste, collapse = " -> ", character(1)))
csv(lin_summary, "CK/CK_lineages.csv")
lid <- if (is.null(LINEAGE_ID)) which.max(colSums(is.finite(ptmat))) else as.integer(LINEAGE_ID)
if (length(lid) != 1 || is.na(lid) || lid < 1 || lid > ncol(ptmat)) stop("LINEAGE_ID is out of range.")
note(paste("Selected lineage ", lid, "; total lineages ", ncol(ptmat), ". Pseudotimes are not averaged across branches."))
if (ncol(ptmat) > 1) note("Branches detected: default analysis uses the lineage with the most valid cells. See CK_all_lineage_pseudotimes.csv and rerun with LINEAGE_ID if appropriate.")
pt <- ptmat[, lid]; valid_pt <- is.finite(pt)
if (sum(valid_pt) < MIN_CK_CELLS) stop("Too few valid cells in the selected lineage.")
if (diff(range(pt[valid_pt])) <= 0) stop("Pseudotime has no variation.")
ck$CK_pseudotime <- pt
state_names <- PERIOD_IDS
period_display <- setNames(if (BIOLOGICAL_DIRECTION_CONFIRMED)
  paste0(state_names, " (", c("Early", "Middle", "Late"), ")") else state_names, state_names)
csv(data.frame(period = state_names, label = unname(period_display),
  reference_order = 1:3, direction_confirmed = BIOLOGICAL_DIRECTION_CONFIRMED), "UMAP/period_label_key.csv")
if (!BIOLOGICAL_DIRECTION_CONFIRMED) note(paste(paste(state_names, collapse = "/"), "denotes reference intervals along the algorithmic direction. Direction is unvalidated; do not interpret as early/middle/late differentiation."))
ck_state <- rep("Other_lineage", ncol(ck)); cuts <- NULL
if (STATE_MODE == "pseudotime_bins") {
  cuts <- as.numeric(quantile(pt[valid_pt], c(1/3, 2/3), names = FALSE))
  if (length(unique(cuts)) < 2) stop("Duplicated pseudotime quantiles prevent a reliable three-way split.")
  ck_state[valid_pt] <- as.character(cut(pt[valid_pt], c(-Inf, cuts, Inf), labels = state_names))
  csv(data.frame(boundary = c("one_third", "two_thirds"), pseudotime = cuts), "CK/CK_fixed_state_boundaries.csv")
  note("CK pseudotime quantiles define three intervals, so approximately one-third of CK main-lineage cells in each interval is by construction, not statistical evidence for three states.")
} else if (STATE_MODE == "curated_clusters") {
  if (is.null(CLUSTER_TO_STATE) || is.null(names(CLUSTER_TO_STATE)) || anyDuplicated(names(CLUSTER_TO_STATE)))
    stop("curated_clusters requires a uniquely named cluster-to-PERIOD_IDS vector.")
  if (!all(names(CLUSTER_TO_STATE) %in% ck$CK_cluster)) stop("CLUSTER_TO_STATE contains nonexistent CK clusters.")
  if (all(CLUSTER_TO_STATE %in% c("S1", "S2", "S3"))) {
    CLUSTER_TO_STATE <- setNames(state_names[match(CLUSTER_TO_STATE, c("S1", "S2", "S3"))], names(CLUSTER_TO_STATE))
  }
  if (!all(CLUSTER_TO_STATE %in% state_names)) stop("CLUSTER_TO_STATE values must belong to PERIOD_IDS.")
  z <- unname(CLUSTER_TO_STATE[ck$CK_cluster])
  ck_state[!is.na(z)] <- z[!is.na(z)]
  note("Using manually curated, marker-supported cluster-to-state annotations.")
} else stop("Unknown STATE_MODE.")
if (!all(state_names %in% ck_state)) stop("Not all three states were obtained; inspect the staging definition.")
ck$CK_state <- ck_state
csv(data.frame(cell = ck_cells, cluster = ck$CK_cluster, sample = ck$analysis_sample,
               state = ck_state, pseudotime = pt), "CK/CK_cell_states.csv")
ck_umap$state <- ck_state
ck_umap$pseudotime <- pt
save_original_umap(original_umap_plot(ck_umap, "CK periods - original UMAP"), "CK_periods_originalUMAP")
save_original_umap(original_umap_plot(ck_umap, "CK pseudotime and periods - original UMAP",
  color_col = "pseudotime", kind = "pseudotime"), "CK_pseudotime_P_labels_originalUMAP")
save_original_umap(original_umap_plot(ck_umap, "CK Endodermis periods in the original full UMAP",
  full_context = TRUE), "CK_periods_full_originalUMAP")
csv(period_label_positions(ck_umap, state_names), "UMAP/CK_period_label_anchors.csv")
# Fit and plot Slingshot curves in PCA space; do not overlay PCA curves directly on UMAP.
grDevices::pdf(file.path(OUT, "CK/CK_trajectory_PCA.pdf"), width = 8, height = 7)
pt_colors <- rep("grey80", length(pt))
pt_colors[valid_pt] <- hcl.colors(100, "Viridis")[as.integer(cut(pt[valid_pt], 100, include.lowest = TRUE))]
plot(pc_ck[, 1:2], col = pt_colors, pch = 16, cex = .45, main = "CK Slingshot in PCA space")
for (cv in slingshot::slingCurves(ss)) lines(cv$s[cv$ord, 1:2, drop = FALSE], lwd = 2)
grDevices::dev.off()
for (nm in c("early", "late")) {
  ids <- intersect(if (nm == "early") EARLY_MARKER_IDS else LATE_MARKER_IDS, rownames(expr))
  if (length(ids)) {
    ck[[paste0(nm, "_marker_mean")]] <- as.numeric(Matrix::colMeans(expr[ids, ck_cells, drop = FALSE]))
    ck_umap$marker_mean <- ck[[]][[paste0(nm, "_marker_mean")]]
    save_original_umap(original_umap_plot(ck_umap, paste("CK", nm, "marker score - original UMAP"),
      color_col = "marker_mean", kind = "score"), paste0("CK_", nm, "_marker_score_originalUMAP"))
    csv(data.frame(cell = ck_cells, state = ck_state, pt = pt,
      score = ck[[]][[paste0(nm, "_marker_mean")]]), paste0("CK/CK_", nm, "_marker_scores.csv"))
  }
}
ck_stage <- subset(ck, cells = ck_cells[ck_state %in% state_names])
Idents(ck_stage) <- "CK_state"
stage_markers <- FindAllMarkers(ck_stage, only.pos = TRUE, min.pct = .10,
                                 logfc.threshold = .25, verbose = FALSE)
csv(stage_markers, "CK/CK_state_markers_exploratory.csv")
if (nrow(stage_markers)) {
  fc_col <- intersect(c("avg_log2FC", "avg_logFC"), colnames(stage_markers))[1]
  top <- stage_markers %>% dplyr::group_by(cluster) %>% dplyr::arrange(dplyr::desc(.data[[fc_col]]), .by_group = TRUE) %>% dplyr::slice_head(n = 15) %>% dplyr::ungroup()
  genes <- unique(top$gene)
  avg <- sapply(state_names, function(s) Matrix::rowMeans(expr[genes, ck_cells[ck_state == s], drop = FALSE]))
  rownames(avg) <- genes
  heat(avg, "CK/CK_state_marker_heatmap.pdf", "CK state marker means", zscore = TRUE)
}
# Show CK variable genes along the selected trajectory; compute correlations in blocks rather than densifying all genes by all cells.
ord_cells <- ck_cells[valid_pt][order(pt[valid_pt])]
pt_ordered <- sort(pt[valid_pt]); rho <- setNames(rep(NA_real_, length(hvg)), hvg)
for (idx in split(seq_along(hvg), ceiling(seq_along(hvg) / 100))) {
  block <- as.matrix(expr[hvg[idx], ord_cells, drop = FALSE])
  rho[idx] <- apply(block, 1, function(v) if (sd(v) > 0) cor(v, pt_ordered, method = "spearman") else NA_real_)
}
csv(data.frame(gene = names(rho), spearman_rho = rho), "CK/CK_pseudotime_gene_correlations_descriptive.csv")
g_pt <- head(names(sort(abs(rho), decreasing = TRUE)), 300)
nbins <- min(60L, floor(length(ord_cells) / 5))
bin <- cut(seq_along(ord_cells), breaks = nbins, labels = FALSE)
binned <- sapply(seq_len(nbins), function(b) Matrix::rowMeans(expr[g_pt, ord_cells[bin == b], drop = FALSE]))
rownames(binned) <- g_pt; colnames(binned) <- paste0("B", seq_len(nbins))
heat(binned, "CK/CK_pseudotime_heatmap.pdf", "CK ordered pseudotime bins", zscore = TRUE)

# 4. Map stress-treated cells to the CK reference.
# Fix CK HVGs, means, standard deviations and PCA loadings; query cells do not contribute to reference fitting.
rotation <- pca_fit$rotation
weighted_rotation <- sweep(rotation, 1, sdv, "/")
offset <- as.numeric(mu %*% weighted_rotation)
pc_all <- as.matrix(t(expr[hvg, , drop = FALSE]) %*% weighted_rotation)
pc_all <- sweep(pc_all, 2, offset, "-")
rownames(pc_all) <- colnames(endo); colnames(pc_all) <- colnames(pc_ck)
if (max(abs(pc_all[ck_cells, , drop = FALSE] - pc_ck)) > 1e-4)
  stop("CK PCA projection consistency check failed.")
k <- min(K_NEIGHBORS, nrow(pc_ck) - 1L)
ref_nn <- FNN::get.knn(pc_ck, k = k)    # Leave one out to prevent each cell from being its own nearest neighbor.
map_nn <- FNN::get.knnx(pc_ck, pc_all, k = k)
map_nn$nn.index[match(ck_cells, rownames(pc_all)), ] <- ref_nn$nn.index
map_nn$nn.dist[match(ck_cells, rownames(pc_all)), ] <- ref_nn$nn.dist
ref_radius <- rowMeans(ref_nn$nn.dist)
# Calibrate distances using local neighborhood density to avoid rejecting sparse CK states with a global threshold.
relative_distance <- function(nn) {
  denom <- rowMeans(matrix(ref_radius[nn$nn.index], nrow = nrow(nn$nn.index)))
  rowMeans(nn$nn.dist) / pmax(denom, 1e-8)
}
ref_relative <- relative_distance(ref_nn)
dist_cut <- as.numeric(quantile(ref_relative, MAP_DISTANCE_QUANTILE))
query_relative <- relative_distance(map_nn)
label_levels <- c(state_names, "Other_lineage")
neighbor_labels <- matrix(ck_state[map_nn$nn.index], nrow = nrow(pc_all))
weights <- 1 / pmax(map_nn$nn.dist, 1e-8)
weights <- weights / rowSums(weights)
votes <- sapply(label_levels, function(s) rowSums(weights * (neighbor_labels == s)))
prediction <- label_levels[max.col(votes, ties.method = "first")]
confidence <- apply(votes, 1, max)
accepted <- confidence >= MIN_MAP_SCORE & query_relative <= dist_cut
mapped_state <- ifelse(accepted, prediction, "Unassigned")
# Retain CK reference assignments; report leave-one-out agreement as a diagnostic, not independent accuracy.
ck_idx <- match(ck_cells, colnames(endo))
csv(as.data.frame(table(Reference = ck_state, LOO = mapped_state[ck_idx])), "mapping/CK_leave_one_out_confusion.csv")
note(sprintf("CK leave-one-out neighbor agreement with reference labels (including rejections): %.3f; this is not independent validation accuracy.",
             mean(mapped_state[ck_idx] == ck_state)))
mapped_state[ck_idx] <- ck_state
# Transfer reference pseudotime for display only; it is not developmental time independently inferred in treated cells.
neighbor_pt <- matrix(pt[map_nn$nn.index], nrow = nrow(pc_all))
finite_neighbor <- is.finite(neighbor_pt); neighbor_pt[!finite_neighbor] <- 0
transferred_pt <- rowSums(weights * neighbor_pt) / rowSums(weights * finite_neighbor)
transferred_pt[!is.finite(transferred_pt) | !accepted | !mapped_state %in% state_names] <- NA_real_
transferred_pt[ck_idx] <- pt
endo$CK_reference_state <- mapped_state
endo$CK_transferred_pseudotime <- transferred_pt
endo$CK_mapping_score <- confidence
endo$CK_relative_distance <- query_relative
map_table <- data.frame(cell = colnames(endo), sample = endo$analysis_sample,
  group = endo$analysis_group, state = mapped_state, predicted_state = prediction,
  mapping_score = confidence, relative_distance = query_relative,
  within_CK_distance = query_relative <= dist_cut, transferred_pseudotime = transferred_pt,
  label_source = ifelse(endo$analysis_group == "CK", "reference_definition", "CK_neighbor_mapping"))
csv(map_table, "mapping/all_cells_CK_reference_states.csv")
# Use the treatment cells' existing UMAP coordinates; do not treat projected PCA coordinates as UMAP.
assert_original_umap(endo)
assert_original_umap(ck)
all_umap <- attach_original_xy(map_table)
csv(all_umap, "UMAP/all_cells_originalUMAP_periods.csv")
csv(period_label_positions(all_umap, state_names), "UMAP/all_groups_period_label_anchors.csv")
save_original_umap(original_umap_plot(all_umap, "CK and stress periods - original UMAP",
  facets = TRUE), "all_groups_periods_originalUMAP", 18, 5)
save_original_umap(original_umap_plot(all_umap, "CK-reference pseudotime and periods - original UMAP",
  color_col = "transferred_pseudotime", kind = "pseudotime", facets = TRUE),
  "all_groups_pseudotime_P_labels_originalUMAP", 18, 5)
save_original_umap(original_umap_plot(all_umap, "Endodermis periods in original full UMAP",
  facets = TRUE, full_context = TRUE), "all_groups_periods_full_originalUMAP", 18, 5)
save_original_umap(original_umap_plot(all_umap, "Endodermis periods - original UMAP, all groups"),
  "merged_periods_originalUMAP")
for (g in group_levels) {
  d <- all_umap[all_umap$group == g, , drop = FALSE]
  save_original_umap(original_umap_plot(d, paste(g, "periods - original UMAP")),
    paste0(g, "_periods_originalUMAP"))
  save_original_umap(original_umap_plot(d, paste(g, "pseudotime and periods - original UMAP"),
    color_col = "transferred_pseudotime", kind = "pseudotime"),
    paste0(g, "_pseudotime_P_labels_originalUMAP"))
}
note("All UMAP panels retain original input coordinates; P labels are overlays only. Trajectories are shown separately in fitted PCA space, without invented UMAP curves.")
csv(data.frame(score_cutoff = MIN_MAP_SCORE, relative_distance_cutoff = dist_cut,
  distance_quantile = MAP_DISTANCE_QUANTILE, k = k), "mapping/mapping_thresholds.csv")
map_df <- cbind(map_table, PC1 = pc_all[, 1], PC2 = pc_all[, 2])
map_df$group <- factor(map_df$group, group_levels)
save_plot(ggplot(map_df, aes(PC1, PC2, color = state)) + geom_point(size = .25, alpha = .6) +
  facet_wrap(~group, nrow = 1) + theme_bw() + labs(title = "Projection onto the same CK PCA reference"),
  "mapping/all_groups_CK_PCA.pdf", 16, 5)
save_plot(ggplot(map_df, aes(mapping_score, relative_distance, color = group)) +
  geom_point(size = .4, alpha = .3) + geom_vline(xintercept = MIN_MAP_SCORE, linetype = 2) +
  geom_hline(yintercept = dist_cut, linetype = 2) + theme_bw(), "mapping/mapping_QC.pdf")
note("Unassigned denotes low-confidence or CK-distant cells, potentially including new stress states or batch effects; cells are not forced into three states.")

# 5. State proportions.
all_states <- c(state_names, "Other_lineage", "Unassigned")
composition <- expand.grid(sample = sample_info$analysis_sample, state = all_states, stringsAsFactors = FALSE)
composition$group <- sample_info$analysis_group[match(composition$sample, sample_info$analysis_sample)]
observed <- map_table %>% dplyr::count(sample, state, name = "n_cells")
composition <- dplyr::left_join(composition, observed, by = c("sample", "state"))
composition$n_cells[is.na(composition$n_cells)] <- 0L
composition <- composition %>% dplyr::group_by(sample) %>% dplyr::mutate(total_endodermis = sum(n_cells),
  proportion = n_cells / total_endodermis) %>% dplyr::ungroup()
csv(composition, "composition/state_proportions_by_sample.csv")
composition$state <- factor(composition$state, all_states)
composition$group <- factor(composition$group, group_levels)
save_plot(ggplot(composition, aes(sample, proportion, fill = state)) + geom_col() +
  facet_grid(~group, scales = "free_x", space = "free_x") + theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(y = "Fraction of all sampled Endodermis cells", x = NULL),
  "composition/state_proportions_by_sample.pdf", 12, 6)
comp_summary <- composition %>% dplyr::group_by(group, state) %>% dplyr::summarise(
  mean_proportion = mean(proportion), sd_proportion = sd(proportion), n_samples = dplyr::n(), .groups = "drop")
ck_prop <- comp_summary %>% dplyr::filter(group == "CK") %>% dplyr::select(state, CK_mean = mean_proportion)
comp_summary <- dplyr::left_join(comp_summary, ck_prop, by = "state") %>% dplyr::mutate(delta_vs_CK = mean_proportion - CK_mean)
csv(comp_summary, "composition/state_proportion_changes_vs_CK.csv")
note("Proportion denominators include all Endodermis cells per sample, including Unassigned/other branches. Changes describe sampled composition, not differentiation rates or absolute cell numbers.")

# 6. Sample-by-state pseudobulk.
pb_md <- map_table[map_table$state %in% state_names, c("cell", "sample", "group", "state")]
pb_key <- interaction(factor(pb_md$sample), factor(pb_md$state), drop = TRUE, lex.order = TRUE)
pb_md$pb_id <- paste0("PB", as.integer(pb_key))
pb_counts <- aggregate_counts(counts[, pb_md$cell, drop = FALSE], pb_md$pb_id)
pb_info <- pb_md %>% dplyr::group_by(pb_id, sample, group, state) %>% dplyr::summarise(n_cells = dplyr::n(), .groups = "drop")
pb_info <- as.data.frame(pb_info[match(colnames(pb_counts), pb_info$pb_id), ])
pb_info$eligible <- pb_info$n_cells >= MIN_CELLS_PER_PB
for (co in unique(c(BATCH_COL, BLOCK_COL))) {
  vals_co <- endo[[]][[co]][match(pb_info$sample, endo$analysis_sample)]
  pb_info[[co]] <- as.character(vals_co)
}
csv(pb_info, "DE/pseudobulk_sample_state_info.csv")
saveRDS(list(counts = pb_counts, metadata = pb_info), file.path(OUT, "DE/pseudobulk_counts.rds"))

run_de <- function(tr, st) {
  ix <- which(pb_info$state == st & pb_info$group %in% c("CK", tr) & pb_info$eligible)
  mi <- pb_info[ix, , drop = FALSE]
  n_ck <- sum(mi$group == "CK"); n_tr <- sum(mi$group == tr)
  if (n_ck < 1 || n_tr < 1) {
    note(paste("Skipping", tr, st, ": at least one group lacks pseudobulk samples meeting the cell threshold.")); return(NULL)
  }
  cc <- as.matrix(pb_counts[, ix, drop = FALSE])
  # Compute descriptive log2FC with equal sample weights so samples with more cells do not dominate.
  cpm0 <- sweep(cc, 2, colSums(cc), "/") * 1e6
  a <- rowMeans(cpm0[, mi$group == "CK", drop = FALSE])
  b <- rowMeans(cpm0[, mi$group == tr, drop = FALSE])
  detectable <- rowSums(cc) >= 10
  ans <- data.frame(gene = rownames(cc)[detectable], logFC = log2((b[detectable] + .5) / (a[detectable] + .5)),
    PValue = NA_real_, FDR = NA_real_, method = "descriptive_mean_CPM", stringsAsFactors = FALSE)
  if (BIOLOGICAL_REPLICATES_CONFIRMED && n_ck >= MIN_REPLICATES && n_tr >= MIN_REPLICATES) {
    design_md <- data.frame(treatment = factor(mi$group, levels = c("CK", tr)))
    covs <- character(0)
    for (nm in c("batch", "block")) {
      src <- if (nm == "batch") BATCH_COL else BLOCK_COL
      if (!is.null(src) && length(unique(mi[[src]])) > 1) {
        design_md[[nm]] <- factor(mi[[src]]); covs <- c(covs, nm)
      }
    }
    design <- model.matrix(reformulate(c(covs, "treatment")), design_md)
    if (qr(design)$rank < ncol(design) || nrow(design) <= ncol(design)) {
      note(paste("Falling back to descriptive comparison: ", tr, st, "Confounded design or no residual degrees of freedom; check batch/block."))
    } else {
      tested <- tryCatch({
        y <- edgeR::DGEList(counts = cc)
        keep_gene <- edgeR::filterByExpr(y, design = design)
        if (sum(keep_gene) < 10) stop("Fewer than 10 genes after filtering.")
        y <- y[keep_gene, , keep.lib.sizes = FALSE]
        y <- edgeR::calcNormFactors(y)
        y <- edgeR::estimateDisp(y, design)
        fit <- edgeR::glmQLFit(y, design)
        coef_tr <- grep("^treatment", colnames(design))
        stopifnot(length(coef_tr) == 1)
        qlf <- edgeR::glmQLFTest(fit, coef = coef_tr)
        tt <- edgeR::topTags(qlf, n = Inf, sort.by = "none")$table
        data.frame(gene = rownames(tt), tt, method = "edgeR_QL_pseudobulk", row.names = NULL)
      }, error = function(e) {
        note(paste("edgeR failed; retaining descriptive results: ", tr, st, conditionMessage(e))); NULL
      })
      if (!is.null(tested)) ans <- tested
    }
  } else if (BIOLOGICAL_REPLICATES_CONFIRMED) {
    note(paste("Insufficient independent replicates; descriptive comparison only: ", tr, st, "CK=", n_ck, "treatment=", n_tr))
  }
  ans$treatment <- tr; ans$state <- st; ans$n_CK <- n_ck; ans$n_treatment <- n_tr
  ans
}
de_list <- list()
for (tr in group_levels[-1]) for (st in state_names) {
  ans <- run_de(tr, st)
  if (!is.null(ans)) de_list[[paste(tr, st, sep = "_")]] <- ans
}
de <- dplyr::bind_rows(de_list)
if (nrow(de)) {
  de$FDR_global <- NA_real_
  tested <- is.finite(de$PValue)
  de$FDR_global[tested] <- p.adjust(de$PValue[tested], method = "BH")
  # Retain per-contrast FDR and use global BH adjustment across all successful contrasts for default DEG selection.
  de$significant <- is.finite(de$FDR_global) & de$FDR_global <= FDR_CUTOFF & abs(de$logFC) >= LOG2FC_CUTOFF
  de$direction <- ifelse(de$significant, ifelse(de$logFC > 0, "Up", "Down"),
                         ifelse(is.finite(de$PValue), "Not_significant", "Descriptive_only"))
  csv(de, "DE/all_treatments_states_vs_CK.csv")
  for (tr in group_levels[-1]) for (st in state_names) {
    sub <- de[de$treatment == tr & de$state == st, ]
    if (!nrow(sub)) next
    csv(sub, paste0("DE/", tr, "_", st, "_vs_CK.csv"))
  }
  inference <- de[is.finite(de$FDR_global), ]
  if (nrow(inference)) {
    inference$state <- factor(inference$state, state_names)
    inference$treatment <- factor(inference$treatment, group_levels[-1])
    save_plot(ggplot(inference, aes(logFC, -log10(pmax(FDR_global, 1e-300)), color = direction)) +
      geom_point(size = .35, alpha = .5) + facet_grid(treatment ~ state) + theme_bw() +
      scale_color_manual(values = c(Up = "#E64B35", Down = "#4DBBD5", Not_significant = "grey75")) +
      labs(y = "-log10(global BH FDR)", x = "log2FC (treatment / CK)"), "DE/volcano_all_states.pdf", 12, 10)
    deg_counts <- inference %>% dplyr::group_by(treatment, state) %>% dplyr::summarise(
      Up = sum(direction == "Up"), Down = sum(direction == "Down"), n_tested = dplyr::n(), .groups = "drop")
    csv(deg_counts, "DE/DEG_counts.csv")
  }
}

# 7. Response patterns across the three states.
# Three-column logFC heatmap: effect-size patterns are descriptive and do not establish significant interactions.
classify_pattern <- function(v) {
  if (all(v >= LOG2FC_CUTOFF)) return("All_states_activated")
  if (all(v <= -LOG2FC_CUTOFF)) return("All_states_repressed")
  if (min(v) <= -LOG2FC_CUTOFF && abs(min(v)) >= max(v)) return(paste0(state_names[which.min(v)], "_biased_repressed"))
  if (max(v) >= LOG2FC_CUTOFF) return(paste0(state_names[which.max(v)], "_biased_activated"))
  "Weak_or_mixed"
}
if (nrow(de)) for (tr in group_levels[-1]) {
  sub <- de[de$treatment == tr, ]
  if (!all(state_names %in% sub$state)) next
  common <- Reduce(intersect, lapply(state_names, function(s) sub$gene[sub$state == s]))
  if (length(common) < 2) next
  mat <- sapply(state_names, function(s) {
    tab <- sub[sub$state == s, ]; tab$logFC[match(common, tab$gene)]
  })
  rownames(mat) <- common; colnames(mat) <- state_names
  modes <- vapply(state_names, function(s) all(sub$method[sub$state == s] == "edgeR_QL_pseudobulk"), logical(1))
  formal <- all(modes)
  if (any(modes) && !formal) {
    note(paste(tr, "States mix inferential and descriptive estimates; skipping the three-state effect heatmap to avoid mixed scales.")); next
  }
  selected <- if (formal) unique(sub$gene[sub$significant]) else rownames(mat)[apply(abs(mat), 1, max) >= LOG2FC_CUTOFF]
  selected <- intersect(rownames(mat), selected)
  if (length(selected) < 2) next
  patterns <- apply(mat[selected, , drop = FALSE], 1, classify_pattern)
  pattern_table <- data.frame(gene = selected, mat[selected, , drop = FALSE], pattern = unname(patterns),
    evidence = if (formal) "DEG_in_at_least_one_state" else "descriptive_only", check.names = FALSE)
  csv(pattern_table, paste0("response/", tr, "_three_state_patterns.csv"))
  vis <- head(selected[order(apply(abs(mat[selected, , drop = FALSE]), 1, max), decreasing = TRUE)], 500)
  vis <- vis[order(patterns[vis])]
  heat(mat[vis, , drop = FALSE], paste0("response/", tr, "_three_state_logFC_heatmap.pdf"),
       paste(tr, "/ CK", if (!formal) "(descriptive)" else "(pseudobulk)"), split = factor(patterns[vis]))
  line_data <- dplyr::bind_rows(lapply(state_names, function(s)
    data.frame(gene = selected, state = s, logFC = mat[selected, s], pattern = unname(patterns))))
  profile <- line_data %>% dplyr::group_by(pattern, state) %>% dplyr::summarise(mean_logFC = mean(logFC), n_genes = dplyr::n(), .groups = "drop")
  profile$state <- factor(profile$state, state_names)
  csv(profile, paste0("response/", tr, "_pattern_mean_profiles.csv"))
  save_plot(ggplot(profile, aes(state, mean_logFC, group = pattern, color = pattern)) +
    geom_hline(yintercept = 0, color = "grey70") + geom_line() + geom_point() +
    facet_wrap(~pattern, scales = "free_y", ncol = 3) + theme_bw() +
    theme(legend.position = "none") + labs(title = paste(tr, "response patterns (descriptive classification)")),
    paste0("response/", tr, "_pattern_profiles.pdf"), 12, 8)
}
note("Significance in one state but not another does not establish a significant difference between states. Biased patterns describe logFC; formal specificity requires a treatment-by-state interaction test.")

# 8. Optional GO enrichment.
if (RUN_GO) {
  go_pkgs <- c("clusterProfiler", "AnnotationDbi", ORGDB_PACKAGE)
  go_ready <- all(vapply(go_pkgs, requireNamespace, logical(1), quietly = TRUE))
  if (!go_ready) {
    note("Skipping GO: clusterProfiler, AnnotationDbi or the configured OrgDb is missing; other analyses will be saved.")
  } else {
    orgdb <- getExportedValue(ORGDB_PACKAGE, ORGDB_PACKAGE)
    if (!ORGDB_KEYTYPE %in% AnnotationDbi::keytypes(orgdb)) {
      note(paste("Skipping GO: unsupported keytype ", ORGDB_KEYTYPE))
    } else {
      dbkeys <- AnnotationDbi::keys(orgdb, keytype = ORGDB_KEYTYPE)
      idmap <- if (is.null(GENE_ID_MAP_CSV)) data.frame(gene = rownames(expr), go_id = rownames(expr)) else
        read.csv(GENE_ID_MAP_CSV, stringsAsFactors = FALSE)
      if (!all(c("gene", "go_id") %in% colnames(idmap))) stop("GENE_ID_MAP_CSV must contain gene and go_id.")
      idmap <- unique(idmap[, c("gene", "go_id")]); idmap$gene <- as.character(idmap$gene); idmap$go_id <- as.character(idmap$go_id)
      idmap <- idmap[!is.na(idmap$gene) & !is.na(idmap$go_id) & idmap$go_id %in% dbkeys, ]
      csv(idmap, "GO/gene_ID_mapping_used.csv")
      go_results <- list(); mapping_audit <- list()
      run_go <- function(genes, universe, label) {
        universe <- unique(universe); genes <- intersect(unique(genes), universe)
        universe_ids <- unique(idmap$go_id[idmap$gene %in% universe])
        gene_ids <- unique(idmap$go_id[idmap$gene %in% genes])
        frac <- if (length(universe)) mean(universe %in% idmap$gene) else 0
        mapping_audit[[label]] <<- data.frame(label = label, input_genes = length(genes),
          mapped_genes = length(gene_ids), universe_genes = length(universe), universe_mapping_fraction = frac)
        if (frac < MIN_GO_MAPPING_FRACTION || length(gene_ids) < 5) return(NULL)
        res <- tryCatch(clusterProfiler::enrichGO(gene = gene_ids, universe = universe_ids,
          OrgDb = orgdb, keyType = ORGDB_KEYTYPE, ont = "BP", pAdjustMethod = "BH",
          pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 5, maxGSSize = 500, readable = FALSE),
          error = function(e) {note(paste("GO failed: ", label, conditionMessage(e))); NULL})
        if (is.null(res)) return(NULL)
        tab <- as.data.frame(res)
        if (!nrow(tab)) return(NULL)
        tab$comparison <- label
        go_results[[label]] <<- tab
      }
      # CK marker GO is exploratory, with detectable CK genes as the background.
      ck_universe <- rownames(expr)[Matrix::rowSums(expr[, ck_cells, drop = FALSE] > 0) >= 10]
      if (nrow(stage_markers)) for (s in state_names) {
        gg <- stage_markers$gene[stage_markers$cluster == s & stage_markers$p_val_adj <= FDR_CUTOFF]
        run_go(gg, ck_universe, paste0("CK_", s, "_markers_exploratory"))
      }
      # Treatment GO uses DEGs supported by replicate-level inference and the genes actually tested in each comparison as background.
      if (nrow(de)) for (tr in group_levels[-1]) for (s in state_names) {
        tab <- de[de$treatment == tr & de$state == s & is.finite(de$PValue), ]
        if (!nrow(tab)) next
        for (direction in c("Up", "Down")) run_go(tab$gene[tab$direction == direction], tab$gene,
          paste(tr, s, direction, sep = "_"))
      }
      csv(dplyr::bind_rows(mapping_audit), "GO/GO_mapping_audit.csv")
      all_go <- dplyr::bind_rows(go_results)
      if (nrow(all_go)) {
        all_go$FDR_all_GO_tests <- p.adjust(all_go$pvalue, "BH")
        csv(all_go, "GO/GO_BP_all_tested.csv")
        sig_go <- all_go[is.finite(all_go$p.adjust) & all_go$p.adjust <= FDR_CUTOFF, ]
        csv(sig_go, "GO/GO_BP_significant_within_each_list.csv")
        if (nrow(sig_go)) {
          top_go <- sig_go %>% dplyr::group_by(comparison) %>% dplyr::arrange(p.adjust, .by_group = TRUE) %>% dplyr::slice_head(n = 5) %>% dplyr::ungroup()
          top_go$neglog10FDR <- -log10(pmax(top_go$p.adjust, 1e-300))
          top_go$term <- vapply(top_go$Description, function(v) paste(strwrap(v, 45), collapse = "\n"), character(1))
          save_plot(ggplot(top_go, aes(neglog10FDR, reorder(term, neglog10FDR), size = Count, color = comparison)) +
            geom_point() + facet_wrap(~comparison, scales = "free_y", ncol = 3) + theme_bw() +
            theme(legend.position = "none") + labs(x = "-log10(BH FDR within gene list)", y = NULL),
            "GO/GO_BP_top_terms.pdf", 17, max(7, 3 * ceiling(length(unique(top_go$comparison)) / 3)))
        }
      } else note("No GO results; inspect GO_mapping_audit.csv. ID mapping failure does not establish absence of biological function.")
    }
  }
}

# 9. Save outputs and interpretation notes.
assert_original_umap(endo)
assert_original_umap(ck)
csv(data.frame(object = c("Endodermis", "CK"), original_reduction = UMAP_REDUCTION,
  coordinates_identical = TRUE, tolerance = 0), "UMAP/original_UMAP_preservation_check.csv")
saveRDS(list(pca = pca_fit, hvg = hvg, center = mu, scale = sdv, slingshot = ss,
  lineage_id = lid, state_mode = STATE_MODE, states = state_names, cuts = cuts,
  cluster_to_state = CLUSTER_TO_STATE, root = ROOT_CLUSTER,
  biological_direction_confirmed = BIOLOGICAL_DIRECTION_CONFIRMED,
  mapping_k = k, mapping_score_cutoff = MIN_MAP_SCORE, mapping_distance_cutoff = dist_cut,
  CK_cells = ck_cells, CK_states = ck_state, CK_pseudotime = pt), file.path(OUT, "CK/CK_reference_model.rds"))
if (SAVE_SEURAT) {
  saveRDS(ck, file.path(OUT, "CK/CK_Endodermis_reference.rds"))
  saveRDS(endo, file.path(OUT, "Endodermis_all_groups_CK_states.rds"))
}
capture.output(sessionInfo(), file = file.path(OUT, "sessionInfo.txt"))
config_vars <- c("UMAP_REDUCTION", "UMAP_DIMS", "PERIOD_IDS", "UMAP_POINT_SIZE", "UMAP_LABEL_SIZE", "UMAP_LABEL_REPEL", "SAVE_UMAP_PNG", "MERGED_RDS", "OUTPUT_ROOT", "CELLTYPE_COL", "TREATMENT_COL", "SAMPLE_COL",
  "BIOLOGICAL_REPLICATES_CONFIRMED", "BATCH_COL", "BLOCK_COL", "TREATMENT_LABELS",
  "ENDODERMIS_LABELS", "SEED", "N_HVG", "N_PCS", "RESOLUTION", "RESOLUTION_GRID",
  "ROOT_CLUSTER", "BIOLOGICAL_DIRECTION_CONFIRMED", "LINEAGE_ID", "STATE_MODE", "CLUSTER_TO_STATE",
  "EARLY_MARKER_IDS", "LATE_MARKER_IDS", "EXCLUDE_FROM_PCA", "K_NEIGHBORS", "MIN_MAP_SCORE",
  "MAP_DISTANCE_QUANTILE", "MIN_CK_CELLS", "MIN_CELLS_PER_PB", "MIN_REPLICATES", "FDR_CUTOFF",
  "LOG2FC_CUTOFF", "RUN_GO", "ORGDB_PACKAGE", "ORGDB_KEYTYPE", "GENE_ID_MAP_CSV")
saveRDS(mget(config_vars), file.path(OUT, "analysis_configuration.rds"))
note(paste0("Completed. Output directory: ", OUT))
note("Interpret results with CK subcluster markers, resolution/sample distributions, trajectory and developmental-marker direction, treatment mapping quality, sample-level composition and within-state pseudobulk responses. This script alone cannot establish a true differentiation process.")


cat("\nCompleted. Results saved in:", OUT, "\n")

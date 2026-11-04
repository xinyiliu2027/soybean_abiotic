# Source: Figure 2 document, extracted text lines 1-474
# Standalone script requiring a merged expression RDS and CK RData.
# See README.md and docs/KNOWN_LIMITATIONS.md before running.

suppressMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(muscat)
  library(edgeR)
  library(limma)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(cowplot)
})

## ==============================
# 0. Paths and parameters.
## ==============================
MERGED_RDS <- file.path(Sys.getenv("SOYBEAN_PROJECT_ROOT", unset = "data"), "seurat_obj_with_celltype.rds")
CK_RDATA   <- file.path(Sys.getenv("SOYBEAN_PROJECT_ROOT", unset = "data"), "data_files/Gm_CK_integrated0.8.RData")  # Source of REF cells.

OUT_DIR <- file.path(Sys.getenv("SOYBEAN_PROJECT_ROOT", unset = "data"), "muscat_results_ck_ref_paper")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Thresholds described in the source document's cited paper.
# FDR <= 0.05; fold-change magnitude >= 1.5; detection frequency >= 10% in either condition.
FDR_CUT <- 0.05
LFC_CUT <- log2(1.5)   # ≈ 0.5849625
FRQ_CUT <- 0.10

# pbDS settings: edgeR with pbDS and aggregateData using summed counts.
PBDS_METHOD <- "edgeR"
PBDS_FILTER <- "none"

# Internal pbDS filtering: minimize additional filtering to follow the source protocol.
PBDS_MIN_CELLS_STRICT <- 10
PBDS_MIN_CELLS_LOOSE  <- 1

cat("========================================\n")
cat("Paper-like parameter settings:\n")
cat("FDR_CUT =", FDR_CUT, "\n")
cat("LFC_CUT =", LFC_CUT, "(fold change ≥", 2^LFC_CUT, ")\n")
cat("FRQ_CUT =", FRQ_CUT, "\n")
cat("PBDS_METHOD =", PBDS_METHOD, "\n")
cat("PBDS_FILTER =", PBDS_FILTER, "\n")
cat("REF source = CK RData file\n")
cat("========================================\n\n")

## ==============================
# 1. Compatibility helpers for reading Seurat counts.
## ==============================
get_counts_Seurat <- function(obj, assay = "RNA") {
  out <- tryCatch({
    LayerData(obj, assay = assay, layer = "counts")
  }, error = function(e) {
    tryCatch({
      GetAssayData(obj, assay = assay, layer = "counts")
    }, error = function(e2) {
      GetAssayData(obj, assay = assay, slot = "counts")
    })
  })
  return(out)
}

pick_assay_with_counts <- function(obj) {
  cand <- unique(c("RNA", "SCT", DefaultAssay(obj)))
  for (a in cand) {
    ok <- tryCatch({
      x <- get_counts_Seurat(obj, assay = a)
      !is.null(x) && nrow(x) > 0 && ncol(x) > 0
    }, error = function(e) FALSE)
    if (ok) return(a)
  }
  stop("No usable counts found in RNA, SCT or DefaultAssay.")
}

## ==============================
# 2. Read the merged object for non-REF treatments only.
## ==============================
cat("Reading merged Seurat object: ", MERGED_RDS, "\n")
merged <- readRDS(MERGED_RDS)

stopifnot(all(c("Group1","orig.ident","cell_type") %in% colnames(merged@meta.data)))

cat("merged cells:", ncol(merged), " genes:", nrow(merged), "\n")
cat("\nGroup1 distribution:\n"); print(table(merged$Group1))
cat("\ncell_type distribution:\n"); print(sort(table(merged$cell_type), decreasing = TRUE))

# Discard all REF cells from the merged object; use the separate CK reference instead.
merged_nonref <- subset(merged, subset = Group1 != "REF")
cat("\n[INFO] merged_nonref cells:", ncol(merged_nonref), "\n")
cat("[INFO] merged_nonref Group1:\n"); print(table(merged_nonref$Group1))

merged_nonref$sample_id_muscat <- as.character(merged_nonref$orig.ident)
merged_nonref$cluster_id_use   <- as.character(merged_nonref$cell_type)

## ==============================
# 3. Read CK RData as the source of REF1 and REF2.
## ==============================
cat("\nReading CK RData file: ", CK_RDATA, "\n")

loaded_env <- new.env()
load(CK_RDATA, envir = loaded_env)

cat("[RData] Loaded objects:\n")
print(ls(envir = loaded_env))

# Find the first Seurat object.
obj_names <- ls(envir = loaded_env)
seurat_obj_name <- NULL
for (nm in obj_names) {
  obj <- get(nm, envir = loaded_env)
  if (inherits(obj, "Seurat")) { seurat_obj_name <- nm; break }
}
if (is.null(seurat_obj_name)) stop("No Seurat object found in the RData file.")

cat("[RData] Seurat object found: ", seurat_obj_name, "\n")
ck_seurat <- get(seurat_obj_name, envir = loaded_env)

cat("[RData] Seurat dimensions: ", ncol(ck_seurat), " cells, ", nrow(ck_seurat), " genes\n")
cat("[RData] Metadata columns: ", paste(colnames(ck_seurat@meta.data), collapse = ", "), "\n")

# Sample field.
sample_field <- NULL
for (field in c("orig.ident", "sample", "Sample", "group", "Group", "condition")) {
  if (field %in% colnames(ck_seurat@meta.data)) { sample_field <- field; break }
}
if (is.null(sample_field)) stop("No sample field found (orig.ident/sample/Sample/group/Group/condition).")
cat("[RData] Sample field: ", sample_field, "\n")
print(table(ck_seurat@meta.data[[sample_field]]))

# Cluster field.
cluster_field <- NULL
for (field in c("seurat_clusters", "cluster", "ClusterID", "cell_type", "CellType")) {
  if (field %in% colnames(ck_seurat@meta.data)) { cluster_field <- field; break }
}
if (is.null(cluster_field)) stop("No cluster field found (seurat_clusters/cluster/ClusterID/cell_type/CellType).")
cat("[RData] Cluster field: ", cluster_field, "\n")
print(table(ck_seurat@meta.data[[cluster_field]]))

## ==============================
# 4. Map numeric cluster labels to cell-type annotations when needed.
## ==============================
cluster_info <- as.character(ck_seurat@meta.data[[cluster_field]])
unique_clusters <- unique(cluster_info)
cat("[RData] Unique cluster values: ", paste(sort(unique_clusters), collapse = ", "), "\n")

if (all(grepl("^\\d+$", unique_clusters))) {
  cat("[RData] Numeric clusters detected; mapping cell types...\n")
  cluster_annotation <- c(
    "0"  = "Cortex","1"="Cortex","2"="QC","3"="Cortex","4"="Cortex",
    "5"  = "Procambium","6"="Dividing Cell","7"="Root Cap","8"="Meristem",
    "9"  = "Endodermis","10"="Pericycle","11"="Epidermis","12"="Cortex",
    "13" = "Pericycle","14"="Procambium","15"="Meristem","16"="Cortex",
    "17" = "Phloem","18"="Xylem"
  )
  cluster_mapped <- unname(cluster_annotation[as.character(cluster_info)])
  cluster_mapped[is.na(cluster_mapped)] <- "Unknown"
  print(sort(table(cluster_mapped), decreasing = TRUE))
} else {
  cluster_mapped <- cluster_info
  cat("[RData] Using existing cell-type labels.\n")
}

## ==============================
# 5. Create SCE objects for REF1 and REF2.
## ==============================
cat("\nBuilding SCE objects for REF1/REF2...\n")

sample_info <- as.character(ck_seurat@meta.data[[sample_field]])
unique_samples <- unique(sample_info)
cat("[RData] Unique sample names: ", paste(unique_samples, collapse = ", "), "\n")

ck_assay <- pick_assay_with_counts(ck_seurat)
cat("[RData] Counts assay: ", ck_assay, "\n")
ck_counts <- get_counts_Seurat(ck_seurat, assay = ck_assay)

# Identify CK1/CK2; the original fallback uses the first two samples. Verify their identities.
ck1_pattern <- "^CK1$|^ck1$|^CK_1$|^ck_1$|^control1$|^Control1$"
ck2_pattern <- "^CK2$|^ck2$|^CK_2$|^ck_2$|^control2$|^Control2$"
ck1_idx <- grep(ck1_pattern, sample_info, ignore.case = TRUE)
ck2_idx <- grep(ck2_pattern, sample_info, ignore.case = TRUE)

if (length(ck1_idx) == 0 || length(ck2_idx) == 0) {
  cat("[WARN] Standard CK1/CK2 names not found; using the first two samples as REF1/REF2.\n")
  if (length(unique_samples) < 2) stop("Fewer than two RData samples; cannot construct REF1/REF2.")
  ck1_name <- unique_samples[1]
  ck2_name <- unique_samples[2]
  ck1_idx <- which(sample_info == ck1_name)
  ck2_idx <- which(sample_info == ck2_name)
  cat("[INFO] REF1=", ck1_name, "；REF2=", ck2_name, "\n", sep = "")
}

cat("[INFO] REF1 cells:", length(ck1_idx), "\n")
cat("[INFO] REF2 cells:", length(ck2_idx), "\n")

counts_ref1 <- ck_counts[, ck1_idx, drop = FALSE]
sce_ref1 <- SingleCellExperiment(
  assays = list(counts = counts_ref1),
  colData = data.frame(
    group_id   = factor(rep("REF",  ncol(counts_ref1))),
    sample_id  = factor(rep("REF1", ncol(counts_ref1))),
    cluster_id = factor(cluster_mapped[ck1_idx]),
    row.names  = colnames(counts_ref1),
    stringsAsFactors = FALSE
  )
)

counts_ref2 <- ck_counts[, ck2_idx, drop = FALSE]
sce_ref2 <- SingleCellExperiment(
  assays = list(counts = counts_ref2),
  colData = data.frame(
    group_id   = factor(rep("REF",  ncol(counts_ref2))),
    sample_id  = factor(rep("REF2", ncol(counts_ref2))),
    cluster_id = factor(cluster_mapped[ck2_idx]),
    row.names  = colnames(counts_ref2),
    stringsAsFactors = FALSE
  )
)

cat("\n[CHECK] REF1 cluster_id:\n"); print(sort(table(colData(sce_ref1)$cluster_id), decreasing = TRUE))
cat("\n[CHECK] REF2 cluster_id:\n"); print(sort(table(colData(sce_ref2)$cluster_id), decreasing = TRUE))

## ==============================
## 6) merged_nonref -> SCE
## ==============================
cat("\nBuilding merged_nonref SCE...\n")
m_assay  <- pick_assay_with_counts(merged_nonref)
cat("[MERGED_NONREF] Counts assay:", m_assay, "\n")
m_counts <- get_counts_Seurat(merged_nonref, assay = m_assay)

sce_nonref <- SingleCellExperiment(
  assays = list(counts = m_counts),
  colData = data.frame(
    group_id   = factor(as.character(merged_nonref$Group1)),
    sample_id  = factor(as.character(merged_nonref$sample_id_muscat)),
    cluster_id = factor(as.character(merged_nonref$cluster_id_use)),
    row.names  = colnames(m_counts),
    stringsAsFactors = FALSE
  )
)
cat("[MERGED_NONREF] cells:", ncol(sce_nonref), " genes:", nrow(sce_nonref), "\n")

## ==============================
# 7. Align genes and combine the complete SCE.
## ==============================
cat("\nAligning genes and merging SCE objects...\n")
common_genes <- Reduce(intersect, list(rownames(sce_nonref), rownames(sce_ref1), rownames(sce_ref2)))
cat("Shared genes:", length(common_genes), "\n")
if (length(common_genes) < 1000) stop("Too few shared genes; check gene identifier consistency.")

sce_all <- cbind(
  sce_nonref[common_genes, ],
  sce_ref1[common_genes, ],
  sce_ref2[common_genes, ]
)

cat("Merged SCE: cells =", ncol(sce_all), ", genes =", nrow(sce_all), "\n")
cat("\nMerged group_id:\n"); print(table(colData(sce_all)$group_id))

## experiment_info（sample-level）
ei <- as.data.frame(colData(sce_all)) %>%
  as_tibble(rownames = "cell") %>%
  dplyr::select(sample_id, group_id) %>%
  distinct() %>%
  arrange(group_id, sample_id)

n_cells <- table(colData(sce_all)$sample_id)
ei$n_cells <- as.numeric(n_cells[match(ei$sample_id, names(n_cells))])
metadata(sce_all)$experiment_info <- as.data.frame(ei)

cat("\n[experiment_info]:\n")
print(metadata(sce_all)$experiment_info)

## ==============================
# 8. Guarded pbDS call with automatic retries.
## ==============================
call_pbDS_safe <- function(pb, design, contrast,
                           method = "edgeR",
                           filter = "none",
                           min_cells = NULL,
                           min_counts = NULL,
                           verbose = TRUE) {
  fml <- names(formals(muscat::pbDS))
  args <- list(pb = pb, design = design, contrast = contrast,
               method = method, filter = filter, verbose = verbose)
  if (!is.null(min_cells)  && "min_cells"  %in% fml) args$min_cells  <- min_cells
  if (!is.null(min_counts) && "min_counts" %in% fml) args$min_counts <- min_counts
  do.call(muscat::pbDS, args)
}

run_pbDS_with_retry <- function(pb, mm, contrast,
                                method = PBDS_METHOD,
                                filter = PBDS_FILTER,
                                min_cells_strict = PBDS_MIN_CELLS_STRICT,
                                min_cells_loose  = PBDS_MIN_CELLS_LOOSE) {
  msg <- NULL
  res <- tryCatch({
    call_pbDS_safe(pb, mm, contrast, method = method, filter = filter,
                   min_cells = min_cells_strict, min_counts = 1, verbose = TRUE)
  }, error = function(e) { msg <<- conditionMessage(e); NULL })
  
  if (!is.null(res)) return(list(res = res, used_min_cells = min_cells_strict, used_min_counts = 1))
  
  cat("\n[WARN] pbDS failed with strict thresholds:\n  ", msg, "\n", sep = "")
  cat("[RETRY] min_cells=", min_cells_loose, "；min_counts=0\n", sep = "")
  
  msg2 <- NULL
  res2 <- tryCatch({
    call_pbDS_safe(pb, mm, contrast, method = method, filter = filter,
                   min_cells = min_cells_loose, min_counts = 0, verbose = TRUE)
  }, error = function(e) { msg2 <<- conditionMessage(e); NULL })
  
  if (!is.null(res2)) return(list(res = res2, used_min_cells = min_cells_loose, used_min_counts = 0))
  
  cat("\n[ERROR] pbDS also failed with relaxed thresholds:\n  ", msg2, "\n", sep = "")
  return(list(res = NULL, used_min_cells = NA, used_min_counts = NA))
}

## ==============================
# 9. Run pbDS for each treatment versus REF and apply the specified thresholds.
## ==============================
treatments <- sort(setdiff(unique(as.character(colData(sce_all)$group_id)), "REF"))
cat("\nComparisons to run (treatment versus REF):\n")
print(treatments)

all_results <- list()

for (treatment in treatments) {
  cat("\n========================================\n")
  cat(">>> Analyzing: ", treatment, " vs REF\n", sep = "")
  cat("========================================\n")
  
  # 9.1 Subset the current treatment and REF.
  keep <- colData(sce_all)$group_id %in% c("REF", treatment)
  sce_sub <- sce_all[, keep]
  
  colData(sce_sub)$group_id   <- droplevels(colData(sce_sub)$group_id)
  colData(sce_sub)$sample_id  <- droplevels(colData(sce_sub)$sample_id)
  colData(sce_sub)$cluster_id <- droplevels(colData(sce_sub)$cluster_id)
  
  cat("Subset cells:", ncol(sce_sub), "\n")
  cat("Subset group_id:\n"); print(table(colData(sce_sub)$group_id))
  
  ## 9.2 prepSCE
  cat("\nRunning prepSCE...\n")
  sce_sub <- prepSCE(
    sce_sub,
    kid  = "cluster_id",
    gid  = "group_id",
    sid  = "sample_id",
    drop = TRUE
  )
  
  ei_sub <- metadata(sce_sub)$experiment_info
  cat("\nexperiment_info:\n"); print(ei_sub)
  cat("\nSamples per group:\n"); print(table(ei_sub$group_id))
  
  # 9.3 Sum raw counts into pseudobulk profiles following the source protocol.
  cat("\nAggregating pseudobulk...\n")
  pb <- aggregateData(
    sce_sub,
    assay = "counts",
    fun   = "sum",
    by    = c("cluster_id", "sample_id")
  )
  cat("Pseudobulk complete: ", length(assayNames(pb)), " cell types\n", sep = "")
  
  # Record cluster-by-sample cell counts for QC.
  nc_mat <- table(colData(sce_sub)$cluster_id, colData(sce_sub)$sample_id)
  nc_mat <- nc_mat[assayNames(pb), colnames(pb), drop = FALSE]
  metadata(pb)$n_cells <- nc_mat
  
  # Align group_id and sample_id in pseudobulk colData.
  ei2 <- metadata(sce_sub)$experiment_info
  colData(pb)$group_id  <- ei2$group_id[match(colnames(pb), ei2$sample_id)]
  colData(pb)$sample_id <- colnames(pb)
  
  # 9.4 Construct the design matrix and contrast.
  pb_groups <- factor(as.character(colData(pb)$group_id), levels = c("REF", treatment))
  mm <- model.matrix(~0 + pb_groups)
  colnames(mm) <- levels(pb_groups)
  rownames(mm) <- colData(pb)$sample_id
  
  contrast_name <- paste0(treatment, "-REF")
  contrast <- limma::makeContrasts(contrasts = contrast_name, levels = mm)
  
  cat("\nDesign matrix (mm):\n"); print(mm)
  cat("\nContrast matrix:\n"); print(contrast)
  
  # 9.5 Run edgeR through pbDS following the source protocol.
  cat("\nRunning pbDS with automatic retries...\n")
  pbds_out <- run_pbDS_with_retry(pb, mm, contrast)
  if (is.null(pbds_out$res)) {
    cat("[SKIP] ", treatment, " versus REF: pbDS failed; skipping.\n", sep = "")
    next
  }
  res <- pbds_out$res
  cat("[OK] pbDS succeeded: min_cells=", pbds_out$used_min_cells,
      "；min_counts=", pbds_out$used_min_counts, "\n", sep="")
  
  # 9.6 Extract results, including detection frequencies.
  cat("\nExtracting results...\n")
  tb_full <- resDS(sce_sub, res, bind = "row", cpm = FALSE, frq = TRUE)
  cat("Total result rows:", nrow(tb_full), "\n")
  
  # 9.7 Use FDR-adjusted P values, preferring local adjustment as in the source code.
  if ("p_adj.loc" %in% colnames(tb_full)) {
    p_adj_col <- "p_adj.loc"
    cat("Filtering with p_adj.loc.\n")
  } else if ("p_adj.glb" %in% colnames(tb_full)) {
    p_adj_col <- "p_adj.glb"
    cat("[WARN] p_adj.loc missing; using p_adj.glb.\n")
  } else {
    stop("tb_full contains neither p_adj.loc nor p_adj.glb.")
  }
  
  # 9.8 Apply the specified significance and effect-size thresholds.
  # FDR <= 0.05 and abs(log2FC) >= log2(1.5).
  tb_sig <- tb_full %>%
    dplyr::filter(
      .data[[p_adj_col]] <= FDR_CUT,
      abs(logFC) >= LFC_CUT
    )
  cat("After ", p_adj_col, " + abs(logFC) filtering: ", nrow(tb_sig), " rows\n", sep = "")
  
  # 9.9 Require detection frequency >= 10% in either condition.
  ref_frq_col <- "REF.frq"
  trt_frq_col <- paste0(treatment, ".frq")
  if (ref_frq_col %in% colnames(tb_sig) && trt_frq_col %in% colnames(tb_sig)) {
    tb_sig <- tb_sig %>%
      dplyr::filter(.data[[ref_frq_col]] >= FRQ_CUT | .data[[trt_frq_col]] >= FRQ_CUT)
    cat("After detection-frequency filtering:", nrow(tb_sig), "rows\n")
  } else {
    cat("[WARN] No frq column found; skipping frequency filtering (not recommended).\n")
  }
  
  # 9.10 Label Up/Down; positive logFC means higher expression in treatment than REF.
  tb_sig <- tb_sig %>%
    dplyr::mutate(
      regulation = dplyr::case_when(
        logFC >=  LFC_CUT ~ "Up",
        logFC <= -LFC_CUT ~ "Down",
        TRUE              ~ "Not_DE"
      ),
      comparison = paste0(treatment, "_vs_REF")
    )
  
  cat("\n========================================\n")
  cat("Paper-like filtered results\n")
  cat("========================================\n")
  cat("Significant DEG rows:", nrow(tb_sig), "\n")
  cat("Unique significant DEG genes:", dplyr::n_distinct(tb_sig$gene), "\n")
  
  # 9.11 Summarize DEGs by cell type.
  deg_summary <- tb_sig %>%
    dplyr::filter(regulation %in% c("Up", "Down")) %>%
    dplyr::group_by(cluster_id, regulation) %>%
    dplyr::summarise(n_genes = dplyr::n(), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = regulation, values_from = n_genes, values_fill = 0) %>%
    dplyr::mutate(
      Up = ifelse(is.na(Up), 0, Up),
      Down = ifelse(is.na(Down), 0, Down),
      Total = Up + Down
    ) %>%
    arrange(desc(Total))
  
  cat("\nDEG summary by cell type:\n")
  print(deg_summary)
  
  # 9.12 Save results.
  treatment_dir <- file.path(OUT_DIR, paste0(treatment, "_vs_REF"))
  dir.create(treatment_dir, showWarnings = FALSE, recursive = TRUE)
  
  write.csv(tb_full, file.path(treatment_dir, "DEG_full_results_muscat.csv"), row.names = FALSE)
  write.csv(tb_sig,  file.path(treatment_dir, "DEG_filtered_paper_thresholds.csv"), row.names = FALSE)
  write.csv(deg_summary, file.path(treatment_dir, "DEG_summary_muscat.csv"), row.names = FALSE)
}

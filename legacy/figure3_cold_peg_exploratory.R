# Source: Figure 3 document, extracted text lines 1-284
# Early exploratory fragment ending at the Step 7 heading; not the main workflow.
# See README.md and docs/KNOWN_LIMITATIONS.md before running.

# ============================================================================
# Endodermis single-cell pseudotime analysis restricted to Low_Temp and PEG.
# Planned panels: trajectory UMAP, density, stage UMAP, heatmap/GO, DEG bars and Venn diagram.
# ============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(Matrix)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(ggplot2)
  library(patchwork)
  library(clusterProfiler)
  library(org.G.max.eg.db) # Soybean annotation database.
  library(slingshot)
  library(SingleCellExperiment)
  library(stringr)
  library(scales)
  library(ggridges)
  library(ggVennDiagram)
})

# ============================================================================
# Step 0: Set the output directory.
# ============================================================================
output_dir <- file.path(Sys.getenv("SOYBEAN_PROJECT_ROOT", unset = "data"), "Endodermis_Cold_PEG_Analysis")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
setwd(output_dir)

# ============================================================================
# Step 1: Read data and retain only Endodermis cells under Low_Temp or PEG.
# ============================================================================
seurat_raw <- readRDS(file.path(Sys.getenv("SOYBEAN_PROJECT_ROOT", unset = "data"), "seurat_obj_with_celltype.rds"))
celltype_col <- "cell_type"
treatment_col <- "treatment"

# Define the cells to retain.
target_treatments <- c("Low_Temp", "PEG", "Cold") # Support alternative naming conventions.
available_trts <- unique(seurat_raw@meta.data[[treatment_col]])
valid_trts <- intersect(target_treatments, available_trts)

# Select Endodermis cells from the requested treatment groups.
keep_cells <- colnames(seurat_raw)[
  seurat_raw@meta.data[[celltype_col]] == "Endodermis" &
    seurat_raw@meta.data[[treatment_col]] %in% valid_trts
]

cat("Total cells after filtering:", length(keep_cells), "\n")
seurat_sub <- subset(seurat_raw, cells = keep_cells)

# Check the UMAP embedding.
if (!"umap" %in% names(seurat_sub@reductions)) {
  stop("UMAP embedding missing; run RunUMAP() first.")
}

# ============================================================================
# Step 2: Infer trajectories and pseudotime with Slingshot.
# ============================================================================
cat("Running Slingshot...\n")
sce <- as.SingleCellExperiment(seurat_sub, assay = "RNA")
reducedDim(sce, "UMAP") <- seurat_sub@reductions$umap@cell.embeddings

# Use k-means clusters in UMAP space as Slingshot input because multiple cell-type labels are unavailable.
set.seed(123)
kclusters <- kmeans(reducedDim(sce, "UMAP"), centers = 4)$cluster
colData(sce)$kmeans <- factor(kclusters)

sce_ss <- slingshot(sce, clusterLabels = "kmeans", reducedDim = "UMAP", approx_points = 150)

# Compute combined pseudotime by averaging across lineages.
pt_mat <- slingPseudotime(sce_ss)
pseudotime <- rowMeans(pt_mat, na.rm = TRUE)
pseudotime[is.nan(pseudotime)] <- NA
seurat_sub$pseudotime <- pseudotime[colnames(seurat_sub)]

# Remove cells with missing values.
valid_cells <- colnames(seurat_sub)[!is.na(seurat_sub$pseudotime)]
seurat_clean <- subset(seurat_sub, cells = valid_cells)
cat("Cells with valid pseudotime:", ncol(seurat_clean), "\n")

# ============================================================================
# Step 3: Define stages P1 (early), P2 (middle) and P3 (late).
# ============================================================================
# Split pseudotime at its 33rd and 66th percentiles.
quantiles <- quantile(seurat_clean$pseudotime, probs = c(0, 0.33, 0.66, 1), na.rm = TRUE)
seurat_clean$Period <- cut(seurat_clean$pseudotime, 
                           breaks = quantiles, 
                           labels = c("P1", "P2", "P3"), 
                           include.lowest = TRUE)

# ============================================================================
# Step 4: Plot pseudotime trajectories and stage assignments on UMAP.
# ============================================================================

# 4.1 Pseudotime trajectory plot.
p_pt_umap <- DimPlot(seurat_clean, reduction = "umap", group.by = "pseudotime", 
                     cols = rev(brewer.pal(11, "Spectral")), pt.size = 0.5) +
  scale_color_gradientn(colors = rev(brewer.pal(11, "Spectral")), na.value = "grey80") +
  labs(title = "Endodermis Pseudotime") + theme_minimal()

# Add Slingshot trajectory lines manually.
umap_emb <- seurat_clean@reductions$umap@cell.embeddings
# Extract Slingshot curves.
curves <- slingCurves(sce_ss)
curve_df <- do.call(rbind, lapply(curves, function(crv) {
  data.frame(UMAP_1 = crv$s[crv$ord, 1], UMAP_2 = crv$s[crv$ord, 2], Lineage = crv$lineage)
}))

p_pt_umap_lines <- p_pt_umap + 
  geom_path(data = curve_df, aes(x = UMAP_1, y = UMAP_2, group = Lineage), 
            color = "black", size = 1.2, inherit.aes = FALSE)

ggsave("01_Pseudotime_UMAP.pdf", p_pt_umap_lines, width = 6, height = 5)

# 4.2 UMAP colored by stage (P1, P2, P3).
p_period_umap <- DimPlot(seurat_clean, reduction = "umap", group.by = "Period", 
                         cols = c("#E41A1C", "#377EB8", "#4DAF4A"), pt.size = 0.5) +
  labs(title = "Endodermis Periods (P1-P3)") + theme_minimal()
ggsave("02_Period_UMAP.pdf", p_period_umap, width = 6, height = 5)

# ============================================================================
# Step 5: Plot pseudotime density/ridges.
# ============================================================================
pt_df <- data.frame(
  Pseudotime = seurat_clean$pseudotime,
  Treatment = seurat_clean@meta.data[[treatment_col]]
)

p_ridges <- ggplot(pt_df, aes(x = Pseudotime, y = Treatment, fill = Treatment)) +
  geom_density_ridges(alpha = 0.7, scale = 1.2, rel_min_height = 0.01) +
  scale_fill_manual(values = c("Low_Temp" = "#377EB8", "Cold" = "#377EB8", "PEG" = "#4DAF4A")) +
  labs(title = "Pseudotime Distribution by Treatment", x = "Pseudotime", y = "") +
  theme_ridges() + theme(legend.position = "none")
ggsave("03_Pseudotime_Ridge.pdf", p_ridges, width = 7, height = 4)

# ============================================================================
# Step 6: Generate expression heatmaps and GO enrichment panels.
# ============================================================================
# Adapt the original supplied function with minor changes.
parse_ratio <- function(x) {
  if (is.numeric(x)) return(x)
  vapply(as.character(x), function(s) {
    if (grepl("/", s)) {
      p <- as.numeric(strsplit(s, "/")[[1]])
      if (p[2] != 0) return(p[1] / p[2])
    }
    suppressWarnings(as.numeric(s))
  }, numeric(1))
}

run_heatmap_go <- function(seurat_obj, n_genes=500, n_bins=50, min_cells_expr=10, 
                           smooth_window=5, go_max_terms=5, go_padj_cutoff=0.05) {
  
  pt <- seurat_obj$pseudotime
  cells_ordered <- names(sort(pt[!is.na(pt)]))
  pt_vals <- pt[cells_ordered]
  
  # Extract the expression matrix.
  expr <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")[, cells_ordered, drop=FALSE]
  expressed_genes <- names(which(Matrix::rowSums(expr > 0) >= min_cells_expr))
  
  # Compute Spearman correlations with pseudotime.
  expr_dense <- as.matrix(expr[expressed_genes, , drop=FALSE])
  gene_cors <- fast_spearman_vs_pt(expr_dense, pt_vals)
  
  # Select the top genes.
  top_genes <- head(names(sort(abs(gene_cors), decreasing=TRUE)), n_genes)
  
  # Bin and smooth gene expression.
  bin_idx <- cut(seq_along(cells_ordered), breaks=n_bins, labels=FALSE)
  expr_binned <- t(apply(expr[top_genes,,drop=FALSE], 1, function(x) tapply(x, bin_idx, mean, na.rm=TRUE)))
  expr_smooth <- t(apply(expr_binned, 1, function(x) {
    if (length(x) > smooth_window) {
      y <- as.numeric(stats::filter(x, rep(1/smooth_window, smooth_window), sides=2))
      y[is.na(y)] <- x[is.na(y)]; y
    } else x
  }))
  rownames(expr_smooth) <- top_genes
  expr_scaled <- t(scale(t(expr_smooth)))
  expr_scaled[!is.finite(expr_scaled)] <- 0
  
  # Cluster by peak time (K = 3).
  K <- 3
  peak_time <- apply(expr_scaled, 1, which.max)
  breaks <- c(0, floor(n_bins/3), floor(2*n_bins/3), n_bins)
  cl_raw <- as.integer(as.character(cut(peak_time, breaks=breaks, labels=1:3, include.lowest=TRUE)))
  cl_raw[is.na(cl_raw)] <- 2
  names(cl_raw) <- rownames(expr_scaled)
  
  # Order genes and extract modules.
  ord_df <- data.frame(gene=rownames(expr_scaled), C=cl_raw, pk=peak_time) %>% arrange(C, pk)
  gene_order <- ord_df$gene
  expr_scaled_ord <- expr_scaled[gene_order,,drop=FALSE]
  cluster_factor <- factor(paste0("C", cl_raw[gene_order]), levels=paste0("C",1:K))
  
  # GO enrichment.
  go_df_list <- vector("list", K)
  gid_keys <- tryCatch(keys(org.G.max.eg.db, keytype="GID"), error=function(e) character(0))
  
  for (kk in 1:K) {
    genes_k <- gene_order[cl_raw[gene_order] == kk]
    mapped <- genes_k[genes_k %in% gid_keys]
    if (length(mapped) < 5) next
    ego <- tryCatch(
      enrichGO(gene=mapped, OrgDb=org.G.max.eg.db, keyType="GID", ont="BP",
               pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff=0.2,
               readable=FALSE, minGSSize=5, maxGSSize=500),
      error=function(e) NULL)
    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
      df <- as.data.frame(ego) %>% filter(p.adjust <= go_padj_cutoff) %>%
        arrange(p.adjust) %>% head(go_max_terms) %>%
        mutate(neglog10p=-log10(p.adjust), Cluster=paste0("Cluster ",kk))
      if (nrow(df) > 0) go_df_list[[kk]] <- df
    }
  }
  go_all_df <- bind_rows(go_df_list)
  
  # Draw the heatmap.
  col_fun <- colorRamp2(seq(-2,2,length.out=101),
                        colorRampPalette(c("#313695","#4575b4","#74add1","#abd9e9","#e0f3f8",
                                           "#ffffbf","#fee090","#fdae61","#f46d43","#d73027","#a50026"))(101))
  
  top_ha <- HeatmapAnnotation(Pseudotime = anno_simple(1:n_bins, 
                                                       col=colorRamp2(c(1,n_bins/2,n_bins), c("#440154","#21918c","#fde725")), 
                                                       height=unit(4,"mm")), show_annotation_name=FALSE)
  
  cluster_colors <- c("#E41A1C","#377EB8","#4DAF4A")
  right_ha <- rowAnnotation(Module=anno_block(gp=gpar(fill=cluster_colors, col="white", lwd=2)))
  
  ht <- Heatmap(expr_scaled_ord, name="Z-score", col=col_fun,
                cluster_rows=FALSE, cluster_columns=FALSE,
                show_row_names=FALSE, show_column_names=FALSE,
                row_split=cluster_factor,
                column_title="Endodermis - Dynamic Gene Expression",
                top_annotation=top_ha, right_annotation=right_ha,
                width=unit(6,"cm"), height=unit(12,"cm"),
                use_raster=TRUE)
  
  # Draw the GO dot plot.
  dot_plots <- lapply(1:K, function(kk) {
    df <- go_all_df %>% filter(Cluster == paste0("Cluster ", kk))
    if (nrow(df) == 0) return(ggplot() + theme_void())
    df <- df %>% arrange(desc(neglog10p)) %>% head(go_max_terms) %>% arrange(neglog10p) %>%
      mutate(term = str_wrap(Description, width=40), term = factor(term, levels=unique(term)))
    
    ggplot(df, aes(x=neglog10p, y=term)) +
      geom_point(aes(size=Count, color=neglog10p)) +
      scale_color_gradient(low="#fddbc7", high="#b2182b") +
      labs(title=paste0("C",kk," GO(BP)"), x="-log10(p.adjust)", y=NULL) +
      theme_bw() + theme(axis.text.y=element_text(size=8))
  })
  go_panel <- wrap_plots(dot_plots, ncol=1)
  
  return(list(ht=ht, go_panel=go_panel, go_df=go_all_df))
}

# Helper: compute Spearman correlations.
fast_spearman_vs_pt <- function(expr_dense, pt_vals) {
  ranked <- t(apply(expr_dense, 1, rank))
  rpt_c <- rank(pt_vals) - mean(rank(pt_vals))
  ranked_c <- ranked - rowMeans(ranked)
  num <- as.numeric(ranked_c %*% rpt_c)
  den <- sqrt(rowSums(ranked_c^2) * sum(rpt_c^2))
  cors <- num / den
  cors[!is.finite(cors)] <- 0
  names(cors) <- rownames(expr_dense)
  cors
}

cat("Generating heatmaps and GO enrichment...\n")
res_heatmap <- run_heatmap_go(seurat_clean, n_genes = 500, n_bins = 50)

# Save the heatmap.
pdf("04_Heatmap_Dynamic_Genes.pdf", width = 9, height = 10)
draw(res_heatmap$ht)
dev.off()

# Save the GO dot plot.
ggsave("05_GO_Enrichment.pdf", res_heatmap$go_panel, width = 6, height = 8)

# ============================================================================
# Step 7: DEG bar plot (the source fragment ends at this heading).

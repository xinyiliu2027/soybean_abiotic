# Source: Figure 1 document, extracted text lines 547-744
# Continues the label-transfer session; retains the original diagnostic output.
# See README.md and docs/KNOWN_LIMITATIONS.md before running.

# Inspect and troubleshoot sample merging.
# Check the actual number of samples.
cat("Reference object name: reference\n")
cat("Hormone samples:", names(hormone.list), "\n")
cat("Total samples:", length(hormone.list), "\n")

# Construct sample names from REF and the current query list; check SA inclusion explicitly.
sample_names <- c("REF", names(hormone.list))
cat("Corrected sample names:", sample_names, "\n")
# Inspect metadata in the first sample.
cat("Checking metadata columns in first sample:\n")
print(colnames(hormone.list[[1]]@meta.data))

# Check the available grouping columns.
cat("\nFirst few rows of metadata:\n")
print(head(hormone.list[[1]]@meta.data))

# Look for prediction-related columns.
predicted_cols <- grep("predict", colnames(hormone.list[[1]]@meta.data), value = TRUE, ignore.case = TRUE)
cat("\nColumns containing 'predict':\n")
print(predicted_cols)

# Look for label-related columns.
label_cols <- grep("label", colnames(hormone.list[[1]]@meta.data), value = TRUE, ignore.case = TRUE)
cat("\nColumns containing 'label':\n")
print(label_cols)
library(Seurat)
library(ggplot2)
library(patchwork)

# 1. Create a DimPlot for each sample, including REF.
sample_names <- c("REF", "ABA1", "ABA2", "eBL1", "eBL2", "highTemp1", "highTemp2", 
                  "IAA1", "IAA2", "lowTemp1", "lowTemp2", "MEJA1", "MEJA2", 
                  "PEG1", "PEG2", "SA1", "SA2")

# Create a sample list that includes REF.
all_samples <- c(list(reference), hormone.list)
names(all_samples) <- sample_names

# Group cells using the cell_type column.
group_by_col <- "cell_type"
cat("Using", group_by_col, "for grouping\n")

# Create a list of DimPlots.
dimplot_list <- list()

for(i in 1:length(all_samples)) {
  sample_name <- sample_names[i]
  
  # Check whether this sample has a UMAP embedding.
  if(!"umap" %in% names(all_samples[[i]]@reductions)) {
    cat("Warning: No UMAP found for", sample_name, ", skipping...\n")
    next
  }
  
  # Check whether the cell_type column exists.
  if(!group_by_col %in% colnames(all_samples[[i]]@meta.data)) {
    cat("Warning: No", group_by_col, "found for", sample_name, ", using seurat_clusters\n")
    current_group_col <- "seurat_clusters"
  } else {
    current_group_col <- group_by_col
  }
  
  # Create a DimPlot.
  p <- DimPlot(all_samples[[i]], 
               reduction = "umap", 
               group.by = current_group_col,
               label = TRUE, 
               label.size = 3,
               pt.size = 0.5) +
    ggtitle(sample_name) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 14))
  
  dimplot_list[[i]] <- p
  
  cat("Created dimplot for", sample_name, "\n")
}

# Remove NULL entries, if present.
dimplot_list <- dimplot_list[!sapply(dimplot_list, is.null)]

# 2. Combine all DimPlots.
if(length(dimplot_list) > 0) {
  # Choose the number of rows and columns.
  n_plots <- length(dimplot_list)
  n_cols <- 4
  n_rows <- ceiling(n_plots / n_cols)
  
  combined_plot <- wrap_plots(dimplot_list, ncol = n_cols)
  
  # Save the combined plot.
  ggsave("all_samples_with_ref_dimplot.pdf", combined_plot, 
         width = 5 * n_cols, height = 4 * n_rows)
  ggsave("all_samples_with_ref_dimplot.png", combined_plot, 
         width = 5 * n_cols, height = 4 * n_rows, dpi = 300)
  
  cat("Saved combined dimplot with", n_plots, "samples (including REF)\n")
  
  # Display the plot.
  print(combined_plot)
} else {
  cat("No valid dimplots created\n")
}

# 3. Merge all samples, including REF.
cat("\nMerging all samples including REF...\n")

# Add sample identifiers to each sample.
for(i in 1:length(all_samples)) {
  all_samples[[i]]$sample_id <- sample_names[i]
  all_samples[[i]]$orig.ident <- sample_names[i]
  
  # Ensure that samples have consistent metadata columns.
  if(!"cell_type" %in% colnames(all_samples[[i]]@meta.data) && 
     sample_names[i] == "REF") {
    # If REF lacks cell_type, check another annotation column or cluster labels.
    if("predicted.celltype" %in% colnames(all_samples[[i]]@meta.data)) {
      all_samples[[i]]$cell_type <- all_samples[[i]]$predicted.celltype
    } else if("annotation" %in% colnames(all_samples[[i]]@meta.data)) {
      all_samples[[i]]$cell_type <- all_samples[[i]]$annotation
    } else {
      all_samples[[i]]$cell_type <- paste0("Cluster_", all_samples[[i]]$seurat_clusters)
    }
  }
}

# Merge all samples.
merged_seurat_with_ref <- merge(x = all_samples[[1]], 
                                y = all_samples[2:length(all_samples)], 
                                add.cell.ids = sample_names,
                                project = "HormoneStress_with_REF")

cat("Successfully merged", length(all_samples), "samples (including REF)\n")
cat("Total cells in merged object:", ncol(merged_seurat_with_ref), "\n")

# 4. Recompute dimensionality reduction for the merged object; assess whether this is appropriate.
cat("\nProcessing merged object with REF...\n")

# Normalize and scale expression.
merged_seurat_with_ref <- NormalizeData(merged_seurat_with_ref, verbose = FALSE)
merged_seurat_with_ref <- FindVariableFeatures(merged_seurat_with_ref, verbose = FALSE)
merged_seurat_with_ref <- ScaleData(merged_seurat_with_ref, verbose = FALSE)

# Reduce dimensionality.
merged_seurat_with_ref <- RunPCA(merged_seurat_with_ref, verbose = FALSE)
merged_seurat_with_ref <- RunUMAP(merged_seurat_with_ref, dims = 1:30, verbose = FALSE)

# 5. Create DimPlots of the merged object.
# Group by sample.
p1 <- DimPlot(merged_seurat_with_ref, 
              reduction = "umap", 
              group.by = "orig.ident",
              label = TRUE,
              pt.size = 0.3) +
  ggtitle("All Samples (including REF) - Grouped by Sample") +
  theme(plot.title = element_text(hjust = 0.5, size = 14))

# Group by cell type.
p2 <- DimPlot(merged_seurat_with_ref, 
              reduction = "umap", 
              group.by = "cell_type",
              label = TRUE,
              pt.size = 0.3) +
  ggtitle("All Samples (including REF) - Grouped by Cell Type") +
  theme(plot.title = element_text(hjust = 0.5, size = 14))

# Highlight REF versus treatment samples.
merged_seurat_with_ref$sample_type <- ifelse(merged_seurat_with_ref$orig.ident == "REF", 
                                             "Reference", "Treatment")

p3 <- DimPlot(merged_seurat_with_ref, 
              reduction = "umap", 
              group.by = "sample_type",
              cols = c("Reference" = "red", "Treatment" = "lightblue"),
              pt.size = 0.3) +
  ggtitle("Reference vs Treatment Samples") +
  theme(plot.title = element_text(hjust = 0.5, size = 14))

# Combine plots.
merged_plot_with_ref <- (p1 | p3) / p2

# Save the merged-sample plot.
ggsave("merged_samples_with_ref_dimplot.pdf", merged_plot_with_ref, 
       width = 16, height = 12)
ggsave("merged_samples_with_ref_dimplot.png", merged_plot_with_ref, 
       width = 16, height = 12, dpi = 300)

# Display the plot.
print(merged_plot_with_ref)


# Save the merged Seurat object.
saveRDS(merged_seurat_with_ref, file = "merged_seurat_with_ref.rds")

# Source: Figure 1 document, extracted text lines 1148-1321
# Continues the preceding session with combined plots and coordinate repair.
# See README.md and docs/KNOWN_LIMITATIONS.md before running.

# Revised visualization and nine-panel layout.
# Load the required packages.
library(Seurat)
library(ggplot2)
library(gridExtra)
library(dplyr)

# Define the color palette.
custom_colors <- c(
  "Cortex" = "#a2d2e7",           # Blue.
  "Dividing Cell" = "#67a8cd",    # Light blue.
  "Endodermis" = "#ffc17f",       # Yellow.
  "Epidermis" = "#cf9f88",        # Orange.
  "Lateral Root" = "#6fb3a8",     # Green.
  "Meristem" = "#b3e19b",         # Light green.
  "Pericycle" = "#50aa4b",        # Dark green.
  "Phloem" = "#ff9d9f",           # Pink.
  "QC" = "#f36569",               # Red.
  "Root Cap" = "#3581b7",         # Blue.
  "Xylem" = "#cdb6da"             # Light purple.
)

# Define the grouping column, usually the cell-type annotation.
group_by_col <- "cell_type"  # Alternatively use seurat_clusters if cell_type is unavailable.

# Collect all samples into a lookup list.
all_samples <- c(list("REF" = reference), hormone.list)

# Define the sample pairs to merge.
merge_pairs <- list(
  "ABA" = c("ABA1", "ABA2"),
  "eBL" = c("eBL1", "eBL2"),
  "highTemp" = c("highTemp1", "highTemp2"),
  "IAA" = c("IAA1", "IAA2"),
  "lowTemp" = c("lowTemp1", "lowTemp2"),
  "MEJA" = c("MEJA1", "MEJA2"),
  "PEG" = c("PEG1", "PEG2"),
  "SA" = c("SA1", "SA2")
)

# Repair existing merged samples by adding UMAP coordinates.
cat("Fixing UMAP information for merged samples...\n")

for(merged_name in names(merge_pairs)) {
  if(merged_name %in% names(merged_samples)) {
    sample_pair <- merge_pairs[[merged_name]]
    
    # Retrieve the original samples.
    sample1 <- all_samples[[sample_pair[1]]]
    sample2 <- all_samples[[sample_pair[2]]]
    
    if("umap" %in% names(sample1@reductions) && "umap" %in% names(sample2@reductions)) {
      # Retrieve the merged sample.
      merged_sample <- merged_samples[[merged_name]]
      
      # Combine UMAP coordinates manually.
      umap1 <- sample1@reductions$umap@cell.embeddings
      umap2 <- sample2@reductions$umap@cell.embeddings
      
      # Prefix cell IDs in the UMAP coordinates.
      rownames(umap1) <- paste0(sample_pair[1], "_", rownames(umap1))
      rownames(umap2) <- paste0(sample_pair[2], "_", rownames(umap2))
      
      # Combine UMAP coordinates.
      combined_umap <- rbind(umap1, umap2)
      
      # Match the UMAP coordinate order to cells in the merged sample.
      if(all(colnames(merged_sample) %in% rownames(combined_umap))) {
        combined_umap <- combined_umap[colnames(merged_sample), ]
        
        # Create a new UMAP reduction object.
        merged_sample[["umap"]] <- CreateDimReducObject(
          embeddings = combined_umap,
          key = "UMAP_",
          assay = DefaultAssay(merged_sample)
        )
        
        # Update the merged sample.
        merged_samples[[merged_name]] <- merged_sample
        
        cat("Fixed UMAP for", merged_name, "\n")
      } else {
        cat("Warning: Cell names don't match for", merged_name, "\n")
      }
    }
  }
}

# Create DimPlots.
merged_dimplot_list <- list()

for(sample_name in names(merged_samples)) {
  sample_obj <- merged_samples[[sample_name]]
  
  # Check whether this sample has a UMAP embedding.
  if(!"umap" %in% names(sample_obj@reductions)) {
    cat("Warning: No UMAP found for", sample_name, ", skipping...\n")
    next
  }
  
  # Check whether the cell_type column exists.
  if(!group_by_col %in% colnames(sample_obj@meta.data)) {
    cat("Warning: No", group_by_col, "found for", sample_name, ", using seurat_clusters\n")
    current_group_col <- "seurat_clusters"
  } else {
    current_group_col <- group_by_col
  }
  
  # Create a DimPlot.
  p <- DimPlot(sample_obj, 
               reduction = "umap", 
               group.by = current_group_col,
               label = FALSE,  # Disable labels to reduce crowding.
               pt.size = 0.3,
               cols = custom_colors) +  # Use custom colors.
    ggtitle(sample_name) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 8))
  
  merged_dimplot_list[[sample_name]] <- p
  
  cat("Created dimplot for", sample_name, "\n")
}

# Arrange and display plots.
if(length(merged_dimplot_list) > 0) {
  # Arrange panels in the specified order.
  plot_order <- c("REF", "ABA", "eBL", "highTemp", 
                  "IAA", "lowTemp", "MEJA", "PEG", "SA")
  
  # Reorder panels, retaining only those that exist.
  ordered_plots <- list()
  for(name in plot_order) {
    if(name %in% names(merged_dimplot_list)) {
      ordered_plots[[length(ordered_plots) + 1]] <- merged_dimplot_list[[name]]
    }
  }
  
  # Create a 3 x 3 layout.
  if(length(ordered_plots) > 0) {
    grid_plot <- do.call(grid.arrange, c(ordered_plots, ncol = 3))
    
    # Save the plots.
    ggsave("merged_samples_dimplot.pdf", grid_plot, width = 12, height = 12, dpi = 300)
    ggsave("merged_samples_dimplot.png", grid_plot, width = 12, height = 12, dpi = 300)
    
    print(grid_plot)
    cat("Plots saved as merged_samples_dimplot.pdf and .png\n")
  } else {
    cat("No plots to display\n")
  }
} else {
  cat("No valid plots created\n")
}

# Print the merge summary.
cat("\n=== Merge summary ===\n")
cat("Successfully merged samples:", length(merged_samples), "\n")
for(name in names(merged_samples)) {
  n_cells <- ncol(merged_samples[[name]])
  has_umap <- "umap" %in% names(merged_samples[[name]]@reductions)
  cat(name, ":", n_cells, "cells, UMAP:", has_umap, "\n")
}
cat("Successfully transferred UMAP coordinates to merged object\n")
# Save an RData file.
save(merged_samples, file = "merged_samples.RData")
cat("Merged samples saved as merged_samples.RData.\n")

# Alternatively save individual RDS objects; this may yield smaller files.
saveRDS(merged_samples, file = "merged_samples.rds")
cat("Merged samples saved as merged_samples.rds.\n")

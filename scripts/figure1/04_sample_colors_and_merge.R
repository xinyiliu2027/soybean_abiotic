# Source: Figure 1 document, extracted text lines 745-926
# Continues the preceding session and creates merged_samples.
# See README.md and docs/KNOWN_LIMITATIONS.md before running.

# Customize colors in the split UMAP plots.
library(Seurat)
library(ggplot2)
library(patchwork)

# Custom color configuration.
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

# 1. Create a custom-colored DimPlot for each sample, including REF.
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
  
  # Create a DimPlot with custom colors.
  p <- DimPlot(all_samples[[i]], 
               reduction = "umap", 
               group.by = current_group_col,
               label = TRUE, 
               label.size = 3,
               pt.size = 0.5,
               cols = custom_colors) +  # Apply custom colors.
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
  ggsave("all_samples_with_ref_dimplot_custom_colors.pdf", combined_plot, 
         width = 5 * n_cols, height = 4 * n_rows)
  ggsave("all_samples_with_ref_dimplot_custom_colors.png", combined_plot, 
         width = 5 * n_cols, height = 4 * n_rows, dpi = 300)
  
  cat("Saved individual sample plots with custom colors\n")
}

# 3. Merge samples.
# Revised sample-merging code.
# Define the color palette.
if (exists("umap_plots", inherits = FALSE)) rm(umap_plots)

# Optional faster alternative that skips reclustering.
# This variant performs only the steps needed to generate UMAP plots.

# Optional simplified merging code.
merge_samples_fast <- function() {
  merged_samples_fast <- list()
  
  for(i in 1:length(merge_pairs)) {
    merged_name <- names(merge_pairs)[i]
    sample_pair <- merge_pairs[[i]]
    
    if(all(sample_pair %in% names(all_samples))) {
      sample1 <- all_samples[[sample_pair[1]]]
      sample2 <- all_samples[[sample_pair[2]]]
      
      if("umap" %in% names(sample1@reductions) && "umap" %in% names(sample2@reductions)) {
        # Merge samples.
        merged_sample <- merge(sample1, sample2, 
                               add.cell.ids = sample_pair,
                               project = merged_name)
        
        cat("Processing", merged_name, "with fast workflow...\n")
        
        # Run only the required steps.
        merged_sample <- NormalizeData(merged_sample, verbose = FALSE)
        merged_sample <- FindVariableFeatures(merged_sample, verbose = FALSE)
        merged_sample <- ScaleData(merged_sample, verbose = FALSE)
        merged_sample <- RunPCA(merged_sample, verbose = FALSE)
        merged_sample <- RunUMAP(merged_sample, dims = 1:30, verbose = FALSE)
        
        merged_samples_fast[[merged_name]] <- merged_sample
        cat("Successfully processed", merged_name, "\n")
      }
    }
  }
  
  return(merged_samples_fast)
}

# To use the faster variant, run the following command.
# merged_samples <- merge_samples_fast()
# Define sample pairs to merge.
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

# Create the merged sample objects.
merged_samples <- list()
merged_samples[["REF"]] <- all_samples[["REF"]]  # Add the reference sample.

for(merged_name in names(merge_pairs)) {
  sample_pair <- merge_pairs[[merged_name]]
  
  # Check that both samples exist.
  if(all(sample_pair %in% names(all_samples))) {
    # Retrieve both samples.
    sample1 <- all_samples[[sample_pair[1]]]
    sample2 <- all_samples[[sample_pair[2]]]
    
    # Check that both samples have UMAP embeddings.
    if("umap" %in% names(sample1@reductions) && "umap" %in% names(sample2@reductions)) {
      # Merge samples.
      merged_sample <- merge(sample1, sample2, 
                             add.cell.ids = sample_pair,
                             project = merged_name)
      
      # Uncomment the following line if UMAP should be recomputed after merging.
      # merged_sample <- RunUMAP(merged_sample, dims = 1:30)
      
      # Store the merged sample.
      merged_samples[[merged_name]] <- merged_sample
      
      cat("Successfully merged", paste(sample_pair, collapse = " and "), "into", merged_name, "\n")
    } else {
      cat("Warning: UMAP not found for", paste(sample_pair, collapse = " or "), ", skipping merge\n")
    }
  } else {
    cat("Warning: Not all samples found for", merged_name, ", skipping merge\n")
  }
}

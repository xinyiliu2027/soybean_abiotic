# Source: Figure 1 document, extracted text lines 281-546
# Requires preloaded Gm_* samples and reference_renamed.qs.
# See README.md and docs/KNOWN_LIMITATIONS.md before running.

# Anchor-based label transfer.

# Load qs.
library(qs)

# Read reference_renamed.qs.
reference_renamed <- qread(Sys.getenv("SOYBEAN_REFERENCE_QS", unset = "data/reference_renamed.qs"))
# Inspect objects whose names contain Gm.
gm_objects <- ls(pattern = "Gm")
print("Objects containing 'Gm':")
print(gm_objects)

# Define the object-renaming map.
rename_map <- c(
  "Gm_1031_N327F" = "eBL1",
  "Gm_2213_eBL" = "eBL2", 
  "Gm_2224_8.37temp" = "highTemp1",
  "Gm_2225_9.6Temp" = "highTemp2",
  "Gm_2226_9.6Temp" = "lowTemp1",
  "Gm_2227_9.6Temp" = "lowTemp2",
  "Gm_2214_SA" = "SA1",
  "Gm_2215_SA" = "SA2",
  "Gm_2216_MEJA" = "MEJA1",
  "Gm_2217_MEJA" = "MEJA2",
  "Gm_2218_ABA" = "ABA1",
  "Gm_2219_ABA" = "ABA2",
  "Gm_2222_PEG" = "PEG1",
  "Gm_2223_PEG" = "PEG2",
  "Gm_2220_IAA" = "IAA1",
  "Gm_2221_IAA" = "IAA2"
)

# Process all Gm objects.
for (old_name in gm_objects) {
  # Check whether the CK object should be removed from the session.
  if (old_name %in% c("Gm_2228_CK", "Gm_2229_CK")) {
    rm(list = old_name)
    cat("Removed:", old_name, "\n")
  } 
  # Check whether the object is in the renaming map.
  else if (old_name %in% names(rename_map)) {
    new_name <- rename_map[old_name]
    assign(new_name, get(old_name))
    rm(list = old_name)
    cat("Renamed:", old_name, "->", new_name, "\n")
  }
}

# Inspect the processed objects.
print("\nObjects after processing:")
seurat_objects <- ls()[sapply(ls(), function(x) class(get(x))[1] == "Seurat")]
print(seurat_objects)


# Find all Seurat objects.
seurat_objects <- ls()[sapply(ls(), function(x) class(get(x))[1] == "Seurat")]
print("Current Seurat objects:")
print(seurat_objects)
# Inspect metadata columns for doublet annotations.
for (obj_name in seurat_objects) {
  obj <- get(obj_name)
  
  cat("\n=== Checking", obj_name, "===\n")
  
  # Retrieve all column names.
  col_names <- colnames(obj@meta.data)
  
  # Find columns that may contain doublet annotations.
  doublet_cols <- grep("doublet|DF|DoubletFinder", col_names, ignore.case = TRUE, value = TRUE)
  classification_cols <- grep("classification", col_names, ignore.case = TRUE, value = TRUE)
  singlet_cols <- grep("singlet", col_names, ignore.case = TRUE, value = TRUE)
  
  if (length(doublet_cols) > 0) {
    cat("Found doublet-related columns:", paste(doublet_cols, collapse = ", "), "\n")
  }
  if (length(classification_cols) > 0) {
    cat("Found classification columns:", paste(classification_cols, collapse = ", "), "\n")
  }
  if (length(singlet_cols) > 0) {
    cat("Found singlet-related columns:", paste(singlet_cols, collapse = ", "), "\n")
  }
  
  # Display all column names when the number of columns is small.
  if (length(col_names) <= 20) {
    cat("All meta.data columns:", paste(col_names, collapse = ", "), "\n")
  } else {
    cat("Total columns:", length(col_names), "\n")
    cat("First 10 columns:", paste(head(col_names, 10), collapse = ", "), "\n")
    cat("Last 10 columns:", paste(tail(col_names, 10), collapse = ", "), "\n")
  }
}



# Remove doublets from all samples.

# Remove doublets from all samples using the revised subsetting syntax.
for (obj_name in seurat_objects) {
  obj <- get(obj_name)
  
  cat("\n=== Processing", obj_name, "===\n")
  cat("Original cells:", ncol(obj), "\n")
  
  # Find the DF.classifications column; its name varies between objects.
  df_col <- grep("^DF\\.classifications", colnames(obj@meta.data), value = TRUE)
  
  if (length(df_col) > 1L) stop("Multiple DoubletFinder classification columns in ", obj_name, "; select one explicitly.")
  if (length(df_col) > 0) {
    cat("Found doublet column:", df_col, "\n")
    print(table(obj@meta.data[[df_col]]))
    
    # Retain only cells classified as Singlet.
    cells_to_keep <- rownames(obj@meta.data)[obj@meta.data[[df_col]] == "Singlet"]
    obj_filtered <- obj[, cells_to_keep]
    
    cat("After removing doublets:", ncol(obj_filtered), "cells\n")
    cat("Removed", ncol(obj) - ncol(obj_filtered), "doublets\n")
    
    # Update the object.
    assign(obj_name, obj_filtered)
    
  } else {
    cat("No doublet information found\n")
  }
}

# Display the processed results.
cat("\n=== Summary after doublet removal ===\n")
for (obj_name in seurat_objects) {
  obj <- get(obj_name)
  cat(obj_name, ":", ncol(obj), "cells\n")
}

# Clean temporary session objects.
rm(list = intersect(c("obj", "obj_filtered", "df_col", "cells_to_keep"), ls()))
rm(reference_renamed)

# Load the required packages.
library(Seurat)
library(dplyr)
library(patchwork)

# Load the reference data.
reference_renamed <- qread(Sys.getenv("SOYBEAN_REFERENCE_QS", unset = "data/reference_renamed.qs"))

# Collect the processed samples in a list.
hormone.list <- list()
for (obj_name in seurat_objects) {
  hormone.list[[obj_name]] <- get(obj_name)
}

# Inspect the reference data.
cat("Reference object info:\n")
cat("Cells:", ncol(reference_renamed), "\n")
cat("Features:", nrow(reference_renamed), "\n")
cat("Cell types:", length(unique(reference_renamed$cell_type)), "\n")

# Apply standard preprocessing to each query dataset.
cat("\n=== Preprocessing query datasets ===\n")
for (i in 1:length(hormone.list)) {
  cat("Processing", names(hormone.list)[i], "\n")
  hormone.list[[i]] <- NormalizeData(hormone.list[[i]], verbose = FALSE) %>% 
    FindVariableFeatures(verbose = FALSE) %>% 
    SCTransform(variable.features.n = 3000,
                method = "glmGamPoi",
                verbose = FALSE)
}

# Prepare the reference data.
cat("\n=== Preparing reference dataset ===\n")
DefaultAssay(reference_renamed) <- "RNA"
reference_renamed <- JoinLayers(reference_renamed)

reference <- reference_renamed

reference <- SCTransform(reference,
                         variable.features.n = 3000,
                         method = "glmGamPoi",
                         verbose = FALSE)

# Run reference PCA and UMAP, retaining the UMAP model for projection.
reference <- RunPCA(reference, verbose = FALSE)
reference <- RunUMAP(reference, reduction = "pca", dims = 1:50, 
                     umap.method = "umap-learn", return.model = TRUE,
                     verbose = FALSE)

DefaultAssay(reference) <- "RNA"

cat("Reference preparation completed!\n")


# Removed session-only rm(SA1); the sample list already retains its own objects.
cat("\n=== Performing anchor integration and label transfer ===\n")

# Transfer labels to each query dataset.
# Transfer labels to each query dataset.
for (i in names(hormone.list)) {
  cat("\n--- Processing", i, "---\n")
  
  # Extract the query dataset.
  query <- hormone.list[[i]]
  
  # Set SCT as the default query assay.
  DefaultAssay(query) <- "SCT"
  DefaultAssay(reference) <- "SCT"
  
  # Find anchors between the reference and query using FindTransferAnchors.
  # Use the SCT assay without recomputing the SCT model.
  anchors <- FindTransferAnchors(
    reference = reference,
    query = query,
    reference.assay = "SCT",
    query.assay = "SCT",
    normalization.method = "SCT",
    recompute.residuals = FALSE,
    dims = 1:50,
    reference.reduction = "pca"
  )
  
  # Transfer annotations using TransferData.
  predictions <- TransferData(
    anchorset = anchors,
    refdata = reference$cell_type,  # Reference cell-type labels to transfer.
    dims = 1:50,
    weight.reduction = "pcaproject"
  )
  
  # Add the transferred annotations to query metadata.
  query <- AddMetaData(query, metadata = predictions$predicted.id, col.name = "cell_type")
  
  # UMAP transfer
  map_query <- MapQuery(
    anchorset = anchors, 
    reference = reference, 
    query = query, 
    refdata = list(celltype = "cell_type"), 
    reference.reduction = "pca", 
    reduction.model = "umap"
  )
  
  query@reductions$umap <- map_query@reductions$ref.umap
  
  # Update the Seurat object in the sample list.
  hormone.list[[i]] <- query
  
  cat("Completed annotation for", i, "\n")
}

cat("\nLabel transfer completed for all samples!\n")

# Define a helper for merging the reference and query datasets.
clean_with_joinlayers <- function() {
  cat("\n=== Alternative: Using JoinLayers ===\n")
  
  # Clean each object
  reference_joined <- JoinLayers(reference)
  hormone_joined <- lapply(hormone.list, JoinLayers)
  names(hormone_joined) <- names(hormone.list)
  
  # Try merge
  merged_joined <- merge(x = reference_joined, 
                         y = hormone_joined, 
                         add.cell.ids = c("ref", names(hormone_joined)))
  
  return(merged_joined)
}

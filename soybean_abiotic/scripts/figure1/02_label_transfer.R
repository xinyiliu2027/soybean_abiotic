# Source: 图1.docx, extracted text lines 281-546
# 需要预载 Gm_* 样本和 reference_renamed.qs
# See README.md and docs/KNOWN_LIMITATIONS.md before running.

##转移锚点

# 加载 qs 包
library(qs)

# 读取 reference_renamed.qs 文件
reference_renamed <- qread(Sys.getenv("SOYBEAN_REFERENCE_QS", unset = "data/reference_renamed.qs"))
# 查看包含Gm的对象
gm_objects <- ls(pattern = "Gm")
print("Objects containing 'Gm':")
print(gm_objects)

# 创建重命名映射（移除了最后的多余逗号）
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

# 处理所有Gm对象
for (old_name in gm_objects) {
  # 检查是否需要删除CK对象
  if (old_name %in% c("Gm_2228_CK", "Gm_2229_CK")) {
    rm(list = old_name)
    cat("Removed:", old_name, "\n")
  } 
  # 检查是否在重命名映射中
  else if (old_name %in% names(rename_map)) {
    new_name <- rename_map[old_name]
    assign(new_name, get(old_name))
    rm(list = old_name)
    cat("Renamed:", old_name, "->", new_name, "\n")
  }
}

# 查看处理后的对象
print("\nObjects after processing:")
seurat_objects <- ls()[sapply(ls(), function(x) class(get(x))[1] == "Seurat")]
print(seurat_objects)


# 获取所有Seurat对象
seurat_objects <- ls()[sapply(ls(), function(x) class(get(x))[1] == "Seurat")]
print("Current Seurat objects:")
print(seurat_objects)
# 检查所有对象的meta.data列名，寻找可能的双细胞信息
for (obj_name in seurat_objects) {
  obj <- get(obj_name)
  
  cat("\n=== Checking", obj_name, "===\n")
  
  # 获取所有列名
  col_names <- colnames(obj@meta.data)
  
  # 查找可能包含双细胞信息的列
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
  
  # 显示所有列名（如果数量不多的话）
  if (length(col_names) <= 20) {
    cat("All meta.data columns:", paste(col_names, collapse = ", "), "\n")
  } else {
    cat("Total columns:", length(col_names), "\n")
    cat("First 10 columns:", paste(head(col_names, 10), collapse = ", "), "\n")
    cat("Last 10 columns:", paste(tail(col_names, 10), collapse = ", "), "\n")
  }
}



# 移除所有样本的双细胞

# 移除所有样本的双细胞 - 修正版
for (obj_name in seurat_objects) {
  obj <- get(obj_name)
  
  cat("\n=== Processing", obj_name, "===\n")
  cat("Original cells:", ncol(obj), "\n")
  
  # 查找DF.classifications列（每个对象的列名不同）
  df_col <- grep("^DF\\.classifications", colnames(obj@meta.data), value = TRUE)
  
  if (length(df_col) > 1L) stop("Multiple DoubletFinder classification columns in ", obj_name, "; select one explicitly.")
  if (length(df_col) > 0) {
    cat("Found doublet column:", df_col, "\n")
    print(table(obj@meta.data[[df_col]]))
    
    # 移除双细胞，只保留Singlet - 修正的语法
    cells_to_keep <- rownames(obj@meta.data)[obj@meta.data[[df_col]] == "Singlet"]
    obj_filtered <- obj[, cells_to_keep]
    
    cat("After removing doublets:", ncol(obj_filtered), "cells\n")
    cat("Removed", ncol(obj) - ncol(obj_filtered), "doublets\n")
    
    # 更新对象
    assign(obj_name, obj_filtered)
    
  } else {
    cat("No doublet information found\n")
  }
}

# 显示处理后的结果
cat("\n=== Summary after doublet removal ===\n")
for (obj_name in seurat_objects) {
  obj <- get(obj_name)
  cat(obj_name, ":", ncol(obj), "cells\n")
}

# 清理环境
rm(list = intersect(c("obj", "obj_filtered", "df_col", "cells_to_keep"), ls()))
rm(reference_renamed)

# 加载必要的包
library(Seurat)
library(dplyr)
library(patchwork)

# 加载参考数据
reference_renamed <- qread(Sys.getenv("SOYBEAN_REFERENCE_QS", unset = "data/reference_renamed.qs"))

# 将所有处理后的样本放入列表
hormone.list <- list()
for (obj_name in seurat_objects) {
  hormone.list[[obj_name]] <- get(obj_name)
}

# 检查参考数据
cat("Reference object info:\n")
cat("Cells:", ncol(reference_renamed), "\n")
cat("Features:", nrow(reference_renamed), "\n")
cat("Cell types:", length(unique(reference_renamed$cell_type)), "\n")

# 对每个查询数据集进行标准预处理
cat("\n=== Preprocessing query datasets ===\n")
for (i in 1:length(hormone.list)) {
  cat("Processing", names(hormone.list)[i], "\n")
  hormone.list[[i]] <- NormalizeData(hormone.list[[i]], verbose = FALSE) %>% 
    FindVariableFeatures(verbose = FALSE) %>% 
    SCTransform(variable.features.n = 3000,
                method = "glmGamPoi",
                verbose = FALSE)
}

# 准备参考数据
cat("\n=== Preparing reference dataset ===\n")
DefaultAssay(reference_renamed) <- "RNA"
reference_renamed <- JoinLayers(reference_renamed)

reference <- reference_renamed

reference <- SCTransform(reference,
                         variable.features.n = 3000,
                         method = "glmGamPoi",
                         verbose = FALSE)

# 为参考数据运行PCA和UMAP（返回模型用于映射）
reference <- RunPCA(reference, verbose = FALSE)
reference <- RunUMAP(reference, reduction = "pca", dims = 1:50, 
                     umap.method = "umap-learn", return.model = TRUE,
                     verbose = FALSE)

DefaultAssay(reference) <- "RNA"

cat("Reference preparation completed!\n")


# Removed session-only rm(SA1); the sample list already retains its own objects.
cat("\n=== Performing anchor integration and label transfer ===\n")

# 对每个查询数据集进行标签转移
# 对每个查询数据集进行标签转移
for (i in names(hormone.list)) {
  cat("\n--- Processing", i, "---\n")
  
  # 提取查询数据集
  query <- hormone.list[[i]]
  
  # 确保查询数据集的默认assay是SCT
  DefaultAssay(query) <- "SCT"
  DefaultAssay(reference) <- "SCT"
  
  # 使用 FindTransferAnchors 找到参考和查询数据集之间的锚点
  # 使用SCT assay但避免重新计算SCT模型
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
  
  # 使用 TransferData 迁移注释
  predictions <- TransferData(
    anchorset = anchors,
    refdata = reference$cell_type,  # 迁移参考数据集的细胞类型
    dims = 1:50,
    weight.reduction = "pcaproject"
  )
  
  # 将迁移的注释添加到查询数据集的 metadata 中
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
  
  # 更新列表中的 Seurat 对象
  hormone.list[[i]] <- query
  
  cat("Completed annotation for", i, "\n")
}

cat("\nLabel transfer completed for all samples!\n")

# 首先合并参考数据集和第一个查询数据集
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

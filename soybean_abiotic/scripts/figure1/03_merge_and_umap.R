# Source: 图1.docx, extracted text lines 547-744
# 续接标签转移会话；含原始诊断输出
# See README.md and docs/KNOWN_LIMITATIONS.md before running.

##merge仍然有问题进行尝试
# 检查实际的样本数量
cat("Reference object name: reference\n")
cat("Hormone samples:", names(hormone.list), "\n")
cat("Total samples:", length(hormone.list), "\n")

# 正确的样本名称（去掉SA1, SA2）
sample_names <- c("REF", names(hormone.list))
cat("Corrected sample names:", sample_names, "\n")
# 首先检查第一个样本的metadata结构
cat("Checking metadata columns in first sample:\n")
print(colnames(hormone.list[[1]]@meta.data))

# 检查可用的分组列
cat("\nFirst few rows of metadata:\n")
print(head(hormone.list[[1]]@meta.data))

# 检查是否有预测相关的列
predicted_cols <- grep("predict", colnames(hormone.list[[1]]@meta.data), value = TRUE, ignore.case = TRUE)
cat("\nColumns containing 'predict':\n")
print(predicted_cols)

# 检查是否有标签相关的列
label_cols <- grep("label", colnames(hormone.list[[1]]@meta.data), value = TRUE, ignore.case = TRUE)
cat("\nColumns containing 'label':\n")
print(label_cols)
library(Seurat)
library(ggplot2)
library(patchwork)

# 1. 为每个样本创建dimplot（包含REF）
sample_names <- c("REF", "ABA1", "ABA2", "eBL1", "eBL2", "highTemp1", "highTemp2", 
                  "IAA1", "IAA2", "lowTemp1", "lowTemp2", "MEJA1", "MEJA2", 
                  "PEG1", "PEG2", "SA1", "SA2")

# 创建包含REF的样本列表
all_samples <- c(list(reference), hormone.list)
names(all_samples) <- sample_names

# 使用 cell_type 列进行分组
group_by_col <- "cell_type"
cat("Using", group_by_col, "for grouping\n")

# 创建dimplot列表
dimplot_list <- list()

for(i in 1:length(all_samples)) {
  sample_name <- sample_names[i]
  
  # 检查该样本是否有UMAP
  if(!"umap" %in% names(all_samples[[i]]@reductions)) {
    cat("Warning: No UMAP found for", sample_name, ", skipping...\n")
    next
  }
  
  # 检查是否有cell_type列
  if(!group_by_col %in% colnames(all_samples[[i]]@meta.data)) {
    cat("Warning: No", group_by_col, "found for", sample_name, ", using seurat_clusters\n")
    current_group_col <- "seurat_clusters"
  } else {
    current_group_col <- group_by_col
  }
  
  # 创建dimplot
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

# 移除NULL元素（如果有的话）
dimplot_list <- dimplot_list[!sapply(dimplot_list, is.null)]

# 2. 组合显示所有dimplot
if(length(dimplot_list) > 0) {
  # 计算合适的行列数
  n_plots <- length(dimplot_list)
  n_cols <- 4
  n_rows <- ceiling(n_plots / n_cols)
  
  combined_plot <- wrap_plots(dimplot_list, ncol = n_cols)
  
  # 保存组合图
  ggsave("all_samples_with_ref_dimplot.pdf", combined_plot, 
         width = 5 * n_cols, height = 4 * n_rows)
  ggsave("all_samples_with_ref_dimplot.png", combined_plot, 
         width = 5 * n_cols, height = 4 * n_rows, dpi = 300)
  
  cat("Saved combined dimplot with", n_plots, "samples (including REF)\n")
  
  # 显示图片
  print(combined_plot)
} else {
  cat("No valid dimplots created\n")
}

# 3. 合并所有样本（包含REF）
cat("\nMerging all samples including REF...\n")

# 为每个样本添加样本标识
for(i in 1:length(all_samples)) {
  all_samples[[i]]$sample_id <- sample_names[i]
  all_samples[[i]]$orig.ident <- sample_names[i]
  
  # 确保所有样本都有相同的metadata列
  if(!"cell_type" %in% colnames(all_samples[[i]]@meta.data) && 
     sample_names[i] == "REF") {
    # 如果REF没有cell_type列，可能需要从其他列获取或使用聚类信息
    if("predicted.celltype" %in% colnames(all_samples[[i]]@meta.data)) {
      all_samples[[i]]$cell_type <- all_samples[[i]]$predicted.celltype
    } else if("annotation" %in% colnames(all_samples[[i]]@meta.data)) {
      all_samples[[i]]$cell_type <- all_samples[[i]]$annotation
    } else {
      all_samples[[i]]$cell_type <- paste0("Cluster_", all_samples[[i]]$seurat_clusters)
    }
  }
}

# 合并所有样本
merged_seurat_with_ref <- merge(x = all_samples[[1]], 
                                y = all_samples[2:length(all_samples)], 
                                add.cell.ids = sample_names,
                                project = "HormoneStress_with_REF")

cat("Successfully merged", length(all_samples), "samples (including REF)\n")
cat("Total cells in merged object:", ncol(merged_seurat_with_ref), "\n")

# 4. 为合并后的对象重新计算降维(这一步需要吗)
cat("\nProcessing merged object with REF...\n")

# 标准化和缩放
merged_seurat_with_ref <- NormalizeData(merged_seurat_with_ref, verbose = FALSE)
merged_seurat_with_ref <- FindVariableFeatures(merged_seurat_with_ref, verbose = FALSE)
merged_seurat_with_ref <- ScaleData(merged_seurat_with_ref, verbose = FALSE)

# 降维
merged_seurat_with_ref <- RunPCA(merged_seurat_with_ref, verbose = FALSE)
merged_seurat_with_ref <- RunUMAP(merged_seurat_with_ref, dims = 1:30, verbose = FALSE)

# 5. 创建合并后的dimplot
# 按样本分组
p1 <- DimPlot(merged_seurat_with_ref, 
              reduction = "umap", 
              group.by = "orig.ident",
              label = TRUE,
              pt.size = 0.3) +
  ggtitle("All Samples (including REF) - Grouped by Sample") +
  theme(plot.title = element_text(hjust = 0.5, size = 14))

# 按细胞类型分组
p2 <- DimPlot(merged_seurat_with_ref, 
              reduction = "umap", 
              group.by = "cell_type",
              label = TRUE,
              pt.size = 0.3) +
  ggtitle("All Samples (including REF) - Grouped by Cell Type") +
  theme(plot.title = element_text(hjust = 0.5, size = 14))

# 创建一个突出显示REF vs 其他样本的图
merged_seurat_with_ref$sample_type <- ifelse(merged_seurat_with_ref$orig.ident == "REF", 
                                             "Reference", "Treatment")

p3 <- DimPlot(merged_seurat_with_ref, 
              reduction = "umap", 
              group.by = "sample_type",
              cols = c("Reference" = "red", "Treatment" = "lightblue"),
              pt.size = 0.3) +
  ggtitle("Reference vs Treatment Samples") +
  theme(plot.title = element_text(hjust = 0.5, size = 14))

# 组合图
merged_plot_with_ref <- (p1 | p3) / p2

# 保存合并后的图
ggsave("merged_samples_with_ref_dimplot.pdf", merged_plot_with_ref, 
       width = 16, height = 12)
ggsave("merged_samples_with_ref_dimplot.png", merged_plot_with_ref, 
       width = 16, height = 12, dpi = 300)

# 显示图
print(merged_plot_with_ref)


# 保存合并后的Seurat对象
saveRDS(merged_seurat_with_ref, file = "merged_seurat_with_ref.rds")

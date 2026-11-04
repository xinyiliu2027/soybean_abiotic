# Source: 图1.docx, extracted text lines 927-1147
# 续接上一步的分组绘图版本
# See README.md and docs/KNOWN_LIMITATIONS.md before running.

##修复4的errorrm
if (exists("all_samples", inherits = FALSE)) rm(all_samples)

library(Seurat)
library(ggplot2)
library(patchwork)

# 自定义颜色配置
custom_colors <- c(
  "Cortex" = "#a2d2e7",           # 蓝色
  "Dividing Cell" = "#67a8cd",    # 浅蓝色
  "Endodermis" = "#ffc17f",       # 黄色
  "Epidermis" = "#cf9f88",        # 橙色
  "Lateral Root" = "#6fb3a8",     # 绿色
  "Meristem" = "#b3e19b",         # 浅绿色
  "Pericycle" = "#50aa4b",        # 深绿色
  "Phloem" = "#ff9d9f",           # 粉色
  "QC" = "#f36569",               # 红色
  "Root Cap" = "#3581b7",         # 蓝色
  "Xylem" = "#cdb6da"             # 浅紫色
)

# 1. 为每个样本创建dimplot（包含REF）- 使用自定义颜色
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
  
  # 创建dimplot - 添加自定义颜色
  p <- DimPlot(all_samples[[i]], 
               reduction = "umap", 
               group.by = current_group_col,
               label = TRUE, 
               label.size = 3,
               pt.size = 0.5,
               cols = custom_colors) +  # 添加自定义颜色
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
  ggsave("all_samples_with_ref_dimplot_custom_colors.pdf", combined_plot, 
         width = 5 * n_cols, height = 4 * n_rows)
  ggsave("all_samples_with_ref_dimplot_custom_colors.png", combined_plot, 
         width = 5 * n_cols, height = 4 * n_rows, dpi = 300)
  
  cat("Saved individual sample plots with custom colors\n")
}

# 3. 合并样本功能
# 定义需要合并的样本对
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

# 创建合并后的样本对象
merged_samples <- list()
merged_samples[["REF"]] <- all_samples[["REF"]]  # 添加参考样本

for(merged_name in names(merge_pairs)) {
  sample_pair <- merge_pairs[[merged_name]]
  
  # 检查两个样本是否都存在
  if(all(sample_pair %in% names(all_samples))) {
    # 获取两个样本
    sample1 <- all_samples[[sample_pair[1]]]
    sample2 <- all_samples[[sample_pair[2]]]
    
    # 检查两个样本是否都有UMAP
    if("umap" %in% names(sample1@reductions) && "umap" %in% names(sample2@reductions)) {
      # 合并样本
      merged_sample <- merge(sample1, sample2, 
                             add.cell.ids = sample_pair,
                             project = merged_name)
      
      # 如果合并后的样本需要重新运行UMAP，可以取消下面的注释
      # merged_sample <- RunUMAP(merged_sample, dims = 1:30)
      
      # 保存合并后的样本
      merged_samples[[merged_name]] <- merged_sample
      
      cat("Successfully merged", paste(sample_pair, collapse = " and "), "into", merged_name, "\n")
    } else {
      cat("Warning: UMAP not found for", paste(sample_pair, collapse = " or "), ", skipping merge\n")
    }
  } else {
    cat("Warning: Not all samples found for", merged_name, ", skipping merge\n")
  }
}

# 4. 为合并后的样本创建dimplot
merged_dimplot_list <- list()
merged_sample_names <- names(merged_samples)

for(i in 1:length(merged_samples)) {
  sample_name <- merged_sample_names[i]
  
  # 检查该样本是否有UMAP
  if(!"umap" %in% names(merged_samples[[i]]@reductions)) {
    cat("Warning: No UMAP found for merged", sample_name, ", skipping...\n")
    next
  }
  
  # 检查是否有cell_type列
  if(!group_by_col %in% colnames(merged_samples[[i]]@meta.data)) {
    cat("Warning: No", group_by_col, "found for merged", sample_name, ", using seurat_clusters\n")
    current_group_col <- "seurat_clusters"
  } else {
    current_group_col <- group_by_col
  }
  
  # 创建dimplot
  p <- DimPlot(merged_samples[[i]], 
               reduction = "umap", 
               group.by = current_group_col,
               label = TRUE, 
               label.size = 3,
               pt.size = 0.5,
               cols = custom_colors) +  # 使用自定义颜色
    ggtitle(sample_name) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 14))
  
  merged_dimplot_list[[i]] <- p
  
  cat("Created dimplot for merged", sample_name, "\n")
}

# 移除NULL元素
merged_dimplot_list <- merged_dimplot_list[!sapply(merged_dimplot_list, is.null)]

# 5. 组合显示合并后的样本dimplot
if(length(merged_dimplot_list) > 0) {
  # 计算合适的行列数
  n_plots <- length(merged_dimplot_list)
  n_cols <- 3  # 合并后样本较少，使用3列
  n_rows <- ceiling(n_plots / n_cols)
  
  merged_combined_plot <- wrap_plots(merged_dimplot_list, ncol = n_cols)
  
  # 保存合并样本的组合图
  ggsave("merged_samples_dimplot_custom_colors.pdf", merged_combined_plot, 
         width = 5 * n_cols, height = 4 * n_rows)
  ggsave("merged_samples_dimplot_custom_colors.png", merged_combined_plot, 
         width = 5 * n_cols, height = 4 * n_rows, dpi = 300)
  
  cat("Saved merged sample plots with custom colors\n")
}

# 6. 创建一个带图例的参考图
if("REF" %in% names(all_samples)) {
  if("umap" %in% names(all_samples[["REF"]]@reductions)) {
    legend_plot <- DimPlot(all_samples[["REF"]], 
                           reduction = "umap", 
                           group.by = group_by_col,
                           label = TRUE, 
                           label.size = 3,
                           pt.size = 0.5,
                           cols = custom_colors) +
      ggtitle("Cell Type Legend") +
      theme(plot.title = element_text(hjust = 0.5, size = 16),
            legend.text = element_text(size = 12))
    
    ggsave("cell_type_legend.pdf", legend_plot, width = 10, height = 8)
    ggsave("cell_type_legend.png", legend_plot, width = 10, height = 8, dpi = 300)
    
    cat("Saved legend plot\n")
  }
}

cat("All plots completed!\n")

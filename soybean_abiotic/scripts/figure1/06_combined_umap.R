# Source: 图1.docx, extracted text lines 1148-1321
# 续接上一步的组合图与坐标修复版本
# See README.md and docs/KNOWN_LIMITATIONS.md before running.

###修改更好的可视化 修复9个图问题
# 加载必要的包
library(Seurat)
library(ggplot2)
library(gridExtra)
library(dplyr)

# 定义颜色配置
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

# 定义用于分组的列名（通常是细胞类型）
group_by_col <- "cell_type"  # 或者使用 "seurat_clusters" 如果没有cell_type

# 合并所有样本到一个列表中以便查找
all_samples <- c(list("REF" = reference), hormone.list)

# 定义要合并的样本对
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

# 修复现有的合并样本，添加UMAP信息
cat("Fixing UMAP information for merged samples...\n")

for(merged_name in names(merge_pairs)) {
  if(merged_name %in% names(merged_samples)) {
    sample_pair <- merge_pairs[[merged_name]]
    
    # 获取原始样本
    sample1 <- all_samples[[sample_pair[1]]]
    sample2 <- all_samples[[sample_pair[2]]]
    
    if("umap" %in% names(sample1@reductions) && "umap" %in% names(sample2@reductions)) {
      # 获取合并后的样本
      merged_sample <- merged_samples[[merged_name]]
      
      # 手动合并UMAP坐标
      umap1 <- sample1@reductions$umap@cell.embeddings
      umap2 <- sample2@reductions$umap@cell.embeddings
      
      # 为UMAP坐标添加细胞ID前缀
      rownames(umap1) <- paste0(sample_pair[1], "_", rownames(umap1))
      rownames(umap2) <- paste0(sample_pair[2], "_", rownames(umap2))
      
      # 合并UMAP坐标
      combined_umap <- rbind(umap1, umap2)
      
      # 确保UMAP坐标与合并后样本的细胞顺序匹配
      if(all(colnames(merged_sample) %in% rownames(combined_umap))) {
        combined_umap <- combined_umap[colnames(merged_sample), ]
        
        # 创建新的UMAP reduction对象
        merged_sample[["umap"]] <- CreateDimReducObject(
          embeddings = combined_umap,
          key = "UMAP_",
          assay = DefaultAssay(merged_sample)
        )
        
        # 更新合并样本
        merged_samples[[merged_name]] <- merged_sample
        
        cat("Fixed UMAP for", merged_name, "\n")
      } else {
        cat("Warning: Cell names don't match for", merged_name, "\n")
      }
    }
  }
}

# 创建DimPlot图形
merged_dimplot_list <- list()

for(sample_name in names(merged_samples)) {
  sample_obj <- merged_samples[[sample_name]]
  
  # 检查该样本是否有UMAP
  if(!"umap" %in% names(sample_obj@reductions)) {
    cat("Warning: No UMAP found for", sample_name, ", skipping...\n")
    next
  }
  
  # 检查是否有cell_type列
  if(!group_by_col %in% colnames(sample_obj@meta.data)) {
    cat("Warning: No", group_by_col, "found for", sample_name, ", using seurat_clusters\n")
    current_group_col <- "seurat_clusters"
  } else {
    current_group_col <- group_by_col
  }
  
  # 创建dimplot
  p <- DimPlot(sample_obj, 
               reduction = "umap", 
               group.by = current_group_col,
               label = FALSE,  # 关闭标签避免拥挤
               pt.size = 0.3,
               cols = custom_colors) +  # 使用自定义颜色
    ggtitle(sample_name) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 8))
  
  merged_dimplot_list[[sample_name]] <- p
  
  cat("Created dimplot for", sample_name, "\n")
}

# 排列和显示图形
if(length(merged_dimplot_list) > 0) {
  # 按照指定顺序排列图形
  plot_order <- c("REF", "ABA", "eBL", "highTemp", 
                  "IAA", "lowTemp", "MEJA", "PEG", "SA")
  
  # 重新排序图形 - 只包含存在的图形
  ordered_plots <- list()
  for(name in plot_order) {
    if(name %in% names(merged_dimplot_list)) {
      ordered_plots[[length(ordered_plots) + 1]] <- merged_dimplot_list[[name]]
    }
  }
  
  # 创建3x3网格布局
  if(length(ordered_plots) > 0) {
    grid_plot <- do.call(grid.arrange, c(ordered_plots, ncol = 3))
    
    # 保存图形
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

# 打印合并结果摘要
cat("\n=== 合并结果摘要 ===\n")
cat("成功合并的样本数量:", length(merged_samples), "\n")
for(name in names(merged_samples)) {
  n_cells <- ncol(merged_samples[[name]])
  has_umap <- "umap" %in% names(merged_samples[[name]]@reductions)
  cat(name, ":", n_cells, "个细胞, UMAP:", has_umap, "\n")
}
cat("Successfully transferred UMAP coordinates to merged object\n")
# 保存为RData文件
save(merged_samples, file = "merged_samples.RData")
cat("合并样本已保存为 merged_samples.RData\n")

# 也可以单独保存为RDS文件（更推荐，文件更小）
saveRDS(merged_samples, file = "merged_samples.rds")
cat("合并样本已保存为 merged_samples.rds\n")

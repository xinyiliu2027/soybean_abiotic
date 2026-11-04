# Source: 图1.docx, extracted text lines 40-280
# 需要已有 Gm_CK_integrated 对象
# See README.md and docs/KNOWN_LIMITATIONS.md before running.

##锚点注释
library(Seurat)
library(ggplot2)

# 合并Seurat v5数据层
Gm_CK_integrated <- JoinLayers(Gm_CK_integrated)
cat("已执行JoinLayers，数据层已合并。\n")

# 1. 基因列表（包含所有标记及新增/保留的Procambium）
soy_genes <- c(
  # Root Cap
  "Glyma.08G161300", "Glyma.15G266500", "Glyma.10G022500",
  # Lateral Root Cap
  "Glyma.12G003200", "Glyma.09G233600",
  # Cortex
  "Glyma.05G103600", "Glyma.05G057500", "Glyma.17G140000",
  "Glyma.04G091800", "Glyma.06G093600",
  # Dividing Cell
  "Glyma.15G040200",
  # Lateral Root
  "Glyma.04G252300", "Glyma.06G110200",
  # Endodermis
  "Glyma.13G307000", "Glyma.12G195500", "Glyma.12G195600", "Glyma.13G306900",
  # Meristem
  "Glyma.19G127200", "Glyma.12G021800", "Glyma.03G123600", 
  "Glyma.11G084500", "Glyma.11G095700",
  # Pericycle
  "Glyma.11G036300", "Glyma.08G271600", "Glyma.13G090700",
  # Xylem
  "Glyma.10G205000", "Glyma.01G088300", "Glyma.20G185500",
  # Phloem
  "Glyma.08G163500", "Glyma.15G263700", "Glyma.12G234200",
  # QC
  "Glyma.18G244600", "Glyma.09G248200",
  # Procambium (保留Pro-3，新增Pro-4)
  "Glyma.06G256000",        # Pro-3
  "Glyma.13G091600",        # Pro-4
  # Epidermis marker (原位，显示ID)
  "Glyma.11G147700",
  # GmWRKY72
  "Glyma.06G190800"
)

# 2. 基因映射（Procambium映射，表皮不映射）
gene_mapping <- c(
  "Glyma.08G161300" = "ANACO33-1",
  "Glyma.15G266500" = "ANACO33-2",
  "Glyma.10G022500" = "Glyma.10G022500", 
  "Glyma.12G003200" = "BRN1-1",
  "Glyma.09G233600" = "BRN1-2",
  "Glyma.05G103600" = "Cortex-1",
  "Glyma.05G057500" = "GmAZI5-1",
  "Glyma.17G140000" = "GmAZI5-2",
  "Glyma.04G091800" = "Cortex AED3-2",
  "Glyma.06G093600" = "Cortex AED3-4",
  "Glyma.15G040200" = "GmHTA12",
  "Glyma.04G252300" = "GmLAX3-1",
  "Glyma.06G110200" = "GmLAX3-2",
  "Glyma.13G307000" = "GmPER3-1",
  "Glyma.12G195500" = "GmPER3-2",
  "Glyma.12G195600" = "GmPER39-1",
  "Glyma.13G306900" = "GmPER39-2",
  "Glyma.19G127200" = "GmCYCB1-2-1",
  "Glyma.12G021800" = "GmUBC19-1",
  "Glyma.03G123600" = "GmCYCB1-2-2",
  "Glyma.11G084500" = "GmCYCB2-3",
  "Glyma.11G095700" = "GmUBC19-2",
  "Glyma.11G036300" = "GmNUDT8",
  "Glyma.08G271600" = "Glyma.08G271600",
  "Glyma.13G090700" = "Glyma.13G090700",
  "Glyma.10G205000" = "GmBHLH30-1",
  "Glyma.01G088300" = "GmBHLH30-2",
  "Glyma.20G185500" = "GmBHLH30-3",
  "Glyma.08G163500" = "GmAPL-1",
  "Glyma.15G263700" = "GmAPL-2",
  "Glyma.12G234200" = "GmTPS1",
  "Glyma.18G244600" = "BBM-2",
  "Glyma.09G248200" = "BBM-3",
  # Procambium

  # 表皮不映射，保留ID
  "Glyma.06G190800" = "GmWRKY72"
)

# 原位基因（绿色星号）
in_situ_genes <- c(
  "Glyma.05G103600", "Glyma.15G040200", "Glyma.12G234200",
  "Glyma.08G163500", "Glyma.10G022500", "Glyma.11G147700"
)

# 检查可用基因
available_genes <- soy_genes[soy_genes %in% rownames(Gm_CK_integrated)]
cat("匹配的基因数量:", length(available_genes), "\n")

# 3. 设置cluster
DefaultAssay(Gm_CK_integrated) <- "RNA"
cluster_order <- c("11", "7", "6", "5", "14", "9", "8", "15", "10", "13", 
                   "18", "17", "0", "4", "12", "3", "1", "2", "16")
Gm_CK_integrated@meta.data$seurat_clusters <- factor(
  Gm_CK_integrated@meta.data$seurat_clusters, levels = cluster_order
)

# 4. 细胞类型分组（仅用于分类）
cell_types <- list(
  "QC (静止中心)" = c("Glyma.18G244600", "Glyma.09G248200"),
  "根冠 (Root Cap)" = c("Glyma.08G161300", "Glyma.15G266500", "Glyma.10G022500"),
  "侧根冠 (Lateral Root Cap)" = c("Glyma.12G003200", "Glyma.09G233600"),
  "皮层 (Cortex)" = c("Glyma.05G103600", "Glyma.05G057500", "Glyma.17G140000",
                    "Glyma.04G091800", "Glyma.06G093600"),
  "分裂细胞 (Dividing Cell)" = c("Glyma.15G040200"),
  "侧根 (Lateral Root)" = c("Glyma.04G252300", "Glyma.06G110200"),
  "内皮层 (Endodermis)" = c("Glyma.13G307000", "Glyma.12G195500", "Glyma.12G195600", "Glyma.13G306900"),
  "分生组织 (Meristem)" = c("Glyma.19G127200", "Glyma.12G021800", "Glyma.03G123600", 
                        "Glyma.11G084500", "Glyma.11G095700"),
  "中柱鞘 (Pericycle)" = c("Glyma.11G036300", "Glyma.08G271600", "Glyma.13G090700"),
  "木质部 (Xylem)" = c("Glyma.10G205000", "Glyma.01G088300", "Glyma.20G185500"),
  "韧皮部 (Phloem)" = c("Glyma.08G163500", "Glyma.15G263700","Glyma.12G234200")
)

# 5. 生成初始排序：非皮层在前，皮层在底部
cortex_genes <- cell_types[["皮层 (Cortex)"]]
cortex_avail <- intersect(cortex_genes, available_genes)
non_cortex_avail <- setdiff(available_genes, cortex_avail)
reordered_genes <- c(non_cortex_avail, cortex_avail)

# ========== 调整基因顺序 ==========
# 已有顺序：... GmLAX3-2 (Glyma.06G110200), GmPER3-1 (Glyma.13G307000) ...
# 将 Pro-3 和 Pro-4 插入到 GmLAX3-2 之后，GmPER3-1 之前
lax3_2 <- "Glyma.06G110200"
per3_1 <- "Glyma.13G307000"
pro3 <- "Glyma.06G256000"
pro4 <- "Glyma.13G091600"
target_genes <- c(pro3, pro4)

if (all(c(lax3_2, per3_1, pro3, pro4) %in% reordered_genes)) {
  # 先从列表中移除这两个基因（保持原顺序）
  reordered_genes <- reordered_genes[!reordered_genes %in% target_genes]
  # 找到 lax3_2 的位置
  pos_lax <- which(reordered_genes == lax3_2)
  # 将两个基因插入到 lax3_2 之后（即 pos_lax + 1）
  # 先插入 pro3，再插入 pro4，顺序可自定（此处按pro3, pro4）
  reordered_genes <- append(reordered_genes, pro3, after = pos_lax)
  pos_pro3 <- which(reordered_genes == pro3)
  reordered_genes <- append(reordered_genes, pro4, after = pos_pro3)
}

# 额外保留之前的顺序：BRN1-2 -> GmWRKY72 -> Glyma.11G147700 -> GmHTA12
# 这段已经在上次调整中实现，此处再次执行（以防被覆盖）
brn1_2 <- "Glyma.09G233600"
wrky72 <- "Glyma.06G190800"
epi_marker <- "Glyma.11G147700"
if (all(c(brn1_2, wrky72, epi_marker) %in% reordered_genes)) {
  # 如果它们顺序不对，可重新调整
  # 但因为我们没有改变这些基因，通常顺序依然保持，为了安全，再执行一次
  reordered_genes <- reordered_genes[!reordered_genes %in% c(wrky72, epi_marker)]
  pos_brn <- which(reordered_genes == brn1_2)
  reordered_genes <- append(reordered_genes, wrky72, after = pos_brn)
  pos_wrky <- which(reordered_genes == wrky72)
  reordered_genes <- append(reordered_genes, epi_marker, after = pos_wrky)
}

# 6. 创建基因标签（原位加星号）
gene_labels <- sapply(reordered_genes, function(g) {
  label <- gene_mapping[g]
  if (is.na(label)) label <- g
  if (g %in% in_situ_genes) label <- paste0(label, " *")
  return(label)
})
names(gene_labels) <- reordered_genes

# 原位颜色和字体
label_colors <- sapply(gene_labels, function(l) if(grepl("\\*$", l)) "green" else "black")
label_faces <- sapply(gene_labels, function(l) if(grepl("\\*$", l)) "bold" else "plain")

# 7. 设置细胞类型元数据（用于分组）
cluster_to_celltype <- c(
  "0" = "Cortex", "1" = "Cortex", "2" = "QC", "3" = "Cortex", "4" = "Cortex",
  "5" = "Procambium", "6" = "Dividing Cell", "7" = "Root Cap", "8" = "Meristem",
  "9" = "Endodermis", "10" = "Pericycle", "11" = "Epidermis", "12" = "Cortex",
  "13" = "Pericycle", "14" = "Procambium", "15" = "Meristem", "16" = "Cortex",
  "17" = "Phloem", "18" = "Xylem"
)
Gm_CK_integrated@meta.data$cell_type <- cluster_to_celltype[as.character(Gm_CK_integrated@meta.data$seurat_clusters)]
celltype_order <- c("Root Cap", "Epidermis", "Dividing Cell", "Procambium",
                    "Endodermis", "Meristem", "Pericycle", "Xylem", "Phloem", "QC", "Cortex")
Gm_CK_integrated@meta.data$cell_type <- factor(Gm_CK_integrated@meta.data$cell_type, levels = celltype_order)

# ======================== 横向版本（基因在x轴，细胞类型在y轴） ========================
# 去除 coord_flip()，即为横向
p_horizontal <- DotPlot(Gm_CK_integrated, features = reordered_genes, group.by = "cell_type") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 12, color = label_colors, face = label_faces, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 18, color = "black"),
    plot.title = element_text(size = 20, color = "black", hjust = 0.5),
    legend.title = element_text(size = 16, color = "black"),
    legend.text = element_text(size = 14, color = "black"),
    plot.margin = margin(t = 20, r = 20, b = 40, l = 20, unit = "pt")
  ) +
  scale_color_gradientn(values = seq(0,1,0.2), colours = c("#FEEDB0","#F5A773","#E26253","#B32E5F","#721A60","#2F0F3E")) +
  labs(x = NULL, y = NULL, title = "大豆根尖细胞标记基因 (横向)") +
  guides(size = guide_legend(order = 3)) +
  scale_x_discrete(labels = gene_labels) +
  scale_y_discrete(limits = celltype_order, labels = celltype_order) +
  scale_size(range = c(1, 6))

# 保存横向版本（图幅缩小，宽度12，高度6）
ggsave("DotPlot_horizontal.pdf", plot = p_horizontal, width = 12, height = 6, device = "pdf")
print(p_horizontal)

# ======================== （可选）纵向版本，如需要可取消注释 ========================
 p_vertical <- DotPlot(Gm_CK_integrated, features = reordered_genes, group.by = "cell_type") +
   coord_flip() +
   theme_bw() +
   theme(
     panel.grid = element_blank(),
     axis.text.x = element_text(size = 14, color = "black", angle = 45, hjust = 1),
     axis.text.y = element_text(size = 12, color = label_colors, face = label_faces, margin = margin(r = 5)),
     axis.title = element_text(size = 18, color = "black"),
     plot.title = element_text(size = 20, color = "black", hjust = 0.5),
     legend.title = element_text(size = 16, color = "black"),
     legend.text = element_text(size = 14, color = "black"),
     plot.margin = margin(t = 20, r = 20, b = 40, l = 120, unit = "pt")
   ) +
   scale_color_gradientn(values = seq(0,1,0.2), colours = c("#FEEDB0","#F5A773","#E26253","#B32E5F","#721A60","#2F0F3E")) +
   labs(x = NULL, y = NULL, title = "大豆根尖细胞标记基因 (纵向)") +
   guides(size = guide_legend(order = 3)) +
   scale_x_discrete(labels = gene_labels) +
   scale_y_discrete(limits = celltype_order, labels = celltype_order)
 ggsave("DotPlot_vertical.pdf", plot = p_vertical, width = 8, height = 10, device = "pdf")
 print(p_vertical)

# 统计信息
cat("\n基因顺序（从x轴左到右）:\n")
cat(paste(gene_labels, collapse = " -> "), "\n")
cat("\nProcambium基因已插入到侧根基因 GmLAX3-2 之后、内皮层基因 GmPER3-1 之前。\n")
cat("横向版本已保存为 DotPlot_horizontal.pdf (宽12英寸，高6英寸)\n")

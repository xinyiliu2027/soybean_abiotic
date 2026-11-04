# Source: Figure 1 document, extracted text lines 40-280
# Requires an existing Gm_CK_integrated object.
# See README.md and docs/KNOWN_LIMITATIONS.md before running.

# Marker-based reference annotation.
library(Seurat)
library(ggplot2)

# Join Seurat v5 assay layers.
Gm_CK_integrated <- JoinLayers(Gm_CK_integrated)
cat("JoinLayers completed; assay layers merged.\n")

# 1. Gene list, including retained and added Procambium markers.
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
  # Procambium: retain Pro-3 and add Pro-4.
  "Glyma.06G256000",        # Pro-3
  "Glyma.13G091600",        # Pro-4
  # Epidermis marker: in situ marker displayed using its gene ID.
  "Glyma.11G147700",
  # GmWRKY72
  "Glyma.06G190800"
)

# 2. Gene-label mapping; retain IDs for unmapped Epidermis markers.
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

  # Do not remap the Epidermis marker; retain its ID.
  "Glyma.06G190800" = "GmWRKY72"
)

# In situ markers, indicated with green asterisks.
in_situ_genes <- c(
  "Glyma.05G103600", "Glyma.15G040200", "Glyma.12G234200",
  "Glyma.08G163500", "Glyma.10G022500", "Glyma.11G147700"
)

# Check which genes are available.
available_genes <- soy_genes[soy_genes %in% rownames(Gm_CK_integrated)]
cat("Number of matched genes:", length(available_genes), "\n")

# 3. Set the cluster order.
DefaultAssay(Gm_CK_integrated) <- "RNA"
cluster_order <- c("11", "7", "6", "5", "14", "9", "8", "15", "10", "13", 
                   "18", "17", "0", "4", "12", "3", "1", "2", "16")
Gm_CK_integrated@meta.data$seurat_clusters <- factor(
  Gm_CK_integrated@meta.data$seurat_clusters, levels = cluster_order
)

# 4. Group markers by cell type for classification only.
cell_types <- list(
  "QC (Quiescent Center)" = c("Glyma.18G244600", "Glyma.09G248200"),
  "Root Cap" = c("Glyma.08G161300", "Glyma.15G266500", "Glyma.10G022500"),
  "Lateral Root Cap" = c("Glyma.12G003200", "Glyma.09G233600"),
  "Cortex" = c("Glyma.05G103600", "Glyma.05G057500", "Glyma.17G140000",
                    "Glyma.04G091800", "Glyma.06G093600"),
  "Dividing Cell" = c("Glyma.15G040200"),
  "Lateral Root" = c("Glyma.04G252300", "Glyma.06G110200"),
  "Endodermis" = c("Glyma.13G307000", "Glyma.12G195500", "Glyma.12G195600", "Glyma.13G306900"),
  "Meristem" = c("Glyma.19G127200", "Glyma.12G021800", "Glyma.03G123600", 
                        "Glyma.11G084500", "Glyma.11G095700"),
  "Pericycle" = c("Glyma.11G036300", "Glyma.08G271600", "Glyma.13G090700"),
  "Xylem" = c("Glyma.10G205000", "Glyma.01G088300", "Glyma.20G185500"),
  "Phloem" = c("Glyma.08G163500", "Glyma.15G263700","Glyma.12G234200")
)

# 5. Initial order: non-Cortex markers first and Cortex markers last.
cortex_genes <- cell_types[["Cortex"]]
cortex_avail <- intersect(cortex_genes, available_genes)
non_cortex_avail <- setdiff(available_genes, cortex_avail)
reordered_genes <- c(non_cortex_avail, cortex_avail)

# Adjust the gene order.
# Existing order includes GmLAX3-2 (Glyma.06G110200), then GmPER3-1 (Glyma.13G307000).
# Insert Pro-3 and Pro-4 after GmLAX3-2 and before GmPER3-1.
lax3_2 <- "Glyma.06G110200"
per3_1 <- "Glyma.13G307000"
pro3 <- "Glyma.06G256000"
pro4 <- "Glyma.13G091600"
target_genes <- c(pro3, pro4)

if (all(c(lax3_2, per3_1, pro3, pro4) %in% reordered_genes)) {
  # Remove these two genes first while preserving the order of the remaining genes.
  reordered_genes <- reordered_genes[!reordered_genes %in% target_genes]
  # Find the position of lax3_2.
  pos_lax <- which(reordered_genes == lax3_2)
  # Insert the two genes after lax3_2, starting at pos_lax + 1.
  # Insert pro3 first, then pro4.
  reordered_genes <- append(reordered_genes, pro3, after = pos_lax)
  pos_pro3 <- which(reordered_genes == pro3)
  reordered_genes <- append(reordered_genes, pro4, after = pos_pro3)
}

# Also preserve BRN1-2 -> GmWRKY72 -> Glyma.11G147700 -> GmHTA12.
# Reapply the previous ordering in case it was overwritten.
brn1_2 <- "Glyma.09G233600"
wrky72 <- "Glyma.06G190800"
epi_marker <- "Glyma.11G147700"
if (all(c(brn1_2, wrky72, epi_marker) %in% reordered_genes)) {
  # Reorder these genes if necessary.
  # Their relative order should be unchanged; reapply it explicitly as a safeguard.
  reordered_genes <- reordered_genes[!reordered_genes %in% c(wrky72, epi_marker)]
  pos_brn <- which(reordered_genes == brn1_2)
  reordered_genes <- append(reordered_genes, wrky72, after = pos_brn)
  pos_wrky <- which(reordered_genes == wrky72)
  reordered_genes <- append(reordered_genes, epi_marker, after = pos_wrky)
}

# 6. Create gene labels, adding asterisks to in situ markers.
gene_labels <- sapply(reordered_genes, function(g) {
  label <- gene_mapping[g]
  if (is.na(label)) label <- g
  if (g %in% in_situ_genes) label <- paste0(label, " *")
  return(label)
})
names(gene_labels) <- reordered_genes

# Colors and font styles for in situ markers.
label_colors <- sapply(gene_labels, function(l) if(grepl("\\*$", l)) "green" else "black")
label_faces <- sapply(gene_labels, function(l) if(grepl("\\*$", l)) "bold" else "plain")

# 7. Set cell-type metadata for grouping.
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

# Horizontal layout: genes on the x-axis and cell types on the y-axis.
# Omit coord_flip() for this orientation.
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
  labs(x = NULL, y = NULL, title = "Soybean root-tip cell markers (horizontal)") +
  guides(size = guide_legend(order = 3)) +
  scale_x_discrete(labels = gene_labels) +
  scale_y_discrete(limits = celltype_order, labels = celltype_order) +
  scale_size(range = c(1, 6))

# Save the horizontal layout at 12 x 6 inches.
ggsave("DotPlot_horizontal.pdf", plot = p_horizontal, width = 12, height = 6, device = "pdf")
print(p_horizontal)

# Optional vertical layout; enable the relevant code if needed.
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
   labs(x = NULL, y = NULL, title = "Soybean root-tip cell markers (vertical)") +
   guides(size = guide_legend(order = 3)) +
   scale_x_discrete(labels = gene_labels) +
   scale_y_discrete(limits = celltype_order, labels = celltype_order)
 ggsave("DotPlot_vertical.pdf", plot = p_vertical, width = 8, height = 10, device = "pdf")
 print(p_vertical)

# Report summary information.
cat("\nGene order from left to right:\n")
cat(paste(gene_labels, collapse = " -> "), "\n")
cat("\nProcambium genes inserted after GmLAX3-2 and before GmPER3-1.\n")
cat("Horizontal plot saved as DotPlot_horizontal.pdf (12 x 6 inches).\n")

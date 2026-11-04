# Source: 图4.docx, extracted text lines 1824-1988
# 读取01生成的consensus_network_object.rds
# See README.md and docs/KNOWN_LIMITATIONS.md before running.

##可视化模块基因图
suppressPackageStartupMessages({
  library(Seurat)
  library(hdWGCNA)
  library(ggplot2)
})

OUT_DIR <- file.path(Sys.getenv("SOYBEAN_PROJECT_ROOT", unset = "data"), "CK2_added_v3_Results/Ensemble_PI_Heat_Cold_PEG_with_CK2")
FILE_PREFIX <- "Endo3_v5_"
NETWORK_NAME <- "AllCells_Endo3_consensus_v5"

## ---------------------------------------------------------------------------
## 1. 从 RDS 读取模块表
## ---------------------------------------------------------------------------
net_obj <- readRDS(file.path(OUT_DIR, paste0(FILE_PREFIX, "consensus_network_object.rds")))

mods <- hdWGCNA::GetModules(net_obj, wgcna_name = NETWORK_NAME)
mods <- as.data.frame(mods, stringsAsFactors = FALSE)

# 强制关键列转为普通字符向量，避免 list 列干扰
mods$gene_name <- as.character(mods$gene_name)
mods$module    <- as.character(mods$module)

print(names(mods))
print(table(mods$module, useNA = "ifany"))

## ---------------------------------------------------------------------------
## 2. 统计非 grey 模块基因数（base R，不用 dplyr）
## ---------------------------------------------------------------------------
tab <- table(mods$module, useNA = "ifany")

counts <- data.frame(
  module  = names(tab),
  n_genes = as.integer(tab),
  stringsAsFactors = FALSE
)

# 去掉 grey / NA
counts <- counts[!is.na(counts$module) & counts$module != "grey", , drop = FALSE]

# M1, M2, ... 编号
counts$Module <- paste0("M", as.integer(sub("^StageProgram", "", counts$module)))

# 标出 StageProgram09
counts$Highlight <- counts$module == "StageProgram09"

# 按 M 编号排序
counts <- counts[order(as.integer(sub("^M", "", counts$Module))), , drop = FALSE]
rownames(counts) <- NULL

print(counts)

## ---------------------------------------------------------------------------
## 3. 画模块基因数柱状图，突出 StageProgram09
## ---------------------------------------------------------------------------
p <- ggplot(counts, aes(x = Module, y = n_genes, fill = Highlight)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = n_genes), vjust = -0.35, size = 3.5) +
  scale_fill_manual(
    values = c(`FALSE` = "#555555", `TRUE` = "#E6007E"),
    guide = "none"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(
    x = "Module",
    y = "Number of genes",
    title = "StageProgram module sizes",
    subtitle = "StageProgram09 / M9 highlighted"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text = element_text(color = "black"),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey30")
  )

ggsave(
  file.path(OUT_DIR, paste0(FILE_PREFIX, "Figure_module_gene_counts_StageProgram09.pdf")),
  p, width = 7, height = 4.2, useDingbats = FALSE
)
ggsave(
  file.path(OUT_DIR, paste0(FILE_PREFIX, "Figure_module_gene_counts_StageProgram09.png")),
  p, width = 7, height = 4.2, dpi = 300, bg = "white"
)

## ---------------------------------------------------------------------------
## 4. 导出 StageProgram09 基因
## ---------------------------------------------------------------------------
p9_genes <- unique(mods$gene_name[mods$module == "StageProgram09"])

writeLines(
  p9_genes,
  file.path(OUT_DIR, paste0(FILE_PREFIX, "StageProgram09_genes.txt"))
)

cat("StageProgram09 基因数：", length(p9_genes), "\n")

## ---------------------------------------------------------------------------
## 5. 定位 P1/P2/P3 候选基因的模块归属
## ---------------------------------------------------------------------------
cand <- data.frame(
  Stage = c("P1","P1","P1","P1",
            "P2","P2","P2",
            "P3"),
  Gene = c("Glyma.03G218700","Glyma.04G122301","Glyma.07G125100","Glyma.10G280000",
           "Glyma.09G281000","Glyma.13G321700","Glyma.17G152800",
           "Glyma.13G321700"),
  stringsAsFactors = FALSE
)

# 模块表若同一基因多行，先去重
mods_small <- unique(mods[, c("gene_name", "module")])

loc <- merge(
  cand,
  mods_small,
  by.x = "Gene",
  by.y = "gene_name",
  all.x = TRUE,
  sort = FALSE
)

loc$In_network        <- !is.na(loc$module)
loc$In_non_grey       <- !is.na(loc$module) & loc$module != "grey"
loc$In_StageProgram09 <- !is.na(loc$module) & loc$module == "StageProgram09"

loc <- loc[order(loc$Stage, loc$Gene), , drop = FALSE]
rownames(loc) <- NULL

write.csv(
  loc,
  file.path(OUT_DIR, paste0(FILE_PREFIX, "candidate_gene_module_location.csv")),
  row.names = FALSE
)

print(loc)

## ---------------------------------------------------------------------------
## 6. 候选基因与 StageProgram09 的交集
## ---------------------------------------------------------------------------
overlap_p9 <- data.frame(
  Gene = unique(cand$Gene),
  In_StageProgram09 = unique(cand$Gene) %in% p9_genes,
  stringsAsFactors = FALSE
)

write.csv(
  overlap_p9,
  file.path(OUT_DIR, paste0(FILE_PREFIX, "StageProgram09_candidate_overlap.csv")),
  row.names = FALSE
)

print(overlap_p9)

## ---------------------------------------------------------------------------
## 7. 若某基因 module 为 NA，检查是否进入网络输入基因
## ---------------------------------------------------------------------------
input_file <- file.path(OUT_DIR, paste0(FILE_PREFIX, "network_input_genes.csv"))

if (file.exists(input_file)) {
  input_genes <- read.csv(input_file, stringsAsFactors = FALSE)
  cand$In_network_input <- cand$Gene %in% input_genes$Gene
  print(cand)
}

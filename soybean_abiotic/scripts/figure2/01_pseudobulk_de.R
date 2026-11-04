# Source: 图2.docx, extracted text lines 1-474
# 独立运行；合并表达RDS和CK RData
# See README.md and docs/KNOWN_LIMITATIONS.md before running.

suppressMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(muscat)
  library(edgeR)
  library(limma)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(cowplot)
})

## ==============================
## 0) 路径与参数
## ==============================
MERGED_RDS <- file.path(Sys.getenv("SOYBEAN_PROJECT_ROOT", unset = "data"), "seurat_obj_with_celltype.rds")
CK_RDATA   <- file.path(Sys.getenv("SOYBEAN_PROJECT_ROOT", unset = "data"), "data_files/Gm_CK_integrated0.8.RData")  # REF 来源

OUT_DIR <- file.path(Sys.getenv("SOYBEAN_PROJECT_ROOT", unset = "data"), "muscat_results_ck_ref_paper")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

## ===== 按文献阈值（paper）=====
## 文献：FDR <= 0.05；|FC| >= 1.5；detection frequency >= 10%（任一条件）
FDR_CUT <- 0.05
LFC_CUT <- log2(1.5)   # ≈ 0.5849625
FRQ_CUT <- 0.10

## pbDS 设置（文献：edgeR + pbDS + aggregateData sum counts）
PBDS_METHOD <- "edgeR"
PBDS_FILTER <- "none"

## pbDS 内部过滤（尽量不额外过滤，和文献一致）
PBDS_MIN_CELLS_STRICT <- 10
PBDS_MIN_CELLS_LOOSE  <- 1

cat("========================================\n")
cat("Paper-like 参数设置:\n")
cat("FDR_CUT =", FDR_CUT, "\n")
cat("LFC_CUT =", LFC_CUT, "(fold change ≥", 2^LFC_CUT, ")\n")
cat("FRQ_CUT =", FRQ_CUT, "\n")
cat("PBDS_METHOD =", PBDS_METHOD, "\n")
cat("PBDS_FILTER =", PBDS_FILTER, "\n")
cat("REF 来源 = CK RData 文件\n")
cat("========================================\n\n")

## ==============================
## 1) Seurat counts 读取兼容函数
## ==============================
get_counts_Seurat <- function(obj, assay = "RNA") {
  out <- tryCatch({
    LayerData(obj, assay = assay, layer = "counts")
  }, error = function(e) {
    tryCatch({
      GetAssayData(obj, assay = assay, layer = "counts")
    }, error = function(e2) {
      GetAssayData(obj, assay = assay, slot = "counts")
    })
  })
  return(out)
}

pick_assay_with_counts <- function(obj) {
  cand <- unique(c("RNA", "SCT", DefaultAssay(obj)))
  for (a in cand) {
    ok <- tryCatch({
      x <- get_counts_Seurat(obj, assay = a)
      !is.null(x) && nrow(x) > 0 && ncol(x) > 0
    }, error = function(e) FALSE)
    if (ok) return(a)
  }
  stop("找不到可用的 counts（尝试 RNA/SCT/DefaultAssay 都失败）")
}

## ==============================
## 2) 读取 merged（只用于非 REF 处理组）
## ==============================
cat("读取 merged Seurat: ", MERGED_RDS, "\n")
merged <- readRDS(MERGED_RDS)

stopifnot(all(c("Group1","orig.ident","cell_type") %in% colnames(merged@meta.data)))

cat("merged cells:", ncol(merged), " genes:", nrow(merged), "\n")
cat("\nGroup1 分布:\n"); print(table(merged$Group1))
cat("\ncell_type 分布:\n"); print(sort(table(merged$cell_type), decreasing = TRUE))

## 丢弃 merged 中所有 REF 细胞（完全不用 merged 的 REF）
merged_nonref <- subset(merged, subset = Group1 != "REF")
cat("\n[INFO] merged_nonref cells:", ncol(merged_nonref), "\n")
cat("[INFO] merged_nonref Group1:\n"); print(table(merged_nonref$Group1))

merged_nonref$sample_id_muscat <- as.character(merged_nonref$orig.ident)
merged_nonref$cluster_id_use   <- as.character(merged_nonref$cell_type)

## ==============================
## 3) 读取 CK RData 文件（作为 REF1/REF2）
## ==============================
cat("\n读取 CK RData 文件: ", CK_RDATA, "\n")

loaded_env <- new.env()
load(CK_RDATA, envir = loaded_env)

cat("[RData] 加载的对象:\n")
print(ls(envir = loaded_env))

## 找到第一个 Seurat 对象
obj_names <- ls(envir = loaded_env)
seurat_obj_name <- NULL
for (nm in obj_names) {
  obj <- get(nm, envir = loaded_env)
  if (inherits(obj, "Seurat")) { seurat_obj_name <- nm; break }
}
if (is.null(seurat_obj_name)) stop("RData 文件中未找到 Seurat 对象")

cat("[RData] 找到 Seurat 对象: ", seurat_obj_name, "\n")
ck_seurat <- get(seurat_obj_name, envir = loaded_env)

cat("[RData] Seurat 对象维度: ", ncol(ck_seurat), " cells, ", nrow(ck_seurat), " genes\n")
cat("[RData] 元数据列名: ", paste(colnames(ck_seurat@meta.data), collapse = ", "), "\n")

## 样本字段
sample_field <- NULL
for (field in c("orig.ident", "sample", "Sample", "group", "Group", "condition")) {
  if (field %in% colnames(ck_seurat@meta.data)) { sample_field <- field; break }
}
if (is.null(sample_field)) stop("未找到样本信息字段（orig.ident/sample/Sample/group/Group/condition）")
cat("[RData] 使用样本字段: ", sample_field, "\n")
print(table(ck_seurat@meta.data[[sample_field]]))

## 聚类字段
cluster_field <- NULL
for (field in c("seurat_clusters", "cluster", "ClusterID", "cell_type", "CellType")) {
  if (field %in% colnames(ck_seurat@meta.data)) { cluster_field <- field; break }
}
if (is.null(cluster_field)) stop("未找到聚类信息字段（seurat_clusters/cluster/ClusterID/cell_type/CellType）")
cat("[RData] 使用聚类字段: ", cluster_field, "\n")
print(table(ck_seurat@meta.data[[cluster_field]]))

## ==============================
## 4) 细胞类型注释映射（如为数字聚类）
## ==============================
cluster_info <- as.character(ck_seurat@meta.data[[cluster_field]])
unique_clusters <- unique(cluster_info)
cat("[RData] 独特聚类值: ", paste(sort(unique_clusters), collapse = ", "), "\n")

if (all(grepl("^\\d+$", unique_clusters))) {
  cat("[RData] 检测到数字聚类，进行细胞类型映射...\n")
  cluster_annotation <- c(
    "0"  = "Cortex","1"="Cortex","2"="QC","3"="Cortex","4"="Cortex",
    "5"  = "Procambium","6"="Dividing Cell","7"="Root Cap","8"="Meristem",
    "9"  = "Endodermis","10"="Pericycle","11"="Epidermis","12"="Cortex",
    "13" = "Pericycle","14"="Procambium","15"="Meristem","16"="Cortex",
    "17" = "Phloem","18"="Xylem"
  )
  cluster_mapped <- unname(cluster_annotation[as.character(cluster_info)])
  cluster_mapped[is.na(cluster_mapped)] <- "Unknown"
  print(sort(table(cluster_mapped), decreasing = TRUE))
} else {
  cluster_mapped <- cluster_info
  cat("[RData] 使用已有的细胞类型标签\n")
}

## ==============================
## 5) 创建 REF1/REF2 的 SCE
## ==============================
cat("\n构建 REF1/REF2 的 SCE...\n")

sample_info <- as.character(ck_seurat@meta.data[[sample_field]])
unique_samples <- unique(sample_info)
cat("[RData] 独特样本名称: ", paste(unique_samples, collapse = ", "), "\n")

ck_assay <- pick_assay_with_counts(ck_seurat)
cat("[RData] 使用 counts assay: ", ck_assay, "\n")
ck_counts <- get_counts_Seurat(ck_seurat, assay = ck_assay)

## 识别 CK1/CK2，找不到就用前两个样本
ck1_pattern <- "^CK1$|^ck1$|^CK_1$|^ck_1$|^control1$|^Control1$"
ck2_pattern <- "^CK2$|^ck2$|^CK_2$|^ck_2$|^control2$|^Control2$"
ck1_idx <- grep(ck1_pattern, sample_info, ignore.case = TRUE)
ck2_idx <- grep(ck2_pattern, sample_info, ignore.case = TRUE)

if (length(ck1_idx) == 0 || length(ck2_idx) == 0) {
  cat("[WARN] 未找到标准 CK1/CK2 命名，使用前两个样本作为 REF1/REF2\n")
  if (length(unique_samples) < 2) stop("RData 样本数 < 2，无法构建 REF1/REF2")
  ck1_name <- unique_samples[1]
  ck2_name <- unique_samples[2]
  ck1_idx <- which(sample_info == ck1_name)
  ck2_idx <- which(sample_info == ck2_name)
  cat("[INFO] REF1=", ck1_name, "；REF2=", ck2_name, "\n", sep = "")
}

cat("[INFO] REF1 细胞数:", length(ck1_idx), "\n")
cat("[INFO] REF2 细胞数:", length(ck2_idx), "\n")

counts_ref1 <- ck_counts[, ck1_idx, drop = FALSE]
sce_ref1 <- SingleCellExperiment(
  assays = list(counts = counts_ref1),
  colData = data.frame(
    group_id   = factor(rep("REF",  ncol(counts_ref1))),
    sample_id  = factor(rep("REF1", ncol(counts_ref1))),
    cluster_id = factor(cluster_mapped[ck1_idx]),
    row.names  = colnames(counts_ref1),
    stringsAsFactors = FALSE
  )
)

counts_ref2 <- ck_counts[, ck2_idx, drop = FALSE]
sce_ref2 <- SingleCellExperiment(
  assays = list(counts = counts_ref2),
  colData = data.frame(
    group_id   = factor(rep("REF",  ncol(counts_ref2))),
    sample_id  = factor(rep("REF2", ncol(counts_ref2))),
    cluster_id = factor(cluster_mapped[ck2_idx]),
    row.names  = colnames(counts_ref2),
    stringsAsFactors = FALSE
  )
)

cat("\n[CHECK] REF1 cluster_id:\n"); print(sort(table(colData(sce_ref1)$cluster_id), decreasing = TRUE))
cat("\n[CHECK] REF2 cluster_id:\n"); print(sort(table(colData(sce_ref2)$cluster_id), decreasing = TRUE))

## ==============================
## 6) merged_nonref -> SCE
## ==============================
cat("\n构建 merged_nonref 的 SCE...\n")
m_assay  <- pick_assay_with_counts(merged_nonref)
cat("[MERGED_NONREF] 使用 counts assay:", m_assay, "\n")
m_counts <- get_counts_Seurat(merged_nonref, assay = m_assay)

sce_nonref <- SingleCellExperiment(
  assays = list(counts = m_counts),
  colData = data.frame(
    group_id   = factor(as.character(merged_nonref$Group1)),
    sample_id  = factor(as.character(merged_nonref$sample_id_muscat)),
    cluster_id = factor(as.character(merged_nonref$cluster_id_use)),
    row.names  = colnames(m_counts),
    stringsAsFactors = FALSE
  )
)
cat("[MERGED_NONREF] cells:", ncol(sce_nonref), " genes:", nrow(sce_nonref), "\n")

## ==============================
## 7) 合并为全量 SCE（对齐基因）
## ==============================
cat("\n合并 SCE（对齐基因）...\n")
common_genes <- Reduce(intersect, list(rownames(sce_nonref), rownames(sce_ref1), rownames(sce_ref2)))
cat("共同基因数:", length(common_genes), "\n")
if (length(common_genes) < 1000) stop("共同基因太少，请检查基因ID是否一致")

sce_all <- cbind(
  sce_nonref[common_genes, ],
  sce_ref1[common_genes, ],
  sce_ref2[common_genes, ]
)

cat("合并后 SCE: cells =", ncol(sce_all), ", genes =", nrow(sce_all), "\n")
cat("\n合并后 group_id:\n"); print(table(colData(sce_all)$group_id))

## experiment_info（sample-level）
ei <- as.data.frame(colData(sce_all)) %>%
  as_tibble(rownames = "cell") %>%
  dplyr::select(sample_id, group_id) %>%
  distinct() %>%
  arrange(group_id, sample_id)

n_cells <- table(colData(sce_all)$sample_id)
ei$n_cells <- as.numeric(n_cells[match(ei$sample_id, names(n_cells))])
metadata(sce_all)$experiment_info <- as.data.frame(ei)

cat("\n[experiment_info]:\n")
print(metadata(sce_all)$experiment_info)

## ==============================
## 8) pbDS 安全调用 + 自动重试
## ==============================
call_pbDS_safe <- function(pb, design, contrast,
                           method = "edgeR",
                           filter = "none",
                           min_cells = NULL,
                           min_counts = NULL,
                           verbose = TRUE) {
  fml <- names(formals(muscat::pbDS))
  args <- list(pb = pb, design = design, contrast = contrast,
               method = method, filter = filter, verbose = verbose)
  if (!is.null(min_cells)  && "min_cells"  %in% fml) args$min_cells  <- min_cells
  if (!is.null(min_counts) && "min_counts" %in% fml) args$min_counts <- min_counts
  do.call(muscat::pbDS, args)
}

run_pbDS_with_retry <- function(pb, mm, contrast,
                                method = PBDS_METHOD,
                                filter = PBDS_FILTER,
                                min_cells_strict = PBDS_MIN_CELLS_STRICT,
                                min_cells_loose  = PBDS_MIN_CELLS_LOOSE) {
  msg <- NULL
  res <- tryCatch({
    call_pbDS_safe(pb, mm, contrast, method = method, filter = filter,
                   min_cells = min_cells_strict, min_counts = 1, verbose = TRUE)
  }, error = function(e) { msg <<- conditionMessage(e); NULL })
  
  if (!is.null(res)) return(list(res = res, used_min_cells = min_cells_strict, used_min_counts = 1))
  
  cat("\n[WARN] pbDS 严格阈值失败：\n  ", msg, "\n", sep = "")
  cat("[RETRY] min_cells=", min_cells_loose, "；min_counts=0\n", sep = "")
  
  msg2 <- NULL
  res2 <- tryCatch({
    call_pbDS_safe(pb, mm, contrast, method = method, filter = filter,
                   min_cells = min_cells_loose, min_counts = 0, verbose = TRUE)
  }, error = function(e) { msg2 <<- conditionMessage(e); NULL })
  
  if (!is.null(res2)) return(list(res = res2, used_min_cells = min_cells_loose, used_min_counts = 0))
  
  cat("\n[ERROR] pbDS 放宽阈值仍失败：\n  ", msg2, "\n", sep = "")
  return(list(res = NULL, used_min_cells = NA, used_min_counts = NA))
}

## ==============================
## 9) 逐个 treatment vs REF 做 pbDS + 按文献阈值过滤
## ==============================
treatments <- sort(setdiff(unique(as.character(colData(sce_all)$group_id)), "REF"))
cat("\n将进行以下比较（treatment vs REF）：\n")
print(treatments)

all_results <- list()

for (treatment in treatments) {
  cat("\n========================================\n")
  cat(">>> 分析: ", treatment, " vs REF\n", sep = "")
  cat("========================================\n")
  
  ## 9.1 子集（当前 treatment + REF）
  keep <- colData(sce_all)$group_id %in% c("REF", treatment)
  sce_sub <- sce_all[, keep]
  
  colData(sce_sub)$group_id   <- droplevels(colData(sce_sub)$group_id)
  colData(sce_sub)$sample_id  <- droplevels(colData(sce_sub)$sample_id)
  colData(sce_sub)$cluster_id <- droplevels(colData(sce_sub)$cluster_id)
  
  cat("子集细胞数:", ncol(sce_sub), "\n")
  cat("子集 group_id:\n"); print(table(colData(sce_sub)$group_id))
  
  ## 9.2 prepSCE
  cat("\n运行 prepSCE...\n")
  sce_sub <- prepSCE(
    sce_sub,
    kid  = "cluster_id",
    gid  = "group_id",
    sid  = "sample_id",
    drop = TRUE
  )
  
  ei_sub <- metadata(sce_sub)$experiment_info
  cat("\nexperiment_info:\n"); print(ei_sub)
  cat("\n每组样本数:\n"); print(table(ei_sub$group_id))
  
  ## 9.3 pseudobulk（sum raw counts）——文献一致
  cat("\n聚合到 pseudobulk...\n")
  pb <- aggregateData(
    sce_sub,
    assay = "counts",
    fun   = "sum",
    by    = c("cluster_id", "sample_id")
  )
  cat("pseudobulk 完成: ", length(assayNames(pb)), " 个细胞类型\n", sep = "")
  
  ## 记录 cluster×sample 细胞数（可用于 QC）
  nc_mat <- table(colData(sce_sub)$cluster_id, colData(sce_sub)$sample_id)
  nc_mat <- nc_mat[assayNames(pb), colnames(pb), drop = FALSE]
  metadata(pb)$n_cells <- nc_mat
  
  ## pb colData 对齐 group_id/sample_id
  ei2 <- metadata(sce_sub)$experiment_info
  colData(pb)$group_id  <- ei2$group_id[match(colnames(pb), ei2$sample_id)]
  colData(pb)$sample_id <- colnames(pb)
  
  ## 9.4 设计矩阵 + contrast
  pb_groups <- factor(as.character(colData(pb)$group_id), levels = c("REF", treatment))
  mm <- model.matrix(~0 + pb_groups)
  colnames(mm) <- levels(pb_groups)
  rownames(mm) <- colData(pb)$sample_id
  
  contrast_name <- paste0(treatment, "-REF")
  contrast <- limma::makeContrasts(contrasts = contrast_name, levels = mm)
  
  cat("\n设计矩阵(mm):\n"); print(mm)
  cat("\n对比矩阵(contrast):\n"); print(contrast)
  
  ## 9.5 pbDS（edgeR）——文献一致
  cat("\n运行 pbDS（自动重试）...\n")
  pbds_out <- run_pbDS_with_retry(pb, mm, contrast)
  if (is.null(pbds_out$res)) {
    cat("[SKIP] ", treatment, " vs REF：pbDS 无法运行，跳过\n", sep = "")
    next
  }
  res <- pbds_out$res
  cat("[OK] pbDS 成功：min_cells=", pbds_out$used_min_cells,
      "；min_counts=", pbds_out$used_min_counts, "\n", sep="")
  
  ## 9.6 提取结果（包含 frq）
  cat("\n提取结果...\n")
  tb_full <- resDS(sce_sub, res, bind = "row", cpm = FALSE, frq = TRUE)
  cat("全部结果行数:", nrow(tb_full), "\n")
  
  ## 9.7 按文献：用 FDR-adjusted P（优先 local，更贴合常见写法）
  if ("p_adj.loc" %in% colnames(tb_full)) {
    p_adj_col <- "p_adj.loc"
    cat("使用 p_adj.loc 进行过滤\n")
  } else if ("p_adj.glb" %in% colnames(tb_full)) {
    p_adj_col <- "p_adj.glb"
    cat("[WARN] p_adj.loc 不存在，使用 p_adj.glb\n")
  } else {
    stop("tb_full 中既没有 p_adj.loc 也没有 p_adj.glb")
  }
  
  ## 9.8 文献阈值过滤：
  ## FDR <= 0.05；|FC|>=1.5（即 |log2FC| >= log2(1.5)）
  tb_sig <- tb_full %>%
    dplyr::filter(
      .data[[p_adj_col]] <= FDR_CUT,
      abs(logFC) >= LFC_CUT
    )
  cat("通过 ", p_adj_col, " + |logFC| 过滤后: ", nrow(tb_sig), " 行\n", sep = "")
  
  ## 9.9 文献阈值：detection frequency >= 10%（任一条件）
  ref_frq_col <- "REF.frq"
  trt_frq_col <- paste0(treatment, ".frq")
  if (ref_frq_col %in% colnames(tb_sig) && trt_frq_col %in% colnames(tb_sig)) {
    tb_sig <- tb_sig %>%
      dplyr::filter(.data[[ref_frq_col]] >= FRQ_CUT | .data[[trt_frq_col]] >= FRQ_CUT)
    cat("通过频率过滤后:", nrow(tb_sig), "行\n")
  } else {
    cat("[WARN] 找不到 frq 列，跳过频率过滤（不建议）\n")
  }
  
  ## 9.10 Up/Down 标签（logFC>0 表示 treatment 高于 REF）
  tb_sig <- tb_sig %>%
    dplyr::mutate(
      regulation = dplyr::case_when(
        logFC >=  LFC_CUT ~ "Up",
        logFC <= -LFC_CUT ~ "Down",
        TRUE              ~ "Not_DE"
      ),
      comparison = paste0(treatment, "_vs_REF")
    )
  
  cat("\n========================================\n")
  cat("★ Paper-like 过滤后结果 ★\n")
  cat("========================================\n")
  cat("显著 DEG 行数:", nrow(tb_sig), "\n")
  cat("显著 DEG 独特基因数:", dplyr::n_distinct(tb_sig$gene), "\n")
  
  ## 9.11 每个 cell type DEG 统计
  deg_summary <- tb_sig %>%
    dplyr::filter(regulation %in% c("Up", "Down")) %>%
    dplyr::group_by(cluster_id, regulation) %>%
    dplyr::summarise(n_genes = dplyr::n(), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = regulation, values_from = n_genes, values_fill = 0) %>%
    dplyr::mutate(
      Up = ifelse(is.na(Up), 0, Up),
      Down = ifelse(is.na(Down), 0, Down),
      Total = Up + Down
    ) %>%
    arrange(desc(Total))
  
  cat("\n各细胞类型 DEG 统计:\n")
  print(deg_summary)
  
  ## 9.12 保存
  treatment_dir <- file.path(OUT_DIR, paste0(treatment, "_vs_REF"))
  dir.create(treatment_dir, showWarnings = FALSE, recursive = TRUE)
  
  write.csv(tb_full, file.path(treatment_dir, "DEG_full_results_muscat.csv"), row.names = FALSE)
  write.csv(tb_sig,  file.path(treatment_dir, "DEG_filtered_paper_thresholds.csv"), row.names = FALSE)
  write.csv(deg_summary, file.path(treatment_dir, "DEG_summary_muscat.csv"), row.names = FALSE)
}

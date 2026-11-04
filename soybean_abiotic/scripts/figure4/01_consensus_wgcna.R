# Source: 图4.docx, extracted text lines 1-1823
# 独立运行；CK2补全RDS；旧模块表用于可选比较
# See README.md and docs/KNOWN_LIMITATIONS.md before running.

# 完整 v5：CK固定三个表达聚类，定向命名P1/P2/P3；全细胞共识WGCNA；论文式灰色渐变UMAP及紧凑大字GO图。
# 全细胞共识WGCNA，Endodermis沿CK参考轨迹划分P1/P2/P3。
# 单个脚本即可运行；无需source旧脚本。原输入RDS不写回。
# P1/P2/P3由CK RNA PCA上的k-means(k=3)确定，使用早期参照命名，不使用module评分定义P1。
MERGED_RDS <- file.path(Sys.getenv("SOYBEAN_PROJECT_ROOT", unset = "data"), "CK2_added_v3_Results/seurat_obj_with_CK2_v3.rds")
OUT_BASE <- file.path(dirname(MERGED_RDS),'Ensemble_PI_Heat_Cold_PEG_with_CK2')
# 固定使用本次指定的原输出目录，不继承旧会话中可能残留的OUT_DIR。
OUT_DIR <- OUT_BASE
# 同目录独立文件前缀。重新运行会更新本版同名前缀的结果/缓存，不碰Nature_旧结果。
FILE_PREFIX <- 'Endo3_v5_'
NETWORK_NAME <- 'AllCells_Endo3_consensus_v5'
PIPELINE_VERSION <- 'Endo3_v5_three_clusters_publication_plots'
SCRIPT_VERSION <- PIPELINE_VERSION
CELLTYPE_COL <- 'cell_type'
SAMPLE_COL <- NULL
CONDITION_COL <- NULL
ASSAY <- 'RNA'
UMAP_REDUCTION <- 'umap'
GRAPH_REDUCTION <- NULL
HARMONY_COL <- NULL
HARMONIZE_MES <- FALSE
CONDITIONS <- c('CK','Heat','Cold','PEG')
TREATMENTS <- CONDITIONS[-1]
EXPECTED_REPLICATES <- 2L
SEED <- 20261026L
N_THREADS <- 4L
RESUME <- TRUE

# 内皮层重新聚类/发育方向。只用CK拟合参考轨迹，三个处理投射到共同参考。
ENDO_CELLTYPE <- 'Endodermis'
EARLY_ANCHOR_TYPES <- c('Meristem','QC')
ROOT_CLUSTER <- NULL # 新三个cluster中的1/2/3；默认用CK Meristem/QC表达接近性命名P1，不继承旧cluster5
ROOT_CELLS_FILE <- NULL # 可选CSV含cell列；已验证的早期内皮层条码，优先于自动锚定
MIN_ROOT_CELLS <- 10L
REFERENCE_MAX_CELLS_PER_CK <- 1200L # CK1/CK2等量抽样拟合PCA/聚类/轨迹
TRAJECTORY_N_HVG <- 2500L
TRAJECTORY_N_PCS <- 20L
PROJECTION_DISTANCE_QUANTILE <- 0.99
PROJECTION_DISTANCE_MULTIPLIER <- 1.5
MIN_CELLS <- 20L
MIN_METACELL_GROUP <- 40L
MAX_TRAIN_CELLS_PER_SAMPLE_CT <- 500L
N_HVG <- 4000L
USE_INPUT_HVG <- FALSE # 细分内皮层后重新选择HVG
MIN_DETECT_FRACTION <- 0.10
EXCLUDE_GENES <- character()
METACELL_K <- 20L
METACELL_MAX_SHARED <- 5L
TARGET_METACELLS <- 100L
MIN_METACELLS_PER_SAMPLE <- 20L
MIN_MODULE_SIZE <- 30L
SOFT_POWERS <- c(1:10,seq(12,20,2))
SOFT_R2 <- 0.80
SOFT_MIN_K <- 5
SOFT_POWER_FIXED <- 6L # 沿用旧设定，但重新检查全部8样本，失败时不会偷偷放宽
SAVE_NETWORK_RDS <- TRUE
COMPUTE_MES <- TRUE
COMPUTE_HUBS <- TRUE
SEURAT_CTRL <- 100L
SEURAT_NBIN <- 24L
UCell_MAX_RANK <- 1500L
SCORE_CHUNK <- 2000L
MIN_ABS_DELTA <- 0
FDR_CUT <- 0.05
MAX_UMAP_CELLS_PER_CONDITION <- 12000L
PSEUDOTIME_BINS <- 10L
CELLTYPE_Z_CUT <- 0.75
MAX_LABEL_CELLTYPES <- 3L

RUN_GO <- TRUE
ORGDB_PACKAGES <- c('org.Gmax.eg.db','org.G.max.eg.db')
KEYTYPE_OVERRIDE <- NULL
MIN_DB_MAPPING_RATE <- 0.20
GENE_ID_MAP_FILE <- NULL # 列Gene,DB_ID,Package,Keytype；显式跨版本映射，不猜ID
MIN_GS_SIZE <- 10L
MAX_GS_SIZE <- 500L
MIN_QUERY_GENES <- 5L
OLD_MODULE_FILE <- file.path(OUT_DIR,'EndoDev_v4_network_modules.csv') # 上一轮已完成的11个模块，含DevProgram02


out  <- function(x) file.path(OUT_DIR, paste0(FILE_PREFIX, x))
csv  <- function(x, name) utils::write.csv(x, out(name), row.names = FALSE,
                                           na = '', fileEncoding = 'UTF-8')
bind <- function(x) if (length(x)) do.call(rbind, x) else data.frame()
canon <- function(x) gsub('[[:space:]_.-]+', '', tolower(trimws(as.character(x))))

note <- function(...) {
  message(...)
  cat(format(Sys.time()), paste0(..., collapse = ''), '\n',
      file = out('run.log'), append = TRUE)
}

need <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop('缺少R包：', paste(missing, collapse = ', '),
         '\n请先安装；limma/GO.db/AnnotationDbi/UCell用BiocManager，',
         'hdWGCNA按官方文档安装；优先使用你本地的org.Gmax.eg.db；注释包缺失仅影响GO。')
  }
}

pick_ct <- function(md) {
  if (!is.null(CELLTYPE_COL)) {
    if (!CELLTYPE_COL %in% names(md)) stop('CELLTYPE_COL不存在')
    return(CELLTYPE_COL)
  }
  hits <- names(md)[canon(names(md)) %in%
                      canon(c('cell_type','celltype','celltype_final','annotation','celltype_annotated',
                              'cell_type_annotation','seurat_annotations'))]
  if (!length(hits)) stop('无法识别细胞类型列，请设置CELLTYPE_COL。')
  if (length(hits) > 1 &&
      !all(vapply(hits[-1], function(k)
        identical(as.character(md[[k]]), as.character(md[[hits[1]]])), logical(1)))) {
    stop('多个注释列不一致，请明确设置CELLTYPE_COL：', paste(hits, collapse = ', '))
  }
  hits[1]
}

condition_key <- function(x) {
  aliases <- list(
    CK   = c('CK','Control','CTRL','REF','Reference','对照'),
    Heat = c('Heat','High_Temp','HighTemp','HighTemperature','HT','高温'),
    Cold = c('Cold','Low_Temp','LowTemp','LowTemperature','LT','低温'),
    PEG  = c('PEG','PEG6000','PEG8000')
  )
  keys <- unlist(lapply(aliases, canon), use.names = FALSE)
  rep(names(aliases), lengths(aliases))[match(canon(x), keys)]
}

save_plot <- function(p, name, w = 12, h = 8) {
  ggplot2::ggsave(out(paste0(name, '.pdf')), p, width = w, height = h,
                  limitsize = FALSE, useDingbats = FALSE)
  ggplot2::ggsave(out(paste0(name, '.png')), p, width = w, height = h, dpi = 300,
                  limitsize = FALSE, bg = 'white')
}

cache <- function(name, key, calculate) {
  path <- out(paste0(name, '.rds'))
  if (RESUME && file.exists(path)) {
    z <- readRDS(path)
    if (identical(z$key, key)) {
      note('复用本版缓存：', name)
      return(z$value)
    }
  }
  z <- calculate()
  saveRDS(list(key = key, value = z), path)
  z
}


bind_rows <- bind
read_csv <- function(path) utils::read.csv(path, check.names=FALSE, stringsAsFactors=FALSE)

assess_fixed_soft_power <- function(f) {
  required <- c('Power','Signed_R2','mean.k.')
  if (!all(required %in% names(f))) stop('软阈值诊断表缺少所需列。')
  if (length(SOFT_POWER_FIXED)!=1L || !is.finite(SOFT_POWER_FIXED) || SOFT_POWER_FIXED<=0)
    stop('SOFT_POWER_FIXED必须为单个正数。')
  f$Chosen <- !is.na(f$Power) & f$Power==SOFT_POWER_FIXED
  if (sum(f$Chosen)!=1L) stop('诊断表没有唯一的固定Power行：',SOFT_POWER_FIXED)
  f$Pass_R2 <- is.finite(f$Signed_R2) & f$Signed_R2>=SOFT_R2
  f$Pass_connectivity <- is.finite(f$mean.k.) & f$mean.k.>=SOFT_MIN_K
  f$QC_pass <- f$Pass_R2 & f$Pass_connectivity
  f
}



celltype_key <- function(x) {
  z <- trimws(as.character(x))
  hit <- !is.na(z) & canon(z) %in% c('lateralroot', 'procambium')
  z[hit] <- 'Procambium'
  z
}

# 样本别名仅映射明确的重复身份。REF/Control 代表处理，不能单独认定为 CK1。
sample_key <- function(x) {
  all_samples <- c('CK1','CK2','ABA1','ABA2','eBL1','eBL2',
                   'highTemp1','highTemp2','IAA1','IAA2','lowTemp1','lowTemp2',
                   'MEJA1','MEJA2','PEG1','PEG2','SA1','SA2')
  aliases <- setNames(lapply(all_samples, function(nm) c(nm, paste0('Gm_', nm))), all_samples)
  aliases$CK1 <- c(aliases$CK1, '2228_CK', 'Gm_2228_CK')
  aliases$CK2 <- c(aliases$CK2, '2229_CK', 'Gm_2229_CK')
  aliases$highTemp1 <- c(aliases$highTemp1, 'Heat1', 'High_Temp1', 'Gm_8.37temp1')
  aliases$highTemp2 <- c(aliases$highTemp2, 'Heat2', 'High_Temp2', 'Gm_8.37temp2')
  aliases$lowTemp1 <- c(aliases$lowTemp1, 'Cold1', 'Low_Temp1', 'Gm_9.6temp1')
  aliases$lowTemp2 <- c(aliases$lowTemp2, 'Cold2', 'Low_Temp2', 'Gm_9.6temp2')
  keys <- unlist(lapply(aliases, canon), use.names = FALSE)
  values <- rep(names(aliases), lengths(aliases))
  unname(values[match(canon(x), keys)])
}

resolve_sample_design <- function(md) {
  if (!is.null(SAMPLE_COL) && !SAMPLE_COL %in% names(md)) stop('指定 SAMPLE_COL 不存在。')
  cols <- intersect(unique(c(SAMPLE_COL, 'sample_id_complete', 'sample_id', 'sample', 'orig.ident')), names(md))
  if (!length(cols)) stop('没有可识别的样本列。')
  mapped <- lapply(cols, function(nm) sample_key(md[[nm]]))
  names(mapped) <- cols
  chosen <- rep(NA_character_, nrow(md))
  origin <- rep(NA_character_, nrow(md))
  conflict <- rep(FALSE, nrow(md))
  for (nm in cols) {
    z <- mapped[[nm]]
    conflict <- conflict | (!is.na(chosen) & !is.na(z) & chosen != z)
    use <- is.na(chosen) & !is.na(z)
    chosen[use] <- z[use]; origin[use] <- nm
  }
  raw <- md[, cols, drop = FALSE]
  names(raw) <- paste0('Input__', names(raw))
  audit <- cbind(data.frame(cell = rownames(md), Sample = chosen, Source_column = origin,
                            Conflicting_sample_columns = conflict, stringsAsFactors = FALSE), raw)
  csv(audit, 'input_sample_identity_audit.csv')
  if (any(conflict)) stop('不同 metadata 列给出矛盾的重复身份，见 input_sample_identity_audit.csv。')
  
  # 以独立样本身份推导实验条件，再与已有处理列交叉核对。
  design <- data.frame(
    Sample = c('CK1','CK2','highTemp1','highTemp2','lowTemp1','lowTemp2','PEG1','PEG2'),
    Condition = rep(c('CK','Heat','Cold','PEG'), each = 2),
    Replicate = rep(1:2, 4), stringsAsFactors = FALSE)
  cond <- design$Condition[match(chosen, design$Sample)]
  if (!is.null(CONDITION_COL) && !CONDITION_COL %in% names(md)) stop('指定 CONDITION_COL 不存在。')
  condition_cols <- intersect(unique(c(CONDITION_COL, 'treatment', 'Group1', 'condition', 'Condition')), names(md))
  target_flag <- rep(FALSE, nrow(md))
  condition_conflict <- rep(FALSE, nrow(md))
  for (nm in condition_cols) {
    v <- condition_key(md[[nm]])
    target_flag <- target_flag | !is.na(v)
    condition_conflict <- condition_conflict |
      (!is.na(cond) & !is.na(v) & cond != v) |
      (!is.na(chosen) & is.na(cond) & !is.na(v))
  }
  audit$Condition <- cond
  audit$Conflicting_condition_columns <- condition_conflict
  audit$Unresolved_target_sample <- target_flag & is.na(chosen)
  audit$Selected_for_analysis <- !is.na(cond)
  csv(audit, 'input_sample_identity_audit.csv')
  inventory <- as.data.frame(table(Sample = ifelse(is.na(chosen), 'UNRESOLVED', chosen)),
                             stringsAsFactors = FALSE)
  names(inventory)[2] <- 'Cells'
  inventory$Condition <- design$Condition[match(inventory$Sample, design$Sample)]
  inventory$Selected_for_analysis <- inventory$Sample %in% design$Sample
  csv(inventory, 'input_all_samples_cell_counts.csv')
  print(inventory, row.names = FALSE)
  if (any(condition_conflict)) stop('样本身份与处理列冲突，见 input_sample_identity_audit.csv。')
  if (any(audit$Unresolved_target_sample)) {
    stop('目标条件中有细胞无法区分重复。REF/Control 不能代替 CK1/CK2；见样本身份审计。')
  }
  absent <- setdiff(design$Sample, unique(chosen))
  if (length(absent)) stop('缺少目标生物学样本：', paste(absent, collapse = ', '),
                           '。请确认读取的是含 CK2 的新文件，不要复制或拆分 CK1。')
  list(Sample = chosen, Condition = cond, design = design, audit = audit)
}

# 只在 assay 副本中 JoinLayers，不 subset/改写载入的旧 Seurat，避开旧 reduction 顺序问题。
extract_counts_readonly <- function(obj, cells, assay = ASSAY) {
  if (!assay %in% names(obj@assays)) stop('缺少 assay：', assay)
  a <- obj@assays[[assay]]
  if (inherits(a, 'Assay5')) {
    layer_names <- SeuratObject::Layers(a, search = '^counts')
    if (!length(layer_names)) stop('原始 RNA counts 不存在。')
    a <- SeuratObject::JoinLayers(a, layers = 'counts', new = 'counts')
    counts <- SeuratObject::LayerData(a, layer = 'counts')
  } else counts <- SeuratObject::GetAssayData(a, slot = 'counts')
  if (is.null(rownames(counts)) || is.null(colnames(counts)) ||
      anyDuplicated(rownames(counts)) || anyDuplicated(colnames(counts))) stop('counts 基因/细胞名缺失或重复。')
  if (!all(cells %in% colnames(counts))) stop('counts 缺少选中的细胞。')
  counts <- methods::as(counts[, cells, drop = FALSE], 'dgCMatrix')
  if (any(!is.finite(counts@x)) || any(counts@x < 0) ||
      any(abs(counts@x - round(counts@x)) > 1e-6)) stop('RNA counts 必须是非负原始整数。')
  counts
}

extract_embedding_readonly <- function(obj, reduction, cells) {
  if (!reduction %in% names(obj@reductions)) stop('缺少 reduction：', reduction)
  emb <- obj@reductions[[reduction]]@cell.embeddings
  if (is.null(rownames(emb)) || anyDuplicated(rownames(emb)) ||
      !all(cells %in% rownames(emb))) stop(reduction, ' 未覆盖选中细胞，或存在重复细胞名。')
  emb <- emb[cells, , drop = FALSE]
  if (any(!is.finite(emb))) stop(reduction, ' 中存在非有限数值。')
  emb
}


read_data <- function() {
  obj <- readRDS(MERGED_RDS)
  if (!inherits(obj, 'Seurat')) stop('输入不是 Seurat 对象。')
  md <- obj@meta.data
  if (is.null(rownames(md)) || anyDuplicated(rownames(md))) stop('metadata 细胞名缺失或重复。')
  ct <- pick_ct(md)
  ids <- resolve_sample_design(md)
  m <- data.frame(
    cell = rownames(md), Celltype = celltype_key(md[[ct]]),
    Condition = ids$Condition, Sample = ids$Sample, stringsAsFactors = FALSE)
  selected <- !is.na(m$Condition)
  annotated <- !is.na(m$Celltype) & nzchar(m$Celltype) &
    !tolower(m$Celltype) %in% c('unknown','unassigned','unclassified','undefined','na')
  keep <- selected & annotated
  selection <- data.frame(cell = m$cell, Sample = m$Sample,
                          Condition = m$Condition, Celltype = m$Celltype,
                          Included = keep,
                          Reason = ifelse(!selected, 'Outside_CK_Heat_Cold_PEG',
                                          ifelse(!annotated, 'Missing_or_unknown_annotation', 'Included')))
  csv(selection, 'input_cell_selection_audit.csv')
  original_ct <- as.character(md[[ct]])
  csv(as.data.frame(table(Original = original_ct, Standardized = m$Celltype)),
      'celltype_name_mapping.csv')
  m <- m[keep, , drop = FALSE]
  if (!nrow(m) || anyNA(m$Sample)) stop('没有可用目标细胞或样本ID缺失。')
  pairs <- unique(m[c('Sample','Condition')])
  if (anyDuplicated(pairs$Sample)) stop('同一样本跨多个条件。')
  if (!setequal(pairs$Sample, ids$design$Sample)) stop('过滤缺失注释后丢失了某个目标样本。')
  ns <- table(factor(pairs$Condition, CONDITIONS))
  print(ns)
  if (any(ns != EXPECTED_REPLICATES)) stop('各条件必须各有两个独立生物学重复。')
  if (any(grepl('^[0-9]+$', m$Celltype))) stop('细胞类型列含数字cluster，请使用 cell_type。')
  pairs <- pairs[match(ids$design$Sample, pairs$Sample), , drop = FALSE]
  pairs$Replicate <- ids$design$Replicate
  pairs$Cells <- as.integer(table(m$Sample)[pairs$Sample])
  csv(pairs, 'biological_samples.csv')
  print(pairs, row.names = FALSE)
  csv(data.frame(Input = MERGED_RDS, Total_input_cells = nrow(md),
                 Identified_input_samples = length(unique(na.omit(ids$Sample))),
                 Analysis_cells = nrow(m), Analysis_samples = nrow(pairs)),
      'input_object_summary.csv')
  note('输入细胞 ', nrow(md), '；分析仅保留 CK/Heat/Cold/PEG 的 ', nrow(m), ' 个细胞、8 个独立样本。')
  
  counts <- extract_counts_readonly(obj, m$cell)
  m$nUMI <- as.numeric(Matrix::colSums(counts))
  m$nFeature <- as.numeric(Matrix::colSums(counts > 0))
  if (any(m$nUMI <= 0)) stop('存在零文库细胞。')
  xy <- extract_embedding_readonly(obj, UMAP_REDUCTION, m$cell)
  if (ncol(xy) < 2L) stop('指定 UMAP 少于两维。')
  m$UMAP1 <- xy[, 1]; m$UMAP2 <- xy[, 2]
  rownames(m) <- m$cell
  if (!is.null(HARMONY_COL)) {
    if (!HARMONY_COL %in% names(md)) stop('HARMONY_COL不存在。')
    m$Batch <- as.character(md[m$cell, HARMONY_COL])
  }
  csv(m, 'cell_metadata.csv')
  qc <- expand.grid(Sample = pairs$Sample, Celltype = sort(unique(m$Celltype)),
                    stringsAsFactors = FALSE)
  qc$Condition <- pairs$Condition[match(qc$Sample, pairs$Sample)]
  tab <- table(m$Sample, m$Celltype)
  qc$n_cells <- as.integer(tab[cbind(match(qc$Sample, rownames(tab)),
                                     match(qc$Celltype, colnames(tab)))])
  qc$Score_test_eligible <- qc$n_cells >= MIN_CELLS
  qc$Metacell_eligible <- qc$n_cells >= MIN_METACELL_GROUP
  csv(qc, 'sample_celltype_QC.csv')
  
  candidates <- if (is.null(GRAPH_REDUCTION)) c('integrated.rpca','rpca','harmony') else GRAPH_REDUCTION
  hit <- intersect(candidates, names(obj@reductions))
  # 部分旧降维未覆盖新 CK2 时，自动模式改用 RNA PCA；不使用 UMAP 构网。
  complete <- vapply(hit, function(nm) all(m$cell %in% rownames(obj@reductions[[nm]]@cell.embeddings)), logical(1))
  hit <- hit[complete]
  if (!is.null(GRAPH_REDUCTION) && !length(hit)) stop('指定 GRAPH_REDUCTION 不存在或未覆盖全部目标细胞。')
  graph <- if (length(hit)) extract_embedding_readonly(obj, hit[1], m$cell) else NULL
  list(counts = counts, md = m, qc = qc,
       input_hvg = SeuratObject::VariableFeatures(obj@assays[[ASSAY]]),
       graph = graph, graph_name = if (length(hit)) hit[1] else 'RNA_PCA_new')
}


# -----------------------------------------------------------------------------
# 6. 共识网络构建
# -----------------------------------------------------------------------------
# 验证 metacell 表达矩阵；平均型 metacell 允许小数，不调用原始counts的整数校验。
validate_metacell_matrix <- function(x, cells, tag) {
  if (length(dim(x)) != 2L || !nrow(x) || !ncol(x)) stop(tag, '：metacell矩阵为空。')
  if (is.null(rownames(x)) || is.null(colnames(x)) ||
      anyNA(rownames(x)) || anyNA(colnames(x)) ||
      anyDuplicated(rownames(x)) || anyDuplicated(colnames(x))) stop(tag, '：矩阵行列名缺失或重复。')
  if (!setequal(colnames(x), cells)) stop(tag, '：表达矩阵与metadata的metacell集合不一致。')
  x <- methods::as(x[, cells, drop = FALSE], 'dgCMatrix')
  if (any(!is.finite(x@x)) || any(x@x < 0)) stop(tag, '：表达值包含非有限值或负数。')
  if (any(Matrix::colSums(x) <= 0)) stop(tag, '：存在零表达metacell。')
  x
}

# 取出 hdWGCNA 内部 metacell，合并 RNA counts layers，再转为经典 Assay。
# 保留 metacell对象的metadata、misc和分组，不重新生成metacell，不平均/取整已聚合矩阵。
prepare_metacell_rna <- function(train, nm, expected_samples) {
  mc <- hdWGCNA::GetMetacellObject(train, wgcna_name = nm)
  if (!inherits(mc, 'Seurat')) stop('hdWGCNA未返回有效的metacell Seurat对象。')
  if (!'RNA' %in% names(mc@assays)) stop('metacell对象没有RNA assay。')
  md_before <- mc[[]]
  cells <- rownames(md_before)
  needed <- c('ers_sample','ers_celltype','ers_pool')
  if (!all(needed %in% names(md_before))) stop('metacell对象丢失样本或细胞类型分组。')
  if (!setequal(as.character(md_before$ers_sample), expected_samples)) stop('metacell对象未包含完整8个样本。')
  a <- mc@assays[['RNA']]
  class_before <- class(a)[1]
  layers_before <- if (inherits(a, 'Assay5')) SeuratObject::Layers(a) else c('counts','data','scale.data')
  csv(data.frame(Stage = 'Before_conversion', RNA_class = class_before, Layer = layers_before),
      'metacell_layers_before.csv')
  if (inherits(a, 'Assay5')) {
    count_layers <- SeuratObject::Layers(a, search = '^counts')
    if (!length(count_layers)) stop('metacell对象缺少counts层。')
    # 只对独立 assay 操作；旧 data.* 不复用，下一步统一重新做log-normalization。
    a <- SeuratObject::JoinLayers(a, layers = 'counts', new = 'counts')
    counts <- SeuratObject::LayerData(a, layer = 'counts')
  } else counts <- SeuratObject::GetAssayData(a, slot = 'counts')
  counts <- validate_metacell_matrix(counts, cells, 'metacell counts')
  classic <- SeuratObject::CreateAssayObject(counts = counts, min.cells = 0, min.features = 0)
  mc[['RNA']] <- classic
  SeuratObject::DefaultAssay(mc) <- 'RNA'
  if (inherits(mc[['RNA']], 'Assay5')) stop('metacell RNA未转换为经典Assay。')
  actual <- SeuratObject::GetAssayData(mc[['RNA']], slot = 'counts')
  if (!identical(rownames(actual), rownames(counts)) ||
      !setequal(colnames(actual), colnames(counts)) ||
      Matrix::nnzero(actual[, cells, drop = FALSE] - counts) != 0) {
    stop('metacell转换过程中表达值或基因名发生变化。')
  }
  keep_cols <- intersect(c(needed, 'cells_merged'), names(md_before))
  md_after <- mc[[]]
  for (co in keep_cols) {
    if (!identical(as.character(md_before[[co]]), as.character(md_after[cells, co]))) {
      stop('metacell转换过程中分组/成员信息发生变化：', co)
    }
  }
  methods::validObject(mc)
  csv(data.frame(RNA_class_before = class_before, RNA_class_after = class(mc[['RNA']])[1],
                 Metacells = ncol(counts), Genes = nrow(counts),
                 Values_unchanged = TRUE,
                 Fractional_values_preserved = any(abs(counts@x - round(counts@x)) > 1e-6)),
      'metacell_RNA_conversion.csv')
  train <- hdWGCNA::SetMetacellObject(train, metacell_obj = mc, wgcna_name = nm)
  note('内部metacell RNA已整理为经典Assay；保留全部 ', ncol(counts), ' 个metacell及其数值。')
  train
}

# 在 SetMultiExpr 前复现它所需的逐样本 GetAssayData 访问，尽早发现层/样本缺失问题。
check_metacell_expression <- function(train, nm, expected_samples) {
  mc <- hdWGCNA::GetMetacellObject(train, wgcna_name = nm)
  if (inherits(mc[['RNA']], 'Assay5')) stop('归一化后metacell RNA重新变成Assay5。')
  md <- mc[[]]
  checks <- lapply(expected_samples, function(s) {
    cells <- rownames(md)[as.character(md$ers_sample) == s]
    if (!length(cells)) stop('metacell缺少样本：', s)
    one <- subset(mc, cells = cells)
    expr <- SeuratObject::GetAssayData(one, assay = 'RNA', slot = 'data')
    expr <- validate_metacell_matrix(expr, cells, paste('normalized metacells', s))
    data.frame(Sample = s, Metacells = ncol(expr), Genes = nrow(expr),
               RNA_class = class(one[['RNA']])[1], Normalized_data_readable = TRUE)
  })
  csv(bind(checks), 'metacell_expression_before_SetMultiExpr.csv')
  invisible(TRUE)
}

build_network <- function(dat) {
  note('建立跨样本共识网络；首次运行可能需要较长时间。')
  m <- dat$md
  
  # 6.1 每样本 × 细胞类型抽样，构建平衡训练集
  groups <- split(seq_len(nrow(m)), interaction(m$Sample, m$Celltype, drop = TRUE))
  set.seed(SEED)
  ix <- unlist(lapply(groups, function(i) {
    if (length(i) < MIN_METACELL_GROUP || any(!m$Network_eligible[i])) integer()
    else sample(i, min(length(i), MAX_TRAIN_CELLS_PER_SAMPLE_CT))
  }), use.names = FALSE)
  if (!length(ix)) stop('没有足够细胞创建metacell。')
  
  mt <- m[ix, ]
  if (!setequal(mt$Sample, m$Sample)) stop('部分样本没有任何可构网的细胞类型。')
  
  train <- SeuratObject::CreateSeuratObject(
    counts = SeuratObject::CreateAssayObject(
      counts = dat$counts[, mt$cell, drop = FALSE]),
    meta.data = mt)
  train$ers_celltype <- mt$Celltype
  train$ers_sample   <- mt$Sample
  train$ers_pool     <- 'All'
  train <- Seurat::NormalizeData(train, verbose = FALSE)
  train <- Seurat::FindVariableFeatures(train, nfeatures = N_HVG, verbose = FALSE)
  
  # 6.2 基因过滤：至少一种细胞类型中检出率达到阈值（原脚本使用 OR）。
  eligible <- rep(FALSE, nrow(dat$counts))
  names(eligible) <- rownames(dat$counts)
  for (ct in unique(mt$Celltype)) {
    cs <- mt$cell[mt$Celltype == ct]
    eligible <- eligible |
      Matrix::rowMeans(dat$counts[, cs, drop = FALSE] > 0) >= MIN_DETECT_FRACTION
  }
  
  hvgs  <- if (USE_INPUT_HVG && length(dat$input_hvg)) dat$input_hvg
  else SeuratObject::VariableFeatures(train)
  genes <- setdiff(intersect(hvgs, names(eligible)[eligible]), EXCLUDE_GENES)
  if (length(genes) < 200) stop('可用网络基因少于200，检查表达/过滤。')
  csv(data.frame(Gene = genes), 'network_input_genes.csv')
  
  # 6.3 PCA 或外部降维
  train <- Seurat::ScaleData(train, features = genes, verbose = FALSE)
  npcs  <- min(30L, length(genes) - 1L, ncol(train) - 1L)
  train <- Seurat::RunPCA(train, features = genes, npcs = npcs,
                          verbose = FALSE, seed.use = SEED)
  
  nm <- NETWORK_NAME
  red <- 'pca'; ndims <- npcs
  if (!is.null(dat$graph)) {
    emb <- dat$graph[mt$cell, , drop = FALSE]
    colnames(emb) <- paste0('networkInput_', seq_len(ncol(emb)))
    train[['networkInput']] <- SeuratObject::CreateDimReducObject(
      embeddings = emb, key = 'networkInput_', assay = 'RNA')
    red <- 'networkInput'; ndims <- min(30L, ncol(emb))
  }
  csv(data.frame(
    Reduction = dat$graph_name,
    Dimensions = ndims,
    Input_HVG_used = USE_INPUT_HVG && length(dat$input_hvg) > 0),
    'network_reduction.csv')
  
  # 6.4 hdWGCNA setup + metacell
  train <- hdWGCNA::SetupForWGCNA(train, wgcna_name = nm,
                                  gene_select = 'custom', gene_list = genes)
  train <- hdWGCNA::MetacellsByGroups(
    train,
    group.by = c('ers_celltype','ers_sample','ers_pool'),
    ident.group = 'ers_celltype',
    k = METACELL_K, max_shared = METACELL_MAX_SHARED,
    min_cells = MIN_METACELL_GROUP, target_metacells = TARGET_METACELLS,
    reduction = red, dims = seq_len(ndims),
    assay = 'RNA', slot = 'counts', wgcna_name = nm)
  # 关键修复：必须整理内部 metacell RNA，单独整理 train[['RNA']] 不足以修复报错。
  train <- prepare_metacell_rna(train, nm, sort(unique(mt$Sample)))
  train <- hdWGCNA::NormalizeMetacells(train, wgcna_name = nm, assay = 'RNA', verbose = FALSE)
  check_metacell_expression(train, nm, sort(unique(mt$Sample)))
  
  # 6.5 构造 multiExpr（每个样本一个数据集）
  train <- hdWGCNA::SetMultiExpr(
    train,
    group_name = 'All', group.by = 'ers_pool',
    multi.group.by = 'ers_sample',
    multi_groups = sort(unique(mt$Sample)),
    assay = 'RNA', slot = 'data', wgcna_name = nm)
  
  multi <- hdWGCNA::GetMultiExpr(train, wgcna_name = nm)
  if (length(multi) != length(unique(mt$Sample))) stop('共识表达矩阵缺少样本。')
  if (is.null(names(multi))) stop('multiExpr没有样本名称，无法审计顺序。')
  if (!setequal(names(multi), unique(mt$Sample)))
    stop('multiExpr样本名不匹配，请检查hdWGCNA版本。')
  
  bg <- colnames(multi[[1]]$data)
  if (any(vapply(multi, function(x) !identical(colnames(x$data), bg), logical(1))))
    stop('样本间网络基因顺序不一致。')
  
  qc <- data.frame(
    Sample    = names(multi),
    Metacells = vapply(multi, function(x) nrow(x$data), integer(1)),
    Genes     = length(bg))
  csv(qc, 'metacell_QC.csv')
  if (any(qc$Metacells < MIN_METACELLS_PER_SAMPLE))
    stop('某样本metacell太少；见ERS_metacell_QC.csv。')
  
  # 6.6 逐样本诊断；共同固定Power=6
  fits <- list(); powers <- numeric(length(multi))
  for (i in seq_along(multi)) {
    x <- multi[[i]]$data
    if (any(!is.finite(as.matrix(x))) || any(apply(x, 2, stats::sd) == 0))
      stop('共识矩阵包含非有限值/零方差基因；请检查SetMultiExpr过滤。')
    
    f <- WGCNA::pickSoftThreshold(x, powerVector = SOFT_POWERS,
                                  networkType = 'signed', verbose = 0)$fitIndices
    f$Sample    <- names(multi)[i]
    f$Signed_R2 <- -sign(f$slope) * f$SFT.R.sq
    
    f <- assess_fixed_soft_power(f)
    powers[i] <- SOFT_POWER_FIXED
    fits[[i]]  <- f
  }
  fits <- bind(fits)
  csv(fits, 'soft_power_diagnostics.csv')
  
  p <- ggplot2::ggplot(fits, ggplot2::aes(Power, Signed_R2, color = Sample)) +
    ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::geom_hline(yintercept = SOFT_R2, linetype = 2) +
    ggplot2::geom_vline(xintercept = SOFT_POWER_FIXED, linetype = 3) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = 'Consensus network: soft-power diagnostics')
  save_plot(p, 'QC_soft_power', 10, 6)
  selected <- fits[fits$Chosen,,drop=FALSE]
  csv(selected,'soft_power_selected.csv')
  p_k <- ggplot2::ggplot(fits,ggplot2::aes(Power,mean.k.,color=Sample)) +
    ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::geom_vline(xintercept=SOFT_POWER_FIXED,linetype=3) +
    ggplot2::geom_hline(yintercept=SOFT_MIN_K,linetype=2) +
    ggplot2::theme_bw() + ggplot2::labs(y='Mean connectivity',
                                        title='Consensus network: connectivity diagnostics')
  save_plot(p_k,'QC_soft_power_connectivity',10,6)
  if (any(!selected$QC_pass))
    stop('固定Power=',SOFT_POWER_FIXED,'未通过R2/连接度检查：',
         paste(selected$Sample[!selected$QC_pass],collapse=', '),
         '；检查Nature_soft_power_selected.csv。没有自动更换阈值。')
  note('全部样本统一Power=',SOFT_POWER_FIXED,'；R2与平均连接度检查通过。')
  
  # 6.7 共识网络 + 模块
  train <- hdWGCNA::ConstructNetwork(
    train, consensus = TRUE, soft_power = unname(powers),
    tom_outdir = out('TOM'), tom_name = 'ERS_consensus',
    overwrite_tom = TRUE, wgcna_name = nm,
    networkType = 'signed', TOMType = 'signed',
    minModuleSize = MIN_MODULE_SIZE, deepSplit = 2,
    maxBlockSize = max(5000L, length(bg) + 1L), randomSeed = SEED)
  
  mods <- hdWGCNA::GetModules(train, wgcna_name = nm)
  if (!all(c('gene_name','module') %in% names(mods)))
    stop('模块表缺少gene_name/module。')
  
  ids <- sort(unique(as.character(mods$module[
    !tolower(as.character(mods$module)) %in% c('grey','gray','0') &
      !is.na(mods$module)])))
  if (!length(ids)) stop('没有非grey模块；不能强行生成共享程序。')
  
  map <- setNames(sprintf('StageProgram%02d', seq_along(ids)), ids)
  mods$original_module <- as.character(mods$module)
  mods$Program <- unname(map[as.character(mods$module)])
  mods$module  <- factor(ifelse(is.na(mods$Program), 'grey', mods$Program),
                         levels = c(unname(map), 'grey'))
  train <- hdWGCNA::SetModules(train, mods, wgcna_name = nm)
  if (anyDuplicated(mods$gene_name)) stop('网络模块表有重复基因。')
  csv(mods, 'network_modules.csv')
  if (SAVE_NETWORK_RDS) saveRDS(train, out('consensus_network_object.rds'))
  
  list(modules = mods, background = bg, soft = fits, object = train)
}


align_ucell_scores <- function(sc, cells, programs) {
  x <- as.matrix(sc)
  valid_axis <- function(v, expected) {
    !is.null(v) && !anyNA(v) && !anyDuplicated(v) &&
      length(v) == length(expected) && setequal(v, expected)
  }
  if (!valid_axis(rownames(x), cells) && valid_axis(colnames(x), cells))
    x <- t(x)
  if (!valid_axis(rownames(x), cells)) {
    stop('UCell细胞名不匹配；实际维度=', paste(dim(x), collapse='x'),
         '；前几个行名=', paste(head(rownames(x)), collapse=', '),
         '；期望细胞=', paste(head(cells), collapse=', '))
  }
  expected <- paste0(programs, '_UCell')
  if (valid_axis(colnames(x), expected)) {
    x <- x[cells, expected, drop=FALSE]
  } else if (valid_axis(colnames(x), programs)) {
    x <- x[cells, programs, drop=FALSE]
  } else {
    stop('UCell模块名不匹配；实际=', paste(colnames(x), collapse=', '),
         '；期望=', paste(expected, collapse=', '))
  }
  if (!is.numeric(x) || any(!is.finite(x))) stop('UCell存在非数值或非有限值。')
  colnames(x) <- programs
  x
}

score_modules <- function(dat, net) {
  z <- net$modules[!is.na(net$modules$Program), , drop=FALSE]
  sets <- lapply(split(as.character(z$gene_name), as.character(z$Program)), unique)
  if (!length(sets) || any(lengths(sets) == 0L)) stop('没有可评分模块。')
  counts <- dat$counts
  if (is.null(colnames(counts)) || anyDuplicated(colnames(counts)))
    stop('RNA矩阵缺少唯一细胞名。')
  if (any(!unlist(sets, use.names=FALSE) %in% rownames(counts)))
    stop('模块基因不在RNA矩阵。')
  chunks <- split(seq_len(ncol(counts)), ceiling(seq_len(ncol(counts))/SCORE_CHUNK))
  ans <- matrix(NA_real_, nrow=ncol(counts), ncol=length(sets),
                dimnames=list(colnames(counts), names(sets)))
  for (i in seq_along(chunks)) {
    ix <- chunks[[i]]
    note('全细胞UCell：', i, '/', length(chunks))
    sc <- UCell::ScoreSignatures_UCell(
      matrix=counts[, ix, drop=FALSE], features=sets,
      maxRank=min(max(UCell_MAX_RANK, max(lengths(sets))+1L), nrow(counts)),
      name='_UCell', ncores=1)
    if (i == 1L) {
      note('UCell输出类型：', paste(class(sc), collapse='/'),
           '；维度：', paste(dim(sc), collapse=' x '))
    }
    ans[ix, ] <- align_ucell_scores(sc, colnames(counts)[ix], names(sets))
  }
  if (any(!is.finite(ans))) stop('UCell存在非有限值。')
  ans
}


aggregate_scores <- function(scores, md) {
  if (anyDuplicated(md$cell) || anyDuplicated(rownames(scores)) ||
      !all(md$cell %in% rownames(scores))) stop('评分聚合：细胞ID不一致。')
  scores <- scores[md$cell, , drop = FALSE]
  ix <- split(seq_len(nrow(md)), interaction(md$Sample, md$Celltype, drop = TRUE))
  bind(lapply(ix, function(i) data.frame(
    Sample = md$Sample[i[1]],
    Condition = md$Condition[i[1]],
    Celltype = md$Celltype[i[1]],
    Program = colnames(scores),
    Score = colMeans(scores[i, , drop = FALSE]),
    n_cells = length(i),
    mean_logUMI = mean(log1p(md$nUMI[i])),
    mean_nFeature = mean(md$nFeature[i]))))
}

# -----------------------------------------------------------------------------
# 8. 样本层面检验 + 共享程序汇总
# -----------------------------------------------------------------------------
test_scores <- function(s, adjust_depth = FALSE) {
  res <- list(); status <- list()
  for (ct in unique(s$Celltype)) {
    d <- s[s$Celltype == ct & s$n_cells >= MIN_CELLS, ]
    for (tr in TREATMENTS) {
      z  <- d[d$Condition %in% c('CK', tr), ]
      sm <- unique(z[c('Sample','Condition','mean_logUMI')])
      ns <- table(factor(sm$Condition, c('CK', tr)))
      
      reason <- if (any(ns < 2)) 'Too_few_biological_replicates' else 'OK'
      
      if (reason == 'OK') {
        fit_result <- tryCatch({
          sm$treated <- as.integer(sm$Condition == tr)
          design <- cbind(Intercept = 1, Treatment = sm$treated)
          if (adjust_depth) {
            depth <- sm$mean_logUMI - mean(sm$mean_logUMI)
            if (stats::sd(depth) > 1e-8) design <- cbind(design, logUMI = depth)
          }
          if (qr(design)$rank < ncol(design) || nrow(design) - ncol(design) < 1)
            stop('Rank_deficient_or_no_residual_df')
          
          mods <- sort(unique(z$Program))
          mat  <- matrix(NA_real_, length(mods), nrow(sm),
                         dimnames = list(mods, sm$Sample))
          if (anyDuplicated(z[c('Program','Sample')]))
            stop('Duplicate_sample_program')
          mat[cbind(match(z$Program, mods), match(z$Sample, sm$Sample))] <- z$Score
          if (anyNA(mat)) stop('Missing_sample_program_score')
          
          f  <- limma::eBayes(limma::lmFit(mat, design), trend = FALSE)
          tt <- limma::topTable(f, coef = 'Treatment', number = Inf, sort.by = 'none')
          tt <- tt[mods, , drop = FALSE]
          data.frame(
            Celltype = ct, Treatment = tr, Program = mods,
            Delta = tt$logFC, P_value = tt$P.Value,
            N_CK = unname(ns[1]), N_treatment = unname(ns[2]),
            Depth_adjusted = adjust_depth)
        }, error = function(e) e)
        
        if (inherits(fit_result, 'error')) reason <- conditionMessage(fit_result)
        else res[[length(res) + 1]] <- fit_result
      }
      
      status[[length(status) + 1]] <- data.frame(
        Celltype = ct, Treatment = tr,
        Depth_adjusted = adjust_depth, Status = reason)
    }
  }
  
  result <- bind(res)
  if (nrow(result)) {
    result$FDR <- p.adjust(result$P_value, 'BH')
    result$Significant <- result$FDR <= FDR_CUT & abs(result$Delta) > MIN_ABS_DELTA
  }
  list(results = result, status = bind(status))
}


make_full_object <- function(dat, net) {
  x <- SeuratObject::CreateSeuratObject(
    counts = SeuratObject::CreateAssayObject(counts = dat$counts),
    meta.data = dat$md)
  x <- Seurat::NormalizeData(x, verbose = FALSE)
  
  x$ers_pool     <- 'All'
  x$ers_sample   <- dat$md$Sample
  x$ers_celltype <- dat$md$Celltype
  
  xy <- as.matrix(dat$md[, c('UMAP1','UMAP2')])
  rownames(xy) <- dat$md$cell; colnames(xy) <- c('UMAP_1', 'UMAP_2')
  x[['umap']] <- SeuratObject::CreateDimReducObject(
    embeddings = xy, key = 'UMAP_', assay = 'RNA')
  
  x@misc[[NETWORK_NAME]] <- net$object@misc[[NETWORK_NAME]]
  hdWGCNA::SetActiveWGCNA(x, wgcna_name = NETWORK_NAME)
}

get_logdata <- function(x) {
  a <- x@assays[['RNA']]
  if (inherits(a, 'Assay5')) SeuratObject::LayerData(a, layer = 'data')
  else SeuratObject::GetAssayData(a, slot = 'data')
}

checked_scores <- function(x, cells, programs, label) {
  x <- as.matrix(x)
  if (is.null(rownames(x)) || !all(cells %in% rownames(x)) ||
      !all(programs %in% colnames(x)))
    stop(label, '的细胞/模块名不匹配。')
  x <- x[cells, programs, drop = FALSE]
  if (any(!is.finite(x))) stop(label, '含非有限值。')
  x
}

all_gene_seurat_scores <- function(x) {
  mods  <- hdWGCNA::GetModules(x, wgcna_name = NETWORK_NAME)
  progs <- levels(mods$module)
  progs <- progs[progs != 'grey']
  if (!length(progs)) stop('没有非grey模块。')
  set.seed(SEED)
  
  if (length(progs) > 1) {
    x <- hdWGCNA::ModuleExprScore(
      x, n_genes = 'all', method = 'Seurat', wgcna_name = NETWORK_NAME,
      assay = 'RNA', nbin = SEURAT_NBIN, ctrl = SEURAT_CTRL,
      seed = SEED, name = 'NatureScore')
  } else {
    # 单模块避免部分 hdWGCNA 版本 drop=TRUE 导致列名丢失
    tmp <- Seurat::AddModuleScore(
      x,
      features = list(mods$gene_name[mods$module == progs[1]]),
      assay = 'RNA', nbin = SEURAT_NBIN, ctrl = SEURAT_CTRL,
      seed = SEED, name = 'NatureScore')
    one <- tmp[[]][, 'NatureScore1', drop = FALSE]
    colnames(one) <- progs
    x <- hdWGCNA::SetModuleScores(x, one, wgcna_name = NETWORK_NAME)
  }
  
  sc <- checked_scores(hdWGCNA::GetModuleScores(x, wgcna_name = NETWORK_NAME),
                       colnames(x), progs, 'Seurat score')
  
  # 明确对 log-normalized RNA 求均值，规避 AvgModuleExpr 版本差异
  expr <- get_logdata(x)
  av <- vapply(progs, function(p) {
    Matrix::colMeans(expr[mods$gene_name[mods$module == p], , drop = FALSE])
  }, numeric(ncol(x)))
  rownames(av) <- colnames(x); colnames(av) <- progs
  x <- hdWGCNA::SetAvgModuleExpr(x, av, wgcna_name = NETWORK_NAME)
  
  list(object = x, seurat = sc, average = av)
}

module_color <- function(program) {
  palette <- c('#9E1B32','#176C70','#503064','#9A336A','#263A7A','#21714A',
               '#635475','#AA7830','#1D6070','#B03D30','#78923B','#C68C39')
  n <- suppressWarnings(as.integer(gsub('[^0-9]','',program)))
  if(is.na(n)) return('#5B6674')
  palette[(n-1L) %% length(palette)+1L]
}
module_limits <- function(value) {
  value <- value[is.finite(value)]
  if(!length(value)) return(c(0,1))
  lim <- as.numeric(stats::quantile(value,c(0.01,0.99),names=FALSE))
  if(diff(lim)<1e-10) lim <- range(value)
  if(diff(lim)<1e-10) lim <- lim+c(-.5,.5)
  lim
}
plot_atlas <- function(mat,md,prefix,label,labels=NULL) {
  mat <- as.matrix(mat)
  programs <- setdiff(colnames(mat),c('grey','MEgrey'))
  if(!length(programs)) return(invisible(NULL))
  set.seed(SEED)
  ix <- unlist(lapply(split(seq_len(nrow(md)),md$Condition),function(i)
    i[sample.int(length(i),min(length(i),MAX_UMAP_CELLS_PER_CONDITION))]),use.names=FALSE)
  q <- md[ix,,drop=FALSE]; q$Condition <- factor(q$Condition,CONDITIONS)
  # 色阶从该评分矩阵的全部细胞确定，所有处理及内皮层放大图共用同一模块色阶。
  limits <- setNames(lapply(programs,function(p) module_limits(mat[,p])),programs)
  csv(bind(lapply(programs,function(p) data.frame(Program=p,Low=limits[[p]][1],High=limits[[p]][2],
                                                  Low_color='#E2E2E2',High_color=module_color(p),Scale='Pooled cells 1%-99%; display clipping only'))),
      paste0(prefix,'_color_scales.csv'))
  title_for <- function(p) {
    lab <- if(is.null(labels)) character() else labels$Label[match(p,labels$Program)]
    text <- if(!length(lab) || is.na(lab))p else paste(p,lab,sep=' | ')
    paste(strwrap(text,width=100),collapse='\n')
  }
  panel <- function(p,reference_only=FALSE) {
    d <- q
    if(reference_only) d <- d[d$Condition=='CK',,drop=FALSE]
    d$Value <- pmin(limits[[p]][2],pmax(limits[[p]][1],mat[d$cell,p]))
    d <- d[order(d$Value),,drop=FALSE]
    g <- ggplot2::ggplot(d,ggplot2::aes(UMAP1,UMAP2,color=Value))+
      ggplot2::geom_point(size=if(reference_only).2 else .16)+ggplot2::coord_equal()+
      ggplot2::scale_color_gradient(low='#E2E2E2',high=module_color(p),limits=limits[[p]])+
      ggplot2::theme_void(base_size=12)+
      ggplot2::labs(title=if(reference_only)p else title_for(p),color=label)+
      ggplot2::theme(plot.title=ggplot2::element_text(size=12,face='bold'),
                     strip.text=ggplot2::element_text(size=12),legend.title=ggplot2::element_text(size=10),
                     legend.text=ggplot2::element_text(size=9),plot.margin=grid::unit(c(3,3,3,3),'pt'))
    if(!reference_only) g <- g+ggplot2::facet_wrap(~Condition,nrow=1)
    g
  }
  pages <- split(programs,ceiling(seq_along(programs)/4))
  for(i in seq_along(pages)) {
    plots <- lapply(pages[[i]],panel)
    save_plot(patchwork::wrap_plots(plots,ncol=1),paste0(prefix,'_page',i),14,2.8*length(plots))
  }
  if(any(q$Condition=='CK')) {
    pages <- split(programs,ceiling(seq_along(programs)/6))
    for(i in seq_along(pages)) save_plot(patchwork::wrap_plots(lapply(pages[[i]],panel,reference_only=TRUE),ncol=3),
                                         paste0(prefix,'_CK_overview_page',i),11,3.2*ceiling(length(pages[[i]])/3))
  }
}
plot_go_compact <- function(go) {
  if(is.null(go) || !nrow(go$significant)) return(invisible(NULL))
  q <- bind(lapply(split(go$significant,go$significant$Program),function(a)
    head(a[order(a$FDR_global,-a$Count),],5)))
  csv(q,'GO_BP_terms_plotted.csv')
  panels <- lapply(split(q,q$Program),function(a) {
    a <- a[order(a$GeneRatio),,drop=FALSE]
    a$Term <- vapply(paste(a$Description,a$GO,sep=' | '),function(s)
      paste(strwrap(s,width=48),collapse='\n'),character(1))
    a$Term <- factor(a$Term,levels=unique(a$Term))
    ggplot2::ggplot(a,ggplot2::aes(GeneRatio,Term,size=Count,color=-log10(pmax(FDR_global,1e-300))))+
      ggplot2::geom_point()+ggplot2::scale_size_continuous(range=c(2,5))+
      ggplot2::scale_color_gradient(low='#F6C879',high='#A81930')+
      ggplot2::theme_bw(base_size=13)+
      ggplot2::theme(panel.grid.minor=ggplot2::element_blank(),
                     axis.text.y=ggplot2::element_text(size=12,color='black'),axis.text.x=ggplot2::element_text(size=11),
                     plot.title=ggplot2::element_text(size=14,face='bold'),legend.title=ggplot2::element_text(size=10),
                     legend.text=ggplot2::element_text(size=10),plot.margin=grid::unit(c(5,5,5,5),'pt'))+
      ggplot2::labs(title=paste(a$Program[1],'GO biological process'),x='Gene ratio',y=NULL,
                    color='-log10(FDR)',size='Genes')
  })
  for(p in names(panels)) save_plot(panels[[p]],paste0('Figure14_GO_',p),9,4.1)
  pages <- split(names(panels),ceiling(seq_along(panels)/2))
  for(i in seq_along(pages)) save_plot(patchwork::wrap_plots(panels[pages[[i]]],ncol=1),
                                       paste0('Figure14_module_GO_compact_page',i),9,4.1*length(pages[[i]]))
}
p1_module_evidence <- function(samples,prefix) {
  d <- samples[samples$Celltype %in% paste0('Endodermis_P',1:3) & samples$n_cells>=MIN_CELLS,,drop=FALSE]
  details <- bind(lapply(split(d,interaction(d$Program,d$Sample,drop=TRUE)),function(a) {
    if(!setequal(a$Celltype,paste0('Endodermis_P',1:3)) || nrow(a)!=3) return(NULL)
    v <- setNames(a$Score,a$Celltype)
    data.frame(Program=a$Program[1],Sample=a$Sample[1],Condition=a$Condition[1],
               P1=v[['Endodermis_P1']],P2=v[['Endodermis_P2']],P3=v[['Endodermis_P3']],
               Delta_P1_minus_P2=v[['Endodermis_P1']]-v[['Endodermis_P2']],
               Delta_P1_minus_P3=v[['Endodermis_P1']]-v[['Endodermis_P3']],
               P1_is_highest=v[['Endodermis_P1']]>max(v[c('Endodermis_P2','Endodermis_P3')]))
  }))
  if(is.null(details) || !nrow(details)) {
    csv(data.frame(Status='Insufficient_three_stage_sample_coverage'),paste0(prefix,'_status.csv'))
    return(invisible(NULL))
  }
  csv(details,paste0(prefix,'_per_sample.csv'))
  evidence <- bind(lapply(split(details,interaction(details$Program,details$Condition,drop=TRUE)),function(a)
    data.frame(Program=a$Program[1],Condition=a$Condition[1],N_samples=nrow(a),
               Mean_P1_minus_P2=mean(a$Delta_P1_minus_P2),Mean_P1_minus_P3=mean(a$Delta_P1_minus_P3),
               Both_replicates_support_P1=nrow(a)==EXPECTED_REPLICATES && all(a$P1_is_highest),
               Basis='Within-Endodermis stage preference; descriptive, not exclusive cell-type specificity')))
  csv(evidence,paste0(prefix,'_evidence.csv'))
}
check_previous_program2 <- function(x,md) {
  if(!file.exists(OLD_MODULE_FILE)) {
    csv(data.frame(Status='Missing_previous_module_table',Path=OLD_MODULE_FILE),'previous_Program02_status.csv')
    return(invisible(NULL))
  }
  old <- read_csv(OLD_MODULE_FILE)
  if(!all(c('Program','gene_name') %in% names(old))) stop('旧模块表缺少Program或gene_name列。')
  genes <- unique(as.character(old$gene_name[old$Program=='DevProgram02' & !is.na(old$Program)]))
  genes <- intersect(genes,rownames(x))
  if(length(genes)<5) {
    csv(data.frame(Status='Insufficient_previous_DevProgram02_genes',Genes=length(genes)),
        'previous_Program02_status.csv'); return(invisible(NULL))
  }
  # 旧Program02基因集合保持固定，新分期独立评估；绝不按评分修改P1身份。
  x <- Seurat::AddModuleScore(x,features=list(genes),assay='RNA',nbin=SEURAT_NBIN,
                              ctrl=SEURAT_CTRL,seed=SEED,name='PreviousDevProgram02Score')
  sc <- as.matrix(x[[]][,'PreviousDevProgram02Score1',drop=FALSE])
  colnames(sc) <- 'Previous_DevProgram02'
  sample_scores <- aggregate_scores(sc,md[md$Analysis_eligible,,drop=FALSE])
  csv(sample_scores,'previous_Program02_sample_scores.csv')
  p1_module_evidence(sample_scores,'previous_Program02_P1')
  plot_atlas(sc,md,'Figure11_previous_Program02','Seurat module score')
  a <- sample_scores[grepl('^Endodermis_P[123]$',sample_scores$Celltype) & sample_scores$n_cells>=MIN_CELLS,]
  if(nrow(a)) {
    a$Period <- factor(sub('Endodermis_','',a$Celltype),c('P1','P2','P3'))
    g <- ggplot2::ggplot(a,ggplot2::aes(Period,Score,color=Sample,group=Sample))+
      ggplot2::geom_line()+ggplot2::geom_point(size=2)+ggplot2::facet_wrap(~Condition,nrow=1)+
      ggplot2::theme_classic(base_size=13)+ggplot2::labs(title='Previous DevProgram02: independent check across new P1/P2/P3',
                                                         subtitle='Each line is one biological replicate; stage labels were defined without module scores')
    save_plot(g,'Figure11_previous_Program02_stage_check',12,4)
  }
  csv(data.frame(Status='Scored_fixed_previous_gene_set',Genes=length(genes)), 'previous_Program02_status.csv')
}


# -----------------------------------------------------------------------------
# 12. ME / Hub
# -----------------------------------------------------------------------------
validate_harmony_design <- function(md) {
  if (!HARMONIZE_MES) return(NULL)
  if (is.null(HARMONY_COL) || !'Batch' %in% names(md))
    stop('需要真实的、跨处理批次列HARMONY_COL。')
  s <- unique(md[, c('Sample','Condition','Batch')])
  if (anyNA(s$Batch) || anyDuplicated(s$Sample) || length(unique(s$Batch)) < 2)
    stop('批次缺失、并非样本级、或仅一个批次。')
  mm <- model.matrix(~ Condition + Batch, s)
  if (qr(mm)$rank < ncol(mm)) stop('批次与处理混淆；不执行Harmony。')
  'Batch'
}

compute_mes_and_hubs <- function(x, md) {
  if (!COMPUTE_MES && !COMPUTE_HUBS) return(x)
  result <- tryCatch({
    group <- validate_harmony_design(md)
    x <- hdWGCNA::ModuleEigengenes(x, group.by.vars = group,
                                   wgcna_name = NETWORK_NAME, assay = 'RNA', npcs = 5)
    
    me <- as.matrix(hdWGCNA::GetMEs(x, harmonized = FALSE, wgcna_name = NETWORK_NAME))
    saveRDS(me, out('module_eigengenes.rds'))
    plot_atlas(me, md, 'Figure09_MEs', 'ME')
    
    if (!is.null(group)) {
      hme <- as.matrix(hdWGCNA::GetMEs(x, harmonized = TRUE, wgcna_name = NETWORK_NAME))
      saveRDS(hme, out('harmonized_module_eigengenes.rds'))
      plot_atlas(hme, md, 'Figure09_hMEs', 'hME')
    }
    
    if (COMPUTE_HUBS) {
      x <- hdWGCNA::ModuleConnectivity(
        x, group_name = 'All', group.by = 'ers_pool',
        assay = 'RNA', wgcna_name = NETWORK_NAME,
        harmonized = HARMONIZE_MES, reassign_modules = FALSE)
      modules <- hdWGCNA::GetModules(x, wgcna_name = NETWORK_NAME)
      hubs <- list()
      for (m in setdiff(levels(modules$module), 'grey')) {
        co <- paste0('kME_', m)
        if (!co %in% names(modules)) stop('缺少kME列：', co)
        d <- modules[modules$module == m & is.finite(modules[[co]]), ]
        d <- d[order(-d[[co]], d$gene_name), ]; d <- head(d, 50)
        hubs[[m]] <- data.frame(Program = m, Rank = seq_len(nrow(d)),
                                Gene = d$gene_name, kME = d[[co]])
      }
      csv(bind(hubs), 'module_Top50_kME.csv')
    }
    
    csv(data.frame(Status = 'OK', Harmony = HARMONIZE_MES), 'ME_status.csv')
    x
  }, error = function(e) {
    csv(data.frame(Status = 'Error', Message = conditionMessage(e)), 'ME_status.csv')
    note('ME/Hub步骤未完成：', conditionMessage(e))
    x
  })
  result
}


annotate_celltypes <- function(samples) {
  evidence <- list(); summaries <- list(); n <- 0L
  for (cond in CONDITIONS) for (pr in sort(unique(samples$Program))) {
    d <- samples[samples$Program==pr & samples$Condition==cond & samples$n_cells>=MIN_CELLS,]
    ns <- unique(samples$Sample[samples$Condition==cond])
    complete <- names(which(table(d$Celltype)==length(ns)))
    d <- d[d$Celltype %in% complete,]
    label <- 'Insufficient cell-type coverage'
    top <- NA_character_; chosen <- character(); ev <- data.frame()
    if (length(ns)==EXPECTED_REPLICATES && length(complete)>=3L) {
      parts <- lapply(ns,function(s) {
        q <- d[d$Sample==s,]; spread <- stats::sd(q$Score)
        q$Z <- if (is.finite(spread) && spread>1e-10) (q$Score-mean(q$Score))/spread else 0
        q$Delta_vs_other_types <- q$Score-(sum(q$Score)-q$Score)/(nrow(q)-1)
        q
      })
      d <- bind_rows(parts)
      ev <- bind_rows(lapply(split(d,d$Celltype),function(q) data.frame(
        Program=pr,Condition=cond,Celltype=q$Celltype[1],
        N_samples=nrow(q),Min_cells=min(q$n_cells),Mean_score=mean(q$Score),
        Mean_Z=mean(q$Z),Min_Z=min(q$Z),
        Mean_delta_vs_other_types=mean(q$Delta_vs_other_types),
        Both_samples_above_others=all(q$Delta_vs_other_types>0),
        Candidate=all(q$Delta_vs_other_types>0 & q$Z>=CELLTYPE_Z_CUT) && mean(q$Score)>0)))
      ev <- ev[order(-ev$Mean_score),]; ev$Rank <- seq_len(nrow(ev))
      chosen <- ev$Celltype[ev$Candidate]; top <- ev$Celltype[1]
      label <- if (!length(chosen)) paste0('No clear preference (top mean: ',top,')')
      else if (length(chosen)>MAX_LABEL_CELLTYPES) 'Multiple cell-type preference'
      else paste0(paste(chosen,collapse=' / '),' preference')
      evidence[[paste(cond,pr)]] <- ev
    }
    n <- n+1L
    summaries[[n]] <- data.frame(Program=pr,Condition=cond,Label=label,
                                 Top_celltype=top,Candidate_celltypes=paste(chosen,collapse=';'),
                                 N_candidate_celltypes=length(chosen),N_testable_celltypes=length(complete),
                                 Basis='Descriptive; sample-balanced; not a cell-type-specificity significance test')
  }
  list(summary=bind_rows(summaries),evidence=bind_rows(evidence))
}


normalize_id <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- sub('^GLYMA[._-]?', 'GLYMA.', x)
  # 只清除标准Glyma位点的转录本后缀；不去除v2/v4等基因组版本。
  sub('^(GLYMA[.][0-9]+G[0-9]+)([.](T)?[0-9]+)$','\\1',x)
}
map_ids <- function(genes, keys, crosswalk=NULL) {
  keys <- unique(as.character(keys)); keys <- keys[!is.na(keys) & nzchar(keys)]
  direct <- split(keys,toupper(trimws(keys)))
  normalized <- split(keys,normalize_id(keys))
  rows <- lapply(genes,function(g) {
    hit <- if (g %in% keys) g else direct[[toupper(trimws(g))]]
    how <- 'Exact_or_case'
    if (!length(hit)) { hit <- normalized[[normalize_id(g)]]; how <- 'Locus_format_or_transcript' }
    if (!length(hit) && !is.null(crosswalk)) {
      hit <- unique(as.character(crosswalk$DB_ID[crosswalk$Gene==g]))
      how <- 'User_crosswalk'
      if (length(hit) && any(!hit %in% keys)) return(data.frame(
        Gene=g,DB_ID=NA_character_,Mapping_status='Crosswalk_key_absent'))
    }
    if (length(hit)!=1L) {
      how <- if (length(hit)>1L) 'Ambiguous' else 'Unmapped'; hit <- NA_character_
    }
    data.frame(Gene=g,DB_ID=hit,Mapping_status=how)
  })
  bind_rows(rows)
}
get_orgdb <- function(pkg) {
  ns <- asNamespace(pkg)
  db <- get0(pkg,envir=ns,inherits=FALSE)
  if (inherits(db,'OrgDb')) return(db)
  possible <- grep('^org.*[.]db$',getNamespaceExports(ns),value=TRUE)
  objects <- lapply(possible,function(n) getExportedValue(pkg,n))
  objects <- objects[vapply(objects,inherits,logical(1),'OrgDb')]
  if (length(objects)!=1L) stop('包中无法唯一识别OrgDb对象：',pkg)
  objects[[1]]
}
bp_annotation <- function(db,ids,keytype) {
  cols <- AnnotationDbi::columns(db)
  if (all(c('GOALL','ONTOLOGYALL') %in% cols)) {
    a <- suppressMessages(AnnotationDbi::select(db,keys=ids,keytype=keytype,columns=c('GOALL','ONTOLOGYALL')))
    a <- a[a$ONTOLOGYALL %in% 'BP' & !is.na(a$GOALL),]
    return(unique(data.frame(DB_ID=a[[keytype]],GO=a$GOALL)))
  }
  if (!'GO' %in% cols) stop('OrgDb没有GO或GOALL列。')
  a <- suppressMessages(AnnotationDbi::select(db,keys=ids,keytype=keytype,columns='GO'))
  a <- unique(data.frame(DB_ID=a[[keytype]],GO=a$GO))
  a <- a[!is.na(a$GO),]
  bpkeys <- AnnotationDbi::mappedkeys(GO.db::GOBPANCESTOR)
  a <- a[a$GO %in% bpkeys,]
  if (!nrow(a)) return(a)
  # 自建OrgDb没有GOALL时，显式补BP祖先，与GOALL分析口径一致。
  gs <- split(a$DB_ID,a$GO)
  expanded <- bind_rows(lapply(names(gs),function(g) {
    ancestors <- GO.db::GOBPANCESTOR[[g]]
    terms <- unique(c(g,ancestors[!is.na(ancestors) & grepl('^GO:',ancestors)]))
    expand.grid(DB_ID=unique(gs[[g]]),GO=terms,stringsAsFactors=FALSE)
  }))
  unique(expanded)
}
choose_db <- function(background,write_table) {
  need(c('AnnotationDbi','GO.db'))
  audits <- list(); candidates <- list(); detail <- list()
  crosswalk <- NULL
  if (!is.null(GENE_ID_MAP_FILE)) {
    crosswalk <- read_csv(GENE_ID_MAP_FILE)
    if (!all(c('Gene','DB_ID','Package','Keytype') %in% names(crosswalk)) || anyNA(crosswalk))
      stop('映射表必须含非缺失Gene/DB_ID/Package/Keytype。')
  }
  for (pkg in ORGDB_PACKAGES) {
    if (!requireNamespace(pkg,quietly=TRUE)) {
      audits[[length(audits)+1L]] <- data.frame(Package=pkg,Version=NA_character_,Keytype=NA_character_,
                                                Total=length(background),Mapped=0L,Mapping_rate=0,BP_mapped=0L,Status='Not_installed')
      next
    }
    db <- tryCatch(get_orgdb(pkg),error=function(e)e)
    if (inherits(db,'error')) {
      audits[[length(audits)+1L]] <- data.frame(Package=pkg,Version=as.character(utils::packageVersion(pkg)),
                                                Keytype=NA_character_,Total=length(background),Mapped=0L,Mapping_rate=0,BP_mapped=0L,
                                                Status=conditionMessage(db)); next
    }
    ks <- AnnotationDbi::keytypes(db)
    kt <- if (!is.null(KEYTYPE_OVERRIDE)) intersect(KEYTYPE_OVERRIDE,ks) else
      intersect(c('GID','SYMBOL','ENSEMBL','ENTREZID','GENEID','ORF','TAIR','REFSEQ'),ks)
    if (!length(kt)) kt <- if (is.null(KEYTYPE_OVERRIDE)) head(ks[!grepl('GO|PATH|EVIDENCE|ONTOLOGY',ks)],1) else character()
    if (!length(kt)) {
      audits[[length(audits)+1L]] <- data.frame(Package=pkg,Version=as.character(utils::packageVersion(pkg)),
                                                Keytype=NA_character_,Total=length(background),Mapped=0L,Mapping_rate=0,BP_mapped=0L,
                                                Status='No_usable_keytype'); next
    }
    for (keytype in kt) {
      mp <- tryCatch(map_ids(background,AnnotationDbi::keys(db,keytype=keytype),
                             if (is.null(crosswalk)) NULL else crosswalk[crosswalk$Package==pkg & crosswalk$Keytype==keytype,]),
                     error=function(e)e)
      if (inherits(mp,'error')) {
        audits[[length(audits)+1L]] <- data.frame(Package=pkg,Version=as.character(utils::packageVersion(pkg)),
                                                  Keytype=keytype,Total=length(background),Mapped=0L,Mapping_rate=0,BP_mapped=0L,
                                                  Status=conditionMessage(mp)); next
      }
      ids <- unique(na.omit(mp$DB_ID))
      bp <- if (length(ids)) tryCatch(bp_annotation(db,ids,keytype),error=function(e)e) else data.frame()
      n_bp <- if (is.data.frame(bp) && nrow(bp)) sum(mp$DB_ID %in% bp$DB_ID) else 0L
      mapped <- sum(!is.na(mp$DB_ID)); rate <- mapped/length(background)
      status <- if (inherits(bp,'error')) conditionMessage(bp) else if (!n_bp) 'No_BP_annotation'
      else if (rate<MIN_DB_MAPPING_RATE) 'Low_mapping_rate' else 'Usable'
      audits[[length(audits)+1L]] <- data.frame(Package=pkg,Version=as.character(utils::packageVersion(pkg)),
                                                Keytype=keytype,Total=length(background),Mapped=mapped,Mapping_rate=rate,BP_mapped=n_bp,Status=status)
      mp$Package <- pkg; mp$Keytype <- keytype
      detail[[paste(pkg,keytype)]] <- mp
      if (status=='Usable') candidates[[paste(pkg,keytype)]] <- list(
        db=db,pkg=pkg,keytype=keytype,mapping=mp,bp=bp,rate=rate,
        exact=sum(mp$Mapping_status=='Exact_or_case'))
    }
  }
  audit <- bind_rows(audits)
  write_table(audit,'OrgDb_mapping_diagnostics.csv')
  write_table(bind_rows(detail),'OrgDb_all_mapping_attempts.csv')
  if (!length(candidates)) stop('两个OrgDb均无可用匹配；查看OrgDb_mapping_diagnostics.csv。不会把映射失败写成无生物学富集。')
  ord <- order(vapply(candidates,function(z)match(z$pkg,ORGDB_PACKAGES),integer(1)),
               -vapply(candidates,`[[`,numeric(1),'rate'),
               -vapply(candidates,`[[`,integer(1),'exact'))
  chosen <- candidates[[ord[1]]]
  chosen$audit <- audit
  chosen
}
go_analysis <- function(mods,background,write_table) {
  z <- choose_db(background,write_table)
  mp <- z$mapping; db <- z$db; kt <- z$keytype
  write_table(data.frame(Package=z$pkg,Version=as.character(utils::packageVersion(z$pkg)),
                         Keytype=kt,Mapping_rate=z$rate,Selection='Package priority, then ID match rate; never GO P values'),
              'GO_database_selected.csv')
  write_table(AnnotationDbi::metadata(db),'GO_database_metadata.csv')
  mp$Description <- NA_character_
  desc_cols <- intersect(c('GENENAME','DESCRIPTION','GENE_NAME','PRODUCT','DEFINITION'),AnnotationDbi::columns(db))
  ids <- unique(na.omit(mp$DB_ID))
  for (co in desc_cols) {
    a <- suppressMessages(AnnotationDbi::select(db,keys=ids,keytype=kt,columns=co))
    text <- vapply(split(as.character(a[[co]]),a[[kt]]),function(v)
      paste(unique(v[!is.na(v) & nzchar(v)]),collapse=' | '),character(1))
    i <- which(is.na(mp$Description)); val <- unname(text[mp$DB_ID[i]])
    val[!is.na(val) & !nzchar(val)] <- NA_character_; mp$Description[i] <- val
  }
  ann <- unique(z$bp)
  universe <- sort(unique(intersect(ids,ann$DB_ID)))
  mp$Has_BP <- mp$DB_ID %in% universe
  mp$Many_input_genes_per_DB_ID <- !is.na(mp$DB_ID) &
    (duplicated(mp$DB_ID) | duplicated(mp$DB_ID,fromLast=TRUE))
  write_table(mp,'gene_ID_mapping.csv')
  annotated <- merge(mods,mp,by.x='gene_name',by.y='Gene',all.x=TRUE,sort=FALSE)
  write_table(annotated,'modules_gene_annotation.csv')
  sets <- lapply(split(ann$DB_ID,ann$GO),function(a) intersect(unique(a),universe))
  sets <- sets[lengths(sets)>=MIN_GS_SIZE & lengths(sets)<=MAX_GS_SIZE]
  term_ids <- intersect(names(sets),AnnotationDbi::keys(GO.db::GO.db,keytype='GOID'))
  terms <- if (length(term_ids)) suppressMessages(AnnotationDbi::select(
    GO.db::GO.db,keys=term_ids,keytype='GOID',columns='TERM')) else data.frame(GOID=character(),TERM=character())
  all <- list(); status <- list()
  for (p in sort(unique(na.omit(mods$Program)))) {
    genes <- mods$gene_name[which(mods$Program==p)]
    m <- mp[match(genes,mp$Gene),]
    q <- intersect(unique(na.omit(m$DB_ID)),universe)
    state <- if (length(q)<MIN_QUERY_GENES) 'Too_few_BP_genes' else if (!length(sets)) 'No_testable_terms' else 'Tested'
    status[[p]] <- data.frame(Program=p,Network_genes=length(genes),
                              Mapped_input_genes=sum(!is.na(m$DB_ID)),Mapped_BP_unique_IDs=length(q),
                              Background_BP_unique_IDs=length(universe),Status=state)
    if (state!='Tested') next
    overlap <- lapply(sets,intersect,y=q); n <- lengths(overlap); K <- lengths(sets)
    pv <- stats::phyper(n-1,K,length(universe)-K,length(q),lower.tail=FALSE)
    label <- terms$TERM[match(names(sets),terms$GOID)]
    label[is.na(label)] <- names(sets)[is.na(label)]
    all[[p]] <- data.frame(Program=p,GO=names(sets),Description=label,
                           Count=n,Query_BP=length(q),Background_BP=length(universe),Term_background=K,
                           GeneRatio=n/length(q),Fold_enrichment=(n/length(q))/(K/length(universe)),
                           P_value=pv,FDR_within=p.adjust(pv,'BH'),
                           DB_IDs=vapply(overlap,paste,collapse=';',character(1)))
  }
  tab <- if (length(all)) bind_rows(all) else data.frame(Program=character(),GO=character(),
                                                         Description=character(),Count=integer(),Query_BP=integer(),Background_BP=integer(),Term_background=integer(),
                                                         GeneRatio=numeric(),Fold_enrichment=numeric(),P_value=numeric(),FDR_within=numeric(),DB_IDs=character())
  tab$FDR_global <- p.adjust(tab$P_value,'BH')
  sig <- tab[tab$Count>0 & tab$FDR_global<=FDR_CUT,]
  st <- bind_rows(status)
  st$Significant_terms <- vapply(st$Program,function(p) sum(sig$Program==p),integer(1))
  st$Status[st$Status=='Tested' & st$Significant_terms==0] <- 'Tested_no_significant_terms'
  write_table(tab,'GO_BP_all_tests.csv'); write_table(sig,'GO_BP_significant.csv')
  write_table(st,'GO_status.csv')
  list(all=tab,significant=sig,status=st,mapping=mp)
}


# ====================== 内皮层参考发育轨迹 ======================
balance_cells <- function(md,types,max_n) {
  ix <- which(md$Condition=='CK' & md$Celltype %in% types)
  g <- split(ix,md$Sample[ix])
  if (!setequal(names(g),c('CK1','CK2'))) stop('指定参考细胞在CK1/CK2中不完整。')
  n <- min(max_n,lengths(g)); if(n<MIN_CELLS) stop('CK参考细胞过少。')
  set.seed(SEED); unlist(lapply(g,sample,size=n),use.names=FALSE)
}
cluster_means <- function(pca,clusters) {
  z <- t(vapply(split(seq_len(nrow(pca)),clusters),function(i)
    colMeans(pca[i,,drop=FALSE]),numeric(ncol(pca))))
  colnames(z) <- colnames(pca); z
}
fit_reference_pca <- function(x,npcs,svd_engine=NULL) {
  # x必须是CK参考细胞×基因；明确记录变量顺序，不依赖SVD输出dimnames。
  x <- as.matrix(x)
  if(!is.numeric(x) || any(!is.finite(x))) stop('CK PCA训练矩阵含非有限值或不是数值。')
  genes <- colnames(x)
  if(is.null(genes) || anyNA(genes) || any(!nzchar(genes)) || anyDuplicated(genes))
    stop('CK PCA训练矩阵基因列名缺失或重复。')
  if(nrow(x)<5L) stop('CK PCA参考细胞不足。')
  mu <- colMeans(x); sigma <- apply(x,2L,stats::sd)
  keep <- is.finite(sigma) & sigma>1e-8
  audit <- data.frame(Gene=genes,Center=as.numeric(mu),Scale=as.numeric(sigma),
                      Retained=keep,Reason=ifelse(keep,'Retained','Zero_or_near_zero_variance'))
  x <- x[,keep,drop=FALSE]; genes <- colnames(x)
  mu <- mu[keep]; sigma <- sigma[keep]
  np <- min(as.integer(npcs),nrow(x)-2L,ncol(x)-2L)
  if(!is.finite(np) || np<3L) stop('过滤零方差基因后，CK参考PCA可用维度不足。')
  z <- sweep(sweep(x,2L,mu,'-'),2L,sigma,'/')
  if(any(!is.finite(z))) stop('CK训练矩阵标准化产生非有限值。')
  # 明确在已标准化矩阵上做截断SVD；不使用prcomp_irlba的center/scale字段。
  if(is.null(svd_engine)) svd_engine <- function(z,n)
    irlba::irlba(z,nv=n,nu=0L,maxit=1000L,tol=1e-6)
  set.seed(SEED)
  s <- svd_engine(z,np)
  rotation <- as.matrix(s$v)
  if(nrow(rotation)!=length(genes) || ncol(rotation)!=np)
    stop('SVD载荷维度异常：预期 ',length(genes),'×',np,'；实际 ',paste(dim(rotation),collapse='×'))
  if(any(!is.finite(rotation)) || length(s$d)!=np || any(!is.finite(s$d)))
    stop('SVD返回无效载荷或奇异值。')
  # SVD的V行顺序与输入变量顺序相同；在尺寸检查后显式附加基因名。
  pcs <- paste0('CKPC_',seq_len(np))
  dimnames(rotation) <- list(genes,pcs)
  names(mu) <- names(sigma) <- genes
  projected <- z %*% rotation
  rownames(projected) <- rownames(x); colnames(projected) <- pcs
  fit <- structure(list(rotation=rotation,sdev=as.numeric(s$d)/sqrt(nrow(x)-1),
                        center=mu,scale=sigma,x=projected,training_genes=genes,
                        projection_version='explicit_standardization_v4.1'),class='prcomp')
  list(fit=fit,audit=audit)
}
project_reference_pca <- function(expr,allcells,fit,chunk_size=1000L) {
  genes <- fit$training_genes; rotation <- fit$rotation
  if(!length(genes) || anyDuplicated(genes) || !identical(genes,rownames(rotation)))
    stop('PCA训练基因与载荷行名不一致。')
  if(!is.numeric(fit$center) || !is.numeric(fit$scale) ||
     length(fit$center)!=length(genes) || length(fit$scale)!=length(genes) ||
     !identical(names(fit$center),genes) || !identical(names(fit$scale),genes) ||
     any(!is.finite(fit$center)) || any(!is.finite(fit$scale)) || any(fit$scale<=0))
    stop('PCA中心/标准差必须是与训练基因同序的有限数值向量。')
  if(is.null(rownames(expr)) || anyDuplicated(rownames(expr)) ||
     !all(genes %in% rownames(expr))) stop('投射RNA矩阵缺少训练基因或基因名重复。')
  if(is.null(colnames(expr)) || anyDuplicated(colnames(expr)) ||
     anyNA(allcells) || anyDuplicated(allcells) || !all(allcells %in% colnames(expr)))
    stop('投射RNA矩阵细胞名称缺失、重复或不匹配。')
  scores <- matrix(NA_real_,length(allcells),ncol(rotation),
                   dimnames=list(allcells,colnames(rotation)))
  chunks <- split(seq_along(allcells),ceiling(seq_along(allcells)/chunk_size))
  for (ix in chunks) {
    y <- t(as.matrix(expr[genes,allcells[ix],drop=FALSE]))
    if(!identical(colnames(y),genes) || ncol(y)!=nrow(rotation))
      stop('PCA投射维度/基因顺序错误：query=',paste(dim(y),collapse='×'),
           '；loadings=',paste(dim(rotation),collapse='×'))
    if(any(!is.finite(y))) stop('待投射表达矩阵含非有限值。')
    y <- sweep(sweep(y,2L,fit$center,'-'),2L,fit$scale,'/')
    scores[ix,] <- y %*% rotation
  }
  if(any(!is.finite(scores))) stop('PCA投射产生非有限结果。')
  scores
}
reference_pca <- function(expr,genes,refcells,allcells) {
  if(is.null(rownames(expr)) || is.null(colnames(expr)) ||
     anyDuplicated(rownames(expr)) || anyDuplicated(colnames(expr)))
    stop('RNA表达矩阵必须有唯一基因名和细胞名。')
  if(anyNA(refcells) || anyDuplicated(refcells) || !all(refcells %in% allcells))
    stop('CK参考细胞重复、缺失或不属于投射细胞集合。')
  genes <- unique(as.character(genes)); genes <- genes[!is.na(genes) & nzchar(genes)]
  absent <- setdiff(genes,rownames(expr))
  genes <- genes[genes %in% rownames(expr)]
  if(length(genes)<5L || !all(refcells %in% colnames(expr))) stop('CK PCA基因或细胞不匹配。')
  x <- t(as.matrix(expr[genes,refcells,drop=FALSE]))
  trained <- fit_reference_pca(x,TRAJECTORY_N_PCS)
  fit <- trained$fit
  scores <- project_reference_pca(expr,allcells,fit)
  # 同一CK矩阵经分块投射，必须与训练时计算的分数一致。
  error <- max(abs(scores[refcells,,drop=FALSE]-fit$x))
  if(!is.finite(error) || error>1e-8*max(1,max(abs(fit$x))))
    stop('CK训练得分与重新投射得分不一致。')
  csv(trained$audit,'trajectory_PCA_gene_scaling.csv')
  csv(data.frame(Gene=absent),'trajectory_PCA_absent_HVG.csv')
  csv(data.frame(Reference_cells=nrow(x),Requested_present_genes=ncol(x),
                 Retained_genes=nrow(fit$rotation),Returned_PCs=ncol(fit$rotation),
                 Projection_cells=nrow(scores),Reference_projection_max_error=error,
                 Method='Explicit CK centering/scaling + truncated SVD; named feature alignment'),
      'trajectory_PCA_projection_QC.csv')
  note('CK PCA检查通过：',nrow(fit$rotation),' 个基因，',ncol(fit$rotation),
       ' 个PC；投射 ',nrow(scores),' 个细胞。')
  list(fit=fit,scores=scores)
}
# 固定三个CK表达聚类；P1/P2/P3不是拟时间三分位，也不使用模块评分决定分组。
nearest_centroid <- function(x,centers) {
  x <- as.matrix(x); centers <- as.matrix(centers)
  if(ncol(x)!=ncol(centers) || any(!is.finite(x)) || any(!is.finite(centers)))
    stop('最近中心分类：输入维度或数值不合法。')
  d2 <- vapply(seq_len(nrow(centers)),function(j)
    rowSums(sweep(x,2,centers[j,],'-')^2),numeric(nrow(x)))
  if(is.null(dim(d2))) d2 <- matrix(d2,nrow=nrow(x))
  best <- max.col(-d2,ties.method='first')
  ord <- t(apply(d2,1,sort))
  list(cluster=rownames(centers)[best],distance=sqrt(d2[cbind(seq_len(nrow(x)),best)]),
       margin=(ord[,2]-ord[,1])/pmax(ord[,2],1e-12))
}
project_fixed_path <- function(x,path) {
  # 固定参考曲线上的弧长投射：起点不随query样本/分块而重置。
  x <- as.matrix(x); path <- as.matrix(path)
  if(ncol(x)!=ncol(path) || any(!is.finite(x)) || any(!is.finite(path))) stop('轨迹投射输入异常。')
  seg <- path[-1,,drop=FALSE]-path[-nrow(path),,drop=FALSE]
  lens <- sqrt(rowSums(seg^2)); good <- which(lens>1e-10)
  if(!length(good)) stop('参考曲线没有长度。')
  starts <- c(0,cumsum(lens))
  best <- rep(Inf,nrow(x)); pt <- rep(NA_real_,nrow(x))
  for(j in good) {
    dx <- sweep(x,2,path[j,],'-')
    u <- pmin(1,pmax(0,as.vector(dx %*% seg[j,])/lens[j]^2))
    residual <- dx-tcrossprod(u,seg[j,])
    dd <- rowSums(residual^2); replace <- dd<best
    best[replace] <- dd[replace]; pt[replace] <- starts[j]+u[replace]*lens[j]
  }
  list(pseudotime=pt,distance=sqrt(best),length=sum(lens))
}
fit_three_clusters <- function(rpc,anchor_centers,root_cells=NULL) {
  set.seed(SEED)
  fit <- stats::kmeans(rpc,centers=3L,nstart=50L,iter.max=200L,algorithm='Hartigan-Wong')
  centers <- fit$centers
  if(length(unique(fit$cluster))!=3L) stop('固定三聚类未形成三个非空组。')
  distance <- vapply(seq_len(3),function(j)
    mean(sqrt(rowSums(sweep(anchor_centers,2,centers[j,],'-')^2))),numeric(1))
  root <- ROOT_CLUSTER; basis <- 'CK_Meristem_QC_expression_proximity'
  if(is.null(root) && length(root_cells)) {
    hit <- intersect(root_cells,rownames(rpc))
    if(length(hit)<MIN_ROOT_CELLS) stop('指定早期条码与CK参考匹配不足。')
    tab <- sort(table(fit$cluster[match(hit,rownames(rpc))]),decreasing=TRUE)
    if(as.numeric(tab[1])/sum(tab)<=0.5) stop('早期条码不能唯一支持一个聚类。')
    root <- names(tab)[1]; basis <- 'Validated_early_barcodes'
  }
  if(is.null(root)) root <- rownames(centers)[which.min(distance)]
  else if(is.null(root_cells)) basis <- 'Explicit_root_cluster'
  root <- as.character(root)
  if(length(root)!=1L || !root %in% rownames(centers)) stop('ROOT_CLUSTER必须是新三聚类的1、2或3，不能沿用旧cluster5。')
  remain <- setdiff(rownames(centers),root)
  d <- rowSums(sweep(centers[remain,,drop=FALSE],2,centers[root,],'-')^2)
  ordered <- c(root,remain[order(d)])
  map <- setNames(c('P1','P2','P3'),ordered)
  audit <- data.frame(Cluster=rownames(centers),Period=unname(map[rownames(centers)]),
                      Distance_to_early_anchors=distance,CK_reference_cells=as.integer(table(factor(fit$cluster,levels=1:3))),
                      Root_basis=basis,Ordering='P1 nearest early anchors; shortest centroid path to other two clusters')
  list(fit=fit,centers=centers,ordered=ordered,map=map,audit=audit,root=root)
}
endo_trajectory <- function(dat) {
  m <- dat$md
  refix <- balance_cells(m,ENDO_CELLTYPE,REFERENCE_MAX_CELLS_PER_CK)
  endo <- which(m$Celltype==ENDO_CELLTYPE)
  anchors <- which(m$Condition=='CK' & m$Celltype %in% EARLY_ANCHOR_TYPES)
  if(!length(endo) || !all(c('CK1','CK2') %in% m$Sample[anchors])) stop('缺少内皮层或CK早期参照细胞。')
  cells <- m$cell[unique(c(endo,anchors))]
  obj <- SeuratObject::CreateSeuratObject(counts=SeuratObject::CreateAssayObject(counts=dat$counts[,cells,drop=FALSE]),
                                          meta.data=m[cells,,drop=FALSE])
  obj <- Seurat::NormalizeData(obj,verbose=FALSE)
  ref <- subset(obj,cells=m$cell[refix])
  ref <- Seurat::FindVariableFeatures(ref,nfeatures=TRAJECTORY_N_HVG,verbose=FALSE)
  expr <- get_logdata(obj)
  p <- reference_pca(expr,SeuratObject::VariableFeatures(ref),m$cell[refix],cells)
  rpc <- p$scores[m$cell[refix],,drop=FALSE]
  apc <- p$scores[m$cell[anchors],,drop=FALSE]
  ac <- cluster_means(apc,paste(m$Sample[anchors],m$Celltype[anchors],sep='_'))
  roots <- NULL
  if(!is.null(ROOT_CELLS_FILE)) {
    root_table <- read_csv(ROOT_CELLS_FILE)
    if(!'cell' %in% names(root_table)) stop('ROOT_CELLS_FILE必须有cell列。')
    roots <- as.character(root_table$cell)
  }
  cl <- fit_three_clusters(rpc,ac,roots)
  csv(cl$audit,'trajectory_root_and_three_cluster_order.csv')
  csv(as.data.frame(table(Sample=m$Sample[refix],Period=cl$map[as.character(cl$fit$cluster)])),
      'CK_three_cluster_sample_support.csv')
  note('CK固定3个表达聚类：',paste(paste(names(cl$map),cl$map,sep='='),collapse=', '),
       '。阶段不是按module评分或拟时间等分决定。')
  epc <- p$scores[m$cell[endo],,drop=FALSE]
  assigned <- nearest_centroid(epc,cl$centers)
  z <- m[m$cell[endo],,drop=FALSE]
  z$Endo_cluster <- assigned$cluster
  z$Period <- unname(cl$map[assigned$cluster])
  z$Cluster_margin <- assigned$margin
  kk <- min(10L,nrow(rpc)-1L)
  rd <- rowMeans(FNN::get.knn(rpc,k=kk)$nn.dist)
  cutoff <- as.numeric(stats::quantile(rd,PROJECTION_DISTANCE_QUANTILE))*PROJECTION_DISTANCE_MULTIPLIER
  z$Reference_distance <- rowMeans(FNN::get.knnx(rpc,epc,k=kk)$nn.dist)
  z$Projection_pass <- z$Reference_distance<=cutoff
  z$Reference_fit_cell <- z$cell %in% rownames(rpc)
  z$Stage_status <- ifelse(z$Projection_pass,'Assigned','Assigned_but_distant_from_CK')
  z$Endo_branch <- 'CK_single_axis'
  # 连续参考轨迹独立于模块评分。拟合不可靠时明确记录并用三中心折线作近似坐标。
  init <- cl$centers[cl$ordered,,drop=FALSE]
  curve <- tryCatch(princurve::principal_curve(rpc,start=init,stretch=1,
                                               approx_points=150,maxit=30),error=function(e)e)
  method <- 'CK_principal_curve'; reason <- 'Converged_and_cluster_medians_ordered'
  path <- init
  if(inherits(curve,'error')) reason <- conditionMessage(curve)
  else if(!isTRUE(curve$converged)) reason <- 'Principal_curve_not_converged'
  else {
    candidate <- unique(curve$s[curve$ord,,drop=FALSE])
    rr <- project_fixed_path(rpc,candidate)$pseudotime
    med <- tapply(rr,cl$map[as.character(cl$fit$cluster)],median)[c('P1','P2','P3')]
    if(med[1]>med[3]) candidate <- candidate[nrow(candidate):1,,drop=FALSE]
    rr <- project_fixed_path(rpc,candidate)$pseudotime
    med <- tapply(rr,cl$map[as.character(cl$fit$cluster)],median)[c('P1','P2','P3')]
    if(all(diff(med)>1e-6)) path <- candidate else reason <- 'Curve_inconsistent_with_three_cluster_order'
  }
  if(reason!='Converged_and_cluster_medians_ordered') {
    method <- 'CK_three_centroid_path_approximation'
    note('连续轨迹使用三中心折线近似；原因：',reason,'。三个表达聚类保持不变，不能将近似坐标声称为已验证的真实发育时间。')
  }
  projected <- project_fixed_path(epc,path)
  z$Pseudotime <- projected$pseudotime
  z$Pseudotime01 <- z$Pseudotime/projected$length
  z$Curve_distance <- projected$distance
  csv(data.frame(Method=method,Reason=reason,Root_cluster=cl$root,
                 Distance_cutoff=cutoff,Note='Three expression clusters; direction provisional; no velocity matrix supplied'),
      'trajectory_model_status.csv')
  csv(z,'Endodermis_three_clusters.csv')
  csv(cl$audit,'trajectory_stage_definition.csv')
  # cluster marker为描述性均值差，不以细胞数充当生物学重复。
  markers <- bind(lapply(c('P1','P2','P3'),function(stage) {
    ids <- m$cell[refix]; stage_ref <- cl$map[as.character(cl$fit$cluster)]
    delta <- Matrix::rowMeans(expr[,ids[stage_ref==stage],drop=FALSE])-
      Matrix::rowMeans(expr[,ids[stage_ref!=stage],drop=FALSE])
    detect <- Matrix::rowMeans(expr[,ids[stage_ref==stage],drop=FALSE]>0)
    a <- data.frame(Period=stage,Gene=rownames(expr),Delta_mean_logRNA=delta,Detected_fraction=detect)
    a <- a[a$Detected_fraction>=.1,]; head(a[order(-a$Delta_mean_logRNA),],50)
  }))
  csv(markers,'CK_three_cluster_markers_descriptive.csv')
  cols <- c(P1='#3B528B',P2='#21918C',P3='#E2C900')
  plot_stages <- function(q,facet,title) {
    q$Period <- factor(q$Period,c('P1','P2','P3'))
    g <- ggplot2::ggplot(q,ggplot2::aes(UMAP1,UMAP2,color=Period))+
      ggplot2::geom_point(size=.55,alpha=.85)+ggplot2::facet_wrap(stats::as.formula(paste('~',facet)),nrow=1)+
      ggplot2::scale_color_manual(values=cols,drop=FALSE)+ggplot2::coord_equal()+ggplot2::theme_classic(base_size=13)+
      ggplot2::labs(title=title,subtitle='Exactly three CK expression clusters; treatment cells assigned to fixed CK centroids')
    g
  }
  save_plot(plot_stages(z[z$Condition=='CK',],'Sample','CK Endodermis: P1 / P2 / P3'),
            'Figure00_CK_reclustering_P1P2P3',10,4.8)
  z$Condition <- factor(z$Condition,CONDITIONS)
  save_plot(plot_stages(z,'Condition','Endodermis: three clusters on original atlas'),
            'Figure01_Endo_three_clusters_atlas',14,4.5)
  eo <- subset(obj,cells=z$cell)
  eo[['ckpca']] <- SeuratObject::CreateDimReducObject(embeddings=epc,key='CKPC_',assay='RNA')
  eo <- Seurat::RunUMAP(eo,reduction='ckpca',dims=seq_len(ncol(epc)),seed.use=SEED,verbose=FALSE)
  xy <- SeuratObject::Embeddings(eo,'umap'); q <- z
  q$UMAP1 <- xy[q$cell,1]; q$UMAP2 <- xy[q$cell,2]
  save_plot(plot_stages(q,'Condition','Endodermis: three clusters on display UMAP'),
            'Figure01_Endo_three_clusters_newUMAP',14,4.5)
  for(k in c('Period','Pseudotime','Pseudotime01','Stage_status','Cluster_margin')) eo[[k]] <- z[[k]]
  saveRDS(eo,out('Endodermis_trajectory_object.rds'))
  saveRDS(list(pca=p$fit,centers=cl$centers,stage_map=cl$map,path=path,method=method),out('CK_three_cluster_reference.rds'))
  long <- z; long$Branch <- 'CK_axis'; long$Lineage <- 'CK_axis'; long$Lineage_weight <- 1
  list(metadata=z,long=long,cutoffs=data.frame(Branch='CK_axis',Status='OK'),
       root=cl$root,lineage='CK_axis',stage_definition=cl$audit,method=method)
}
apply_endo_stages <- function(dat,trajectory) {
  m <- dat$md; z <- trajectory$metadata; m$Original_celltype <- m$Celltype
  ix <- which(m$Celltype==ENDO_CELLTYPE); j <- match(m$cell[ix],z$cell)
  if(anyNA(j) || anyDuplicated(z$cell)) stop('三聚类结果细胞名异常。')
  for(k in c('Period','Pseudotime','Pseudotime01','Endo_branch','Stage_status','Cluster_margin')) {
    m[[k]] <- z[[k]][rep(NA_integer_,nrow(m))]; m[[k]][ix] <- z[[k]][j]
  }
  m$Celltype[ix] <- paste0('Endodermis_',m$Period[ix])
  m$Network_eligible <- TRUE; m$Network_eligible[ix] <- z$Projection_pass[j]
  m$Analysis_eligible <- m$Network_eligible
  if(!setequal(unique(m$Period[ix]),c('P1','P2','P3'))) stop('内皮层未得到完整三个聚类。')
  csv(m,'all_cells_with_Endodermis_stages.csv')
  csv(as.data.frame(table(Sample=m$Sample,Celltype=m$Celltype,Eligible=m$Analysis_eligible)),
      'sample_stage_network_coverage.csv')
  dat$md <- m; dat$input_hvg <- character(); dat['graph'] <- list(NULL)
  dat$graph_name <- 'RNA_PCA_after_three_Endodermis_clusters'; dat
}

aggregate_lineage_scores <- function(scores,long) {
  if(anyDuplicated(paste(long$cell,long$Branch,sep='::'))) stop('统计输入cell×branch重复。')
  z <- long[long$Projection_pass & !is.na(long$Period),,drop=FALSE]
  empty <- data.frame(Branch=character(),Lineage=character(),Sample=character(),
                      Condition=character(),Period=character(),n_cells=integer(),Program=character(),Score=numeric(),
                      mean_logUMI=numeric(),mean_nFeature=numeric())
  if(!nrow(z)) return(empty)
  if(!all(z$cell %in% rownames(scores))) stop('分支细胞没有对应模块评分。')
  groups <- split(seq_len(nrow(z)),interaction(z$Branch,z$Sample,z$Period,drop=TRUE))
  do.call(rbind,lapply(groups,function(ix) {
    a <- z[ix,,drop=FALSE]
    data.frame(Branch=a$Branch[1],Lineage=a$Lineage[1],Sample=a$Sample[1],
               Condition=a$Condition[1],Period=a$Period[1],n_cells=nrow(a),
               Program=colnames(scores),Score=colMeans(scores[a$cell,,drop=FALSE]),
               mean_logUMI=mean(log1p(a$nUMI)),mean_nFeature=mean(a$nFeature),stringsAsFactors=FALSE)
  }))
}


stage_interaction_design <- function(sm) {
  sm$Sample <- factor(sm$Sample)
  X <- stats::model.matrix(~0+Sample,data=sm)
  for(p in c('P2','P3')) X <- cbind(X,as.integer(sm$Period==p))
  colnames(X)[(ncol(X)-1):ncol(X)] <- c('CK_P2_vs_P1','CK_P3_vs_P1')
  for(tr in TREATMENTS) for(p in c('P2','P3')) {
    X <- cbind(X,as.integer(sm$Condition==tr & sm$Period==p))
    colnames(X)[ncol(X)] <- paste0(tr,'_',p,'_interaction')
  }
  X
}
stage_dynamics <- function(scores,md,trajectory) {
  programs <- colnames(scores)
  allscores <- aggregate_lineage_scores(scores,trajectory$long)
  csv(allscores,'Endodermis_sample_stage_module_scores.csv')
  d <- allscores[allscores$n_cells>=MIN_CELLS,,drop=FALSE]
  branches <- trajectory$cutoffs$Branch
  expected <- expand.grid(Sample=unique(md$Sample),Period=paste0('P',1:3),stringsAsFactors=FALSE)
  expected$Key <- paste(expected$Sample,expected$Period,sep='::')
  tests <- bases <- statuses <- missing <- list()
  for(b in branches) {
    a <- d[d$Branch==b,,drop=FALSE]
    sm <- unique(a[,c('Sample','Condition','Period'),drop=FALSE])
    sm$Key <- paste(sm$Sample,sm$Period,sep='::')
    absent <- expected[!expected$Key %in% sm$Key,,drop=FALSE]; absent$Branch <- rep(b,nrow(absent))
    missing[[b]] <- absent
    if(nrow(absent)) {
      statuses[[b]] <- data.frame(Branch=b,Status='Skipped',Observations=nrow(sm),
                                  Independent_samples=length(unique(sm$Sample)),Residual_df=NA_integer_,
                                  Reason='Missing_sample_stage_groups; no zero filling')
      next
    }
    X <- stage_interaction_design(sm)
    if(qr(X)$rank<ncol(X) || nrow(X)-ncol(X)<2L) {
      statuses[[b]] <- data.frame(Branch=b,Status='Skipped',Observations=nrow(sm),
                                  Independent_samples=length(unique(sm$Sample)),Residual_df=nrow(X)-qr(X)$rank,
                                  Reason='Non_estimable_design'); next
    }
    Y <- matrix(NA_real_,length(programs),nrow(sm),dimnames=list(programs,sm$Key))
    Y[cbind(match(a$Program,programs),match(paste(a$Sample,a$Period,sep='::'),sm$Key))] <- a$Score
    good <- apply(Y,1,function(y) all(is.finite(y)) && stats::var(y)>1e-12)
    if(!any(good)) {
      statuses[[b]] <- data.frame(Branch=b,Status='Skipped',Observations=nrow(sm),
                                  Independent_samples=length(unique(sm$Sample)),Residual_df=nrow(X)-ncol(X),
                                  Reason='No_finite_variable_modules'); next
    }
    f <- limma::eBayes(limma::lmFit(Y[good,,drop=FALSE],X))
    base <- limma::topTable(f,coef=c('CK_P2_vs_P1','CK_P3_vs_P1'),number=Inf,sort.by='none')
    bases[[b]] <- data.frame(Branch=b,Program=rownames(base),
                             Delta_P2_vs_P1=base$CK_P2_vs_P1,Delta_P3_vs_P1=base$CK_P3_vs_P1,
                             F_statistic=base$F,P_value=base$P.Value)
    tests[[b]] <- bind(lapply(TREATMENTS,function(tr) {
      co <- paste0(tr,c('_P2_interaction','_P3_interaction'))
      tt <- limma::topTable(f,coef=co,number=Inf,sort.by='none')
      data.frame(Branch=b,Program=rownames(tt),Treatment=tr,Delta_P2_vs_P1=tt[[co[1]]],
                 Delta_P3_vs_P1=tt[[co[2]]],F_statistic=tt$F,P_value=tt$P.Value)
    }))
    statuses[[b]] <- data.frame(Branch=b,Status='OK',Observations=nrow(X),
                                Independent_samples=length(unique(sm$Sample)),Residual_df=nrow(X)-ncol(X),
                                Reason='Sample intercepts + CK stage + treatment x stage; shared cells may occur in other branch models')
  }
  tests <- bind(tests); bases <- bind(bases)
  # 联合校正所有已检验模块×分支×处理；不同分支模型可能共享早期细胞，不能当成独立验证。
  if(nrow(tests)) tests$FDR_global <- p.adjust(tests$P_value,'BH')
  else tests <- data.frame(Branch=character(),Program=character(),Treatment=character(),
                           Delta_P2_vs_P1=numeric(),Delta_P3_vs_P1=numeric(),F_statistic=numeric(),P_value=numeric(),FDR_global=numeric())
  if(nrow(bases)) bases$FDR_global <- p.adjust(bases$P_value,'BH')
  else bases <- data.frame(Branch=character(),Program=character(),Delta_P2_vs_P1=numeric(),
                           Delta_P3_vs_P1=numeric(),F_statistic=numeric(),P_value=numeric(),FDR_global=numeric())
  csv(tests,'stage_condition_interaction_tests.csv'); csv(bases,'CK_development_module_F_tests.csv')
  csv(bind(statuses),'stage_interaction_status.csv'); csv(bind(missing),'stage_interaction_missing_groups.csv')
  # 分支内阶段的2对2处理比较，缺组由test_scores记录，所有比较统一BH。
  if(nrow(allscores)) {
    aa <- allscores; aa$Celltype <- paste0('Endodermis_',aa$Branch,'_',aa$Period)
    st <- test_scores(aa)
    csv(st$results,'lineage_stage_treatment_tests.csv'); csv(st$status,'lineage_stage_treatment_status.csv')
  }
  pages <- split(programs,ceiling(seq_along(programs)/6))
  d$Period <- factor(d$Period,c('P1','P2','P3')); d$Condition <- factor(d$Condition,CONDITIONS)
  for(b in branches) for(i in seq_along(pages)) {
    q <- d[d$Branch==b & d$Program %in% pages[[i]],,drop=FALSE]
    if(!nrow(q)) next
    g <- ggplot2::ggplot(q,ggplot2::aes(Period,Score,color=Condition,group=Sample))+
      ggplot2::geom_line(alpha=.7)+ggplot2::geom_point(size=1.8)+
      ggplot2::facet_wrap(~Program,scales='free_y',ncol=3)+ggplot2::theme_bw()+
      ggplot2::labs(title=paste('Endodermis module profiles:',b),
                    subtitle='Each line is one biological sample; fixed within-branch CK cutoffs',y='Mean Seurat module score')
    save_plot(g,paste0('Figure12_',b,'_stage_profiles_page',i),13,7)
  }
  # 连续趋势也按分支分别展示，不把不同分支拟时间平均为一个伪拟时间。
  m <- trajectory$long
  m <- m[m$Projection_pass & !is.na(m$Period) & is.finite(m$Pseudotime01),,drop=FALSE]
  if(nrow(m)) {
    m$Bin <- cut(m$Pseudotime01,c(-Inf,seq(0,1,length.out=PSEUDOTIME_BINS+1)[-c(1,PSEUDOTIME_BINS+1)],Inf),labels=FALSE)
    groups <- split(seq_len(nrow(m)),interaction(m$Branch,m$Sample,m$Bin,drop=TRUE))
    bs <- bind(lapply(groups,function(ix) data.frame(Branch=m$Branch[ix[1]],
                                                     Sample=m$Sample[ix[1]],Condition=m$Condition[ix[1]],Bin=m$Bin[ix[1]],
                                                     Pseudotime01=mean(m$Pseudotime01[ix]),n_cells=length(ix),Program=programs,
                                                     Score=colMeans(scores[m$cell[ix],,drop=FALSE]))))
    csv(bs,'continuous_pseudotime_binned_module_scores.csv')
    for(b in branches) for(i in seq_along(pages)) {
      q <- bs[bs$Branch==b & bs$Program %in% pages[[i]] & bs$n_cells>=MIN_CELLS,,drop=FALSE]
      if(!nrow(q)) next
      g <- ggplot2::ggplot(q,ggplot2::aes(Pseudotime01,Score,color=Condition,group=Sample))+
        ggplot2::geom_line()+ggplot2::geom_point(size=1)+
        ggplot2::facet_wrap(~Program,scales='free_y',ncol=3)+ggplot2::theme_bw()+
        ggplot2::labs(title=paste('Continuous module dynamics:',b),
                      subtitle='Within-branch CK scale; sample-level bins; no cell-level P values')
      save_plot(g,paste0('Figure13_',b,'_continuous_dynamics_page',i),13,7)
    }
  }
  peak <- bind(lapply(branches,function(b) bind(lapply(programs,function(p) {
    a <- d[d$Branch==b & d$Program==p & d$Condition=='CK',,drop=FALSE]
    if(!all(table(factor(a$Period,c('P1','P2','P3')))==EXPECTED_REPLICATES))
      return(data.frame(Branch=b,Program=p,CK_peak_stage=NA_character_,Status='Insufficient_CK_stages'))
    means <- tapply(a$Score,a$Period,mean)
    data.frame(Branch=b,Program=p,CK_peak_stage=names(means)[which.max(means)],
               Status='Descriptive_peak_not_specificity')
  }))))
  csv(peak,'module_CK_peak_stage.csv')
}

compare_old_modules <- function(net) {
  path <- if(is.null(OLD_MODULE_FILE)) file.path(OUT_DIR,'Nature_network_modules.csv') else OLD_MODULE_FILE
  if(!file.exists(path)) return(invisible(NULL))
  old <- read_csv(path)
  if(!all(c('gene_name','Program') %in% names(old))) return(invisible(NULL))
  old <- old[!is.na(old$Program) & nzchar(old$Program),]
  new <- net$modules[!is.na(net$modules$Program),]
  result <- bind(lapply(sort(unique(old$Program)),function(op) bind(lapply(sort(unique(new$Program)),function(np) {
    a <- unique(old$gene_name[old$Program==op]); b <- unique(new$gene_name[new$Program==np])
    data.frame(Old_program=op,New_program=np,Old_size=length(a),New_size=length(b),
               Intersection=length(intersect(a,b)),Jaccard=length(intersect(a,b))/length(union(a,b)))
  }))))
  csv(result,'old_new_module_overlap_descriptive.csv')
}

main <- function() {
  if(!dir.exists(OUT_DIR) && !dir.create(OUT_DIR,recursive=TRUE)) stop('无法创建输出目录。')
  finished <- FALSE
  csv(data.frame(Status='Running',Time=as.character(Sys.time()),Version=PIPELINE_VERSION),'pipeline_status.csv')
  on.exit(if(!finished) csv(data.frame(Status='Stopped_before_completion',Time=as.character(Sys.time()),
                                       Version=PIPELINE_VERSION,Note='Inspect the console and current run.log; older same-prefix files may exist'),
                            'pipeline_status.csv'),add=TRUE)
  need(c('Seurat','SeuratObject','Matrix','hdWGCNA','WGCNA','ggplot2','patchwork',
         'digest','princurve','irlba','FNN','limma','UCell'))
  op <- options(Seurat.object.assay.version='v3')
  on.exit(options(op),add=TRUE)
  on.exit(capture.output(sessionInfo(),file=out('sessionInfo.txt')),add=TRUE)
  set.seed(SEED); WGCNA::enableWGCNAThreads(nThreads=N_THREADS)
  note('全部目标细胞类型构网；内皮层固定三个表达聚类P1/P2/P3；输出目录=',OUT_DIR)
  dat <- read_data()
  key <- digest::digest(list(PIPELINE_VERSION,normalizePath(MERGED_RDS),file.info(MERGED_RDS)[,c('size','mtime')],
                             dat$md,ROOT_CLUSTER,if(is.null(ROOT_CELLS_FILE)) NULL else read_csv(ROOT_CELLS_FILE),
                             EARLY_ANCHOR_TYPES,REFERENCE_MAX_CELLS_PER_CK,TRAJECTORY_N_HVG,TRAJECTORY_N_PCS,
                             SEED,PROJECTION_DISTANCE_QUANTILE,
                             PROJECTION_DISTANCE_MULTIPLIER,body(endo_trajectory),body(reference_pca),
                             body(fit_reference_pca),body(project_reference_pca),body(fit_three_clusters),body(nearest_centroid),body(project_fixed_path),
                             sapply(c('Seurat','irlba','princurve'),function(p)as.character(utils::packageVersion(p)))))
  trajectory <- cache('trajectory_cache',key,function() endo_trajectory(dat))
  dat <- apply_endo_stages(dat,trajectory)
  nkey <- digest::digest(list(key,dat$md,SEED,N_HVG,MIN_DETECT_FRACTION,EXCLUDE_GENES,
                              MIN_METACELL_GROUP,MAX_TRAIN_CELLS_PER_SAMPLE_CT,METACELL_K,METACELL_MAX_SHARED,TARGET_METACELLS,
                              MIN_MODULE_SIZE,SOFT_POWER_FIXED,SOFT_R2,SOFT_MIN_K,body(build_network),
                              sapply(c('hdWGCNA','WGCNA','Seurat'),function(p)as.character(utils::packageVersion(p)))))
  net <- cache('network_cache',nkey,function()build_network(dat))
  csv(net$modules,'network_modules.csv'); compare_old_modules(net)
  skey <- digest::digest(list(nkey,SEURAT_CTRL,SEURAT_NBIN,body(all_gene_seurat_scores),body(make_full_object)))
  scored <- cache('Seurat_score_object',skey,function()all_gene_seurat_scores(make_full_object(dat,net)))
  scores <- scored$seurat; x <- scored$object
  analysis_md <- dat$md[dat$md$Analysis_eligible,,drop=FALSE]
  samples <- aggregate_scores(scores,analysis_md); csv(samples,'sample_program_scores.csv')
  p1_module_evidence(samples,'new_modules_P1')
  check_previous_program2(x,dat$md)
  labels <- annotate_celltypes(samples[samples$Celltype!='Endodermis_Unassigned',])
  csv(labels$summary,'module_celltype_stage_labels.csv')
  csv(labels$evidence,'module_celltype_stage_label_evidence.csv')
  cklabels <- labels$summary[labels$summary$Condition=='CK',]
  csv(cklabels,'module_celltype_stage_labels_CK.csv')
  le <- labels$evidence
  if(nrow(le)) {
    le$Condition <- factor(le$Condition,CONDITIONS)
    g <- ggplot2::ggplot(le,ggplot2::aes(Celltype,Program,fill=Mean_Z))+
      ggplot2::geom_tile()+ggplot2::geom_point(data=le[le$Candidate,],shape=8,size=1)+
      ggplot2::facet_wrap(~Condition,ncol=2)+ggplot2::scale_fill_gradient2(low='#3565A0',mid='white',high='#B72C43')+
      ggplot2::theme_bw()+ggplot2::theme(axis.text.x=ggplot2::element_text(angle=60,hjust=1))+
      ggplot2::labs(x=NULL,y=NULL,title='Module preference across cell types and Endodermis stages',
                    subtitle='Sample-balanced descriptive labels; stars are not significance tests')
    save_plot(g,'Figure07_module_celltype_stage_heatmap',16,11)
  }
  tests <- test_scores(samples[samples$Celltype!='Endodermis_Unassigned',])
  csv(tests$results,'stage_and_celltype_treatment_tests.csv'); csv(tests$status,'comparison_status.csv')
  depth <- test_scores(samples[samples$Celltype!='Endodermis_Unassigned',],TRUE)
  csv(depth$results,'treatment_tests_depth_adjusted.csv'); csv(depth$status,'depth_test_status.csv')
  ukey <- digest::digest(list(skey,UCell_MAX_RANK,SCORE_CHUNK,body(score_modules)))
  ucell <- cache('UCell_scores',ukey,function()score_modules(dat,net))
  us <- aggregate_scores(ucell,analysis_md); ut <- test_scores(us[us$Celltype!='Endodermis_Unassigned',])
  csv(us,'UCell_sample_scores.csv'); csv(ut$results,'UCell_treatment_tests.csv')
  stage_dynamics(scores,dat$md,trajectory)
  plot_atlas(scores,dat$md,'Figure08_all_modules_Seurat','Seurat module score',cklabels)
  end <- dat$md[dat$md$Original_celltype==ENDO_CELLTYPE,]
  plot_atlas(scores,end,'Figure10_Endodermis_modules','Seurat module score',cklabels)
  x <- compute_mes_and_hubs(x,dat$md)
  x@misc[['Endodermis_three_cluster_trajectory']] <- trajectory
  saveRDS(x,out('all_cells_modules_scored.rds'))
  # 优先org.Gmax.eg.db，自动ID诊断；失败不会改阈值或阻断已经完成的网络结果。
  go <- if(RUN_GO) tryCatch(go_analysis(net$modules,net$background,csv),error=function(e) {
    csv(data.frame(Status='Error',Message=conditionMessage(e)),'GO_error.csv'); note('GO未完成：',conditionMessage(e)); NULL
  }) else NULL
  plot_go_compact(go)
  writeLines(c('Exactly three CK expression clusters; all eight-sample cell types retained.',
               'P1 is assigned using CK Meristem/QC expression proximity; P2/P3 use the shortest ordered centroid path.',
               'This is a fixed three-cluster analysis by user request, not evidence that development has exactly three natural states.',
               'A CK principal curve provides continuous ordering if converged and compatible with cluster order; otherwise a named centroid-path approximation is used.',
               'See trajectory_model_status.csv. No RNA velocity matrix was supplied; biological direction remains to be validated.',
               'All Endodermis cells receive one of P1/P2/P3; cells distant from CK are flagged separately and excluded from network/statistical training.',
               'Other cell types are unchanged. Metacells are grouped within sample x celltype/stage; independent sample datasets remain eight.',
               'The previous DevProgram02 gene set is scored on the new stages independently, without using it to define P1.',
               'New modules are StageProgramXX. They do not inherit old DevProgram numbers. See old_new_module_overlap_descriptive.csv.',
               'UMAP uses grey-to-module-color gradients with pooled 1%-99% display limits; numerical scores are unchanged.',
               'GO panels use larger fonts, wrapped terms, and two modules per page; enrichment thresholds are unchanged.',
               'Output directory unchanged; Endo3_v5_ files are new, old EndoDev_v4_ results and input RDS are preserved.',
               'Check pipeline_status, stage_interaction_status, ME_status, GO_status and trajectory_model_status.'),
             out('README.txt'),useBytes=TRUE)
  note('完成。检查三个聚类图、previous_Program02_P1_evidence、new_modules_P1_evidence、紧凑GO图及pipeline_status。')
  csv(data.frame(Status='Core_pipeline_complete',Time=as.character(Sys.time()),Version=PIPELINE_VERSION,
                 GO=if(!RUN_GO)'Disabled' else if(is.null(go))'Error_see_GO_error' else 'Completed',
                 Note='Inspect stage_interaction_status and ME_status for optional/coverage-dependent results'),
      'pipeline_status.csv')
  finished <- TRUE
  invisible(list(trajectory=trajectory,modules=net$modules,scores=scores,GO=go))
}
if(!isTRUE(getOption('soybean.dev.skip_main',FALSE))) main()

# Source: 图2.docx, extracted text lines 475-1718
# 独立运行；CK2补全RDS和外部汇总DEG表
# See README.md and docs/KNOWN_LIMITATIONS.md before running.

# 完整CK2输入版：原四分量PI流程不变，表达输入改为补回CK2的全细胞RDS。
# 仍只比较Heat/Cold/PEG vs CK，保留全部已注释细胞类型，不做TF或内皮层阶段分析。
# 原DEG文件保持不变；换表达RDS不会自动重算外部DEG。
# 原始输入RDS不覆盖；同目录同PI_UserDEG_前缀的结果会在本次运行时更新。
# 缓存由输入路径、文件信息、配置和函数签名验证，旧RDS缓存不会被误用。
# source本文件会自动运行main一次，不要在source后再次调用main。
# R Console：source(file.path(Sys.getenv("SOYBEAN_PROJECT_ROOT", unset = "data"), "Soybean_UserDEG_PI_Figure4_v3_CK2.R"), encoding="UTF-8")
# Terminal：Rscript /data01/liuxy/Project/Soybean_UserDEG_PI_Figure4_v3_CK2.R

# DEG分量只来自指定已筛选DEG；其余分量为Augur、scDist、Mixscape。
# 不更改原OUT_DIR，不运行DESeq2，不使用Endodermis状态表替代胁迫比较。
# 论文13分量 != 本脚本已有细胞类型4分量；不能复制同一DEG计数伪造四分辨率。
# 所有本次结果使用 PI_UserDEG_ 前缀，保存在原目录，旧PI缓存不被改写。
# 用法 source("Soybean_UserDEG_PI_Figure4_v3_CK2.R", encoding="UTF-8")
# =============================================================================
# ------------------------- 1. 配置 -------------------------------------------
MERGED_RDS <- file.path(Sys.getenv("SOYBEAN_PROJECT_ROOT", unset = "data"), "CK2_added_v3_Results/seurat_obj_with_CK2_v3.rds")
ANALYSIS_MODE <- "annotated" # 本版固定已有细胞类型/4项；不能用同一DEG计数伪造13项
# 固定到原脚本实际输出目录，不随MERGED_RDS的新位置移动。
OUT_DIR <- file.path(Sys.getenv("SOYBEAN_PROJECT_ROOT", unset = "data"), "CK2_added_v3_Results")
# 新版提前检查每个群的两个重复是否均有足够细胞；零细胞组合也写入 QC 表。
# 计算函数改变时使用独立内容签名缓存；无需更换输出目录。
RESUME <- TRUE
CELLTYPE_COL <- "cell_type"   # CK2完整RDS的已知注释列；保持全部细胞类型，不使用P1/P2/P3替代
CONDITION_COL <- "treatment"  # 本次报错修复：明确选择处理列，避免与 Group 冲突
SAMPLE_COL <- "sample_id"     # 你最新日志已确认此列通过检查；样本数仍会核验并导出
EXPECTED_REPLICATES <- 2L     # 你确认每个处理有两个重复；与元数据不一致时停止并报告
PAIR_COL <- NULL              # 有真实配对设计才填写；同一对跨处理共享此 ID
COVARIATE_COLS <- character() # 可选样本级协变量，如 c("batch")；需设计满秩
MIX_NEIGHBOR_BLOCK_COL <- NULL # 可选共同批次/配对区组，每个区组必须包含 CK
CELLTYPE_KEEP <- NULL         # 默认保留所有已注释类型
# 细胞类型重命名：把所有 lateral root 变体统一为 Procambium。
# 使用 canon() 规范化后做前缀匹配，可覆盖 "Lateral root"、"Lateral_root"、
# "Lateral root cells"、"Lateral roots" 等变体。
# 该映射同时作用于表达对象元数据与用户 DEG 表的 Celltype 列，保证两边对齐。
CELLTYPE_RENAME_MAP <- list(
  Procambium = c("Lateral root","Lateral Root","lateral root","Lateral_root","lateral_root",
                 "Lateral root cells","Lateral roots","lateral roots")
)
CONDITION_MAP <- list(
  CK=c("CK","Control","CTRL","REF","Reference","对照"),
  Heat=c("Heat","High_Temp","HighTemp","HighTemperature","HT","高温"),
  Cold=c("Cold","Low_Temp","LowTemp","LowTemperature","LT","低温"),
  PEG=c("PEG","PEG6000","PEG8000")
)
RESOLUTIONS <- c(0.1,0.2,0.4,0.8) # 原文分辨率；不强制复现水稻的群数
ASSAY <- "RNA"
N_HVG <- 3000L
N_PCS <- 30L
GRAPH_REDUCTION <- NULL       # NULL=新算 RNA PCA；或明确填写原 RDS 的 "harmony" 等
UMAP_REDUCTION <- "umap"      # 优先保留原图坐标；缺失时新算 UMAP
SEED <- 20261021L
N_THREADS <- 4L               # Augur 使用；若服务器内存吃紧可降到 2L
MIN_CELLS_PER_SAMPLE_CLUSTER <- 20L
MIN_REPLICATES_PER_CONDITION <- 2L # 至少 2 个，3 个以上更可靠
MIN_PB_LIBRARY <- 1000
AUGUR_N_SUBSAMPLES <- 50L
AUGUR_SUBSAMPLE_SIZE <- 20L
AUGUR_FOLDS <- 3L
SCDIST_D <- 20L
SCDIST_MIN_CELLS <- 20L
MIX_NEIGHBORS <- 20L
MIX_NDIMS <- 20L
MIX_MIN_DE_GENES <- 5L
MIX_LOGFC_THRESHOLD <- 0.25
MIX_ITERATIONS <- 10L
DENSE_MEMORY_BUDGET_GB <- 32 # 主要密集矩阵的估计预算，不是机器内存自动检测
SAVE_SEURAT_WITH_PI <- FALSE # TRUE 会额外保存重建的精简对象；不覆盖原 RDS
PIPELINE_VERSION <- "ensemble_pi_adapted_1.5_augur_future_fix"

# 旧计算配置保留以核验原缓存；下方新 main 不执行 DESeq2，不访问 DESeq2_tables。
# Augur 默认按细胞交叉验证，不是按生物学样本留一验证；PI 不是显著性或损伤概率。

# ------------------------- 2. 基础与检查 -------------------------------------
csv <- function(x,p) utils::write.csv(x,p,row.names=FALSE,na="",fileEncoding="UTF-8")
bind <- function(x) if(length(x)) do.call(rbind,x) else data.frame()
canon <- function(x) gsub("[[:space:]_.-]","",tolower(trimws(as.character(x))))
# 细胞类型重命名：把 lateral root 相关标签统一为 Procambium。
# 先用 canon() 规范化，再按前缀匹配，避免遗漏 "Lateral root cells" 之类的变体。
# 如需增删映射，只改 CELLTYPE_RENAME_MAP；函数本身无需修改。
apply_celltype_rename <- function(x, map=CELLTYPE_RENAME_MAP) {
  if(!length(map)) return(x)
  x_chr <- as.character(x)
  out <- x_chr
  cx <- canon(x_chr)
  for(target in names(map)) {
    aliases <- canon(map[[target]])
    aliases <- aliases[!is.na(aliases) & nzchar(aliases)]
    for(a in aliases) {
      hit <- !is.na(cx) & nzchar(cx) & (cx == a | startsWith(cx, a))
      out[hit] <- target
    }
  }
  out
}
map_condition <- function(x) {
  k<-unlist(lapply(CONDITION_MAP,canon),use.names=FALSE)
  v<-rep(names(CONDITION_MAP),lengths(CONDITION_MAP))
  if(any(vapply(split(v,k),function(a)length(unique(a))>1L,logical(1))))
    stop("CONDITION_MAP 有冲突别名。")
  unname(v[match(canon(x),k)])
}
pick_col <- function(md,explicit,candidates,label,required=TRUE,transform=as.character) {
  if(!is.null(explicit)) {
    if(length(explicit)!=1L || !explicit %in% names(md)) stop(label," 不存在：",explicit)
    return(explicit)
  }
  hit<-names(md)[canon(names(md)) %in% canon(candidates)]
  if(length(hit)>1L) {
    v<-lapply(hit,function(h)transform(md[[h]]))
    if(!all(vapply(v[-1],identical,logical(1),v[[1]])))
      stop(label," 候选列冲突：",paste(hit,collapse=", "),"。请在配置区指定。")
  }
  if(!length(hit)) {
    if(required) stop("无法识别 ",label,"；查看 00_metadata_inventory.csv 并指定。")
    return(NULL)
  }
  hit[1]
}
same_partition <- function(a,b) {
  a<-trimws(as.character(a));b<-trimws(as.character(b))
  identical(match(a,unique(a)),match(b,unique(b)))
}
resolve_sample_column <- function(md,condition,explicit=NULL) {
  # md 已限于本次分析的细胞；不让盐处理或未注释细胞影响样本列选择。
  if(!nrow(md)) stop("目标处理/细胞类型过滤后为空，不能检查样本。")
  candidates<-c("sample_id","sampleID","biological_sample","sample","sample_name")
  hits<-unique(unlist(lapply(candidates,function(k)names(md)[canon(names(md))==canon(k)])))
  if(!is.null(explicit)) {
    if(length(explicit)!=1L || !explicit %in% names(md)) stop("指定 SAMPLE_COL 不存在：",explicit)
    hits<-unique(c(explicit,hits))
  }
  # orig.ident 仅在没有标准样本列时作为备选，不能用它偷偷取代冲突的样本列。
  if(!length(hits) && "orig.ident" %in% names(md)) hits<-"orig.ident"
  if(!length(hits)) stop("没有样本候选列；请查看 00_metadata_inventory.csv。")
  cond<-as.character(condition);groups<-c("CK","Heat","Cold","PEG")
  cross<-data.frame(mapped_condition=cond,stringsAsFactors=FALSE)
  for(k in hits)cross[[k]]<-as.character(md[[k]])
  keys<-do.call(paste,c(lapply(cross,function(x)match(x,unique(x))),sep="/"))
  first<-!duplicated(keys);cross<-cross[first,,drop=FALSE]
  cross$n_cells<-tabulate(match(keys,keys[first]),nbins=sum(first))
  csv(cross,file.path(OUT_DIR,"00_sample_column_crosswalk.csv"))
  audits<-lapply(hits,function(k) {
    x<-trimws(as.character(md[[k]]));missing<-is.na(x)|!nzchar(x)
    pairs<-unique(data.frame(sample=x[!missing],condition=cond[!missing]))
    consistent<-!anyDuplicated(pairs$sample)
    ns<-vapply(groups,function(g)length(unique(x[!missing & cond==g])),integer(1))
    n<-length(unique(x[!missing]));avg<-if(n)sum(!missing)/n else 0
    bad<-c(if(any(missing))"missing_ID",if(!consistent)"same_ID_in_multiple_conditions",
           if(any(ns<MIN_REPLICATES_PER_CONDITION))"too_few_sample_IDs",
           if(avg<MIN_CELLS_PER_SAMPLE_CLUSTER)"too_many_IDs_or_too_few_cells")
    data.frame(column=k,n_sample_IDs=n,n_missing=sum(missing),
               CK=ns[1],Heat=ns[2],Cold=ns[3],PEG=ns[4],
               one_condition_per_ID=as.logical(consistent),mean_cells_per_ID=avg,
               structurally_eligible=!length(bad),reason=if(length(bad))paste(bad,collapse=";") else "passes_structure_checks",
               stringsAsFactors=FALSE,row.names=NULL)
  })
  audit<-bind(audits);csv(audit,file.path(OUT_DIR,"00_sample_column_audit.csv"))
  eligible<-audit$column[audit$structurally_eligible]
  decision<-function(chosen,why) {
    csv(data.frame(selected_sample_column=chosen,reason=why,
                   caveat="Structural consistency does not prove biological independence; technical libraries may still need merging."),
        file.path(OUT_DIR,"00_sample_column_decision.csv"))
    message("样本列：",chosen,"；",why)
    chosen
  }
  if(!is.null(explicit)) {
    if(!explicit %in% eligible) stop("指定的样本列未通过结构检查，见 00_sample_column_audit.csv。")
    return(decision(explicit,"user_explicit_selection"))
  }
  if(!length(eligible)) stop("样本候选列均未通过结构检查；见 00_sample_column_audit.csv 与 00_sample_column_crosswalk.csv。不能伪造生物学重复。")
  if(length(eligible)==1L) return(decision(eligible,"only_structurally_eligible_candidate"))
  equiv<-all(vapply(eligible[-1L],function(k)same_partition(md[[eligible[1L]]],md[[k]]),logical(1)))
  if(equiv) return(decision(eligible[1L],paste0("equivalent_cell_partition: ",paste(eligible,collapse=" = "))))
  csv(data.frame(selected_sample_column=NA_character_,reason="different_eligible_partitions_require_sample_design",
                 caveat="The metadata columns encode different sample groupings; choose using experimental sample records."),
      file.path(OUT_DIR,"00_sample_column_decision.csv"))
  stop("样本列确实把细胞分成不同组，不能仅按列名判断生物学重复。已导出 00_sample_column_crosswalk.csv 和 00_sample_column_audit.csv；请核对样本记录后设置 SAMPLE_COL。")
}
get_metadata <- function(obj) {
  old<-obj[[]]
  inv<-lapply(names(old),function(k) data.frame(column=k,n_unique=length(unique(old[[k]])),
                                                examples=paste(head(unique(as.character(old[[k]])),12),collapse=" | ")))
  csv(bind(inv),file.path(OUT_DIR,"00_metadata_inventory.csv"))
  ct<-pick_col(old,CELLTYPE_COL,c("celltype","cell_type","celltype_final","annotation","celltype_annotated","cell_type_annotation","seurat_annotations"),"CELLTYPE_COL")
  co<-pick_col(old,CONDITION_COL,c("condition","treatment","group","group_id","condition_id","stress"),
               "CONDITION_COL",transform=map_condition)
  mapped<-map_condition(old[[co]])
  au<-as.data.frame(table(as.character(old[[co]]),useNA="ifany"),stringsAsFactors=FALSE)
  names(au)<-c("original_condition","n_cells");au$mapped_condition<-map_condition(au$original_condition)
  csv(au,file.path(OUT_DIR,"00_condition_mapping.csv"))
  # 关键：先把 lateral root 变体重命名为 Procambium，再做 scope/过滤，确保后续所有计数按 Procambium 统计。
  types<-apply_celltype_rename(trimws(as.character(old[[ct]])))
  scope<-!is.na(mapped)&!is.na(types)&nzchar(types)&
    !tolower(types) %in% c("undefined","unknown","unassigned","unclassified","na")
  if(!is.null(CELLTYPE_KEEP))scope<-scope & types %in% CELLTYPE_KEEP
  if(!all(c("CK","Heat","Cold","PEG") %in% mapped[scope]))
    stop("必须识别 CK、高温、低温、PEG；请查看 00_condition_mapping.csv，其他处理不参与分析。")
  sa<-resolve_sample_column(old[scope,,drop=FALSE],mapped[scope],SAMPLE_COL)
  z<-data.frame(cell=rownames(old),pi_celltype=types,
                pi_condition=map_condition(old[[co]]),pi_sample=trimws(as.character(old[[sa]])),
                row.names=rownames(old),stringsAsFactors=FALSE)
  extra<-unique(c(PAIR_COL,COVARIATE_COLS,MIX_NEIGHBOR_BLOCK_COL))
  if(any(!extra %in% names(old))) stop("协变量/配对/区组列不存在：",paste(setdiff(extra,names(old)),collapse=", "))
  if(!is.null(PAIR_COL)) z$pi_pair<-as.character(old[[PAIR_COL]])
  if(length(COVARIATE_COLS)) for(i in seq_along(COVARIATE_COLS))
    z[[paste0("pi_cov",i)]]<-old[[COVARIATE_COLS[i]]]
  if(!is.null(MIX_NEIGHBOR_BLOCK_COL)) z$pi_mixblock<-as.character(old[[MIX_NEIGHBOR_BLOCK_COL]])
  keep<-!is.na(z$pi_condition)&!is.na(z$pi_celltype)&nzchar(z$pi_celltype)&
    !tolower(z$pi_celltype) %in% c("undefined","unknown","unassigned","unclassified","na")
  if(!is.null(CELLTYPE_KEEP)) keep<-keep & z$pi_celltype %in% CELLTYPE_KEEP
  z<-z[keep,,drop=FALSE]
  if(!nrow(z)||anyDuplicated(z$cell)) stop("细胞为空或条形码重复。")
  if(anyNA(z)||any(!nzchar(z$pi_sample))) stop("保留细胞中存在缺失 metadata；请先补齐相关列。")
  if(!all(c("CK","Heat","Cold","PEG") %in% z$pi_condition)) stop("必须识别 CK、Heat、Cold、PEG；检查映射表。")
  cols<-c("pi_sample","pi_condition",if(!is.null(PAIR_COL))"pi_pair",paste0("pi_cov",seq_along(COVARIATE_COLS)))
  # paste0 的零长度行为显式处理。
  cols<-intersect(cols,names(z)); sm<-unique(z[,cols,drop=FALSE])
  csv(sm,file.path(OUT_DIR,"00_sample_mapping.csv"))
  if(anyDuplicated(sm$pi_sample)) stop("一个样本 ID 对应多个处理/协变量/配对值；请修正样本信息。")
  if(nrow(sm)>nrow(z)/MIN_CELLS_PER_SAMPLE_CLUSTER) stop("样本 ID 过多；检查是否误用了单细胞条形码。")
  nr<-table(factor(sm$pi_condition,levels=c("CK","Heat","Cold","PEG")))
  csv(data.frame(condition=names(nr),n_samples=as.integer(nr)),file.path(OUT_DIR,"00_sample_counts.csv"))
  control_counts<-as.data.frame(table(Sample=z$pi_sample[z$pi_condition=="CK"]),stringsAsFactors=FALSE)
  names(control_counts)<-c("Sample","Cells")
  csv(control_counts,file.path(OUT_DIR,"00_CK_control_sample_counts.csv"))
  message("CK对照样本及细胞数：",paste(paste(control_counts$Sample,control_counts$Cells,sep="="),collapse=", "))
  message("全局样本数：",paste(paste(names(nr),as.integer(nr),sep="="),collapse=", "))
  if(!is.null(EXPECTED_REPLICATES) && any(nr!=EXPECTED_REPLICATES))
    stop("sample_id 中的样本数与每组两个重复不一致，见 00_sample_mapping.csv；不会自动合并或拆分样本。")
  if(any(nr<MIN_REPLICATES_PER_CONDITION))
    stop("至少一个条件缺少足够生物学重复，不能可靠生成完整 PI。不会用随机拆分细胞伪造重复。")
  csv(data.frame(role=c("celltype","condition","sample"),column=c(ct,co,sa)),
      file.path(OUT_DIR,"00_resolved_columns.csv"))
  z$pi_condition<-factor(z$pi_condition,levels=c("CK","Heat","Cold","PEG"))
  z
}

raw_counts <- function(obj,cells) {
  if(!ASSAY %in% names(obj@assays)) stop("缺少原始 RNA assay。")
  a<-obj[[ASSAY]]
  if(inherits(a,"Assay5")) {
    ls<-SeuratObject::Layers(a,search=NA);cl<-ls[grepl("^counts($|[.])",ls)]
    if(!length(cl)) stop("RNA 中没有 counts 层。")
    if(length(cl)>1) {
      cc<-unlist(lapply(cl,function(k)colnames(SeuratObject::LayerData(a,layer=k))))
      if(anyDuplicated(cc)) stop("多个 counts 层有重复细胞，不能直接合并。")
      a<-SeuratObject::JoinLayers(a,layers="counts",new="counts");cl<-"counts"
    }
    x<-SeuratObject::LayerData(a,layer=cl)
  } else x<-SeuratObject::GetAssayData(a,slot="counts")
  if(!all(cells %in% colnames(x))) stop("部分细胞缺少 RNA counts。")
  x<-x[,cells,drop=FALSE]
  if(!inherits(x,"sparseMatrix")) x<-methods::as(x,"dgCMatrix")
  if(is.null(rownames(x))||anyDuplicated(rownames(x))) stop("基因 ID 缺失或重复。")
  if(any(!is.finite(x@x))||any(x@x<0)||any(abs(x@x-round(x@x))>1e-6))
    stop("需要原始非负整数 counts，不能用 data/SCT/integrated 值。")
  if(any(Matrix::colSums(x)==0)) stop("存在总 counts 为 0 的细胞。")
  x[Matrix::rowSums(x)>0,,drop=FALSE]
}
resolution_key <- function(r) paste0("r",gsub("[.]","p",format(r,trim=TRUE,scientific=FALSE)))
analysis_groups <- function() {
  if(identical(ANALYSIS_MODE,"annotated"))
    return(data.frame(grouping="annotation",column="pi_celltype",resolution=NA_real_))
  if(!identical(ANALYSIS_MODE,"multi_resolution"))stop("未知 ANALYSIS_MODE。")
  rk<-vapply(RESOLUTIONS,resolution_key,character(1))
  data.frame(grouping=rk,column=paste0("pi_",rk),resolution=RESOLUTIONS)
}
memory_check <- function(ng,nc,multiplier=3) {
  # nrow/ncol 返回 integer；必须在相乘前转 double，避免 >2^31-1 时溢出为 NA。
  dims<-as.double(c(ng,nc,multiplier))
  if(length(dims)!=3L || any(!is.finite(dims)) || any(dims<0))stop("内存估计的维度无效。")
  gb<-prod(dims)*8/1024^3
  if(!is.finite(gb))stop("内存估计超过可表示范围。")
  if(gb>DENSE_MEMORY_BUDGET_GB) stop(sprintf("估计密集矩阵工作量 %.1f GB 超过配置预算 %.1f GB；请调整内存预算/服务器，不能静默删基因。",gb,DENSE_MEMORY_BUDGET_GB))
  invisible(gb)
}
get_expression <- function(obj,assay="RNA",layer="data") {
  if(utils::packageVersion("SeuratObject") >= "5.0.0")
    SeuratObject::GetAssayData(obj,assay=assay,layer=layer)
  else SeuratObject::GetAssayData(obj,assay=assay,slot=layer)
}
# ------------------------- 3. 预处理及共同聚类 -------------------------------
prepare <- function() {
  original<-readRDS(MERGED_RDS)
  if(!inherits(original,"Seurat")) stop("RDS 不是 Seurat 对象。")
  md<-get_metadata(original);counts<-raw_counts(original,md$cell)
  # 重建精简 RNA 对象，避免不明 integrated/data 层被当作表达输入。
  obj<-SeuratObject::CreateSeuratObject(counts=counts,meta.data=md,assay="RNA")
  obj<-Seurat::NormalizeData(obj,normalization.method="LogNormalize",scale.factor=10000,verbose=FALSE)
  obj<-Seurat::FindVariableFeatures(obj,selection.method="vst",nfeatures=N_HVG,verbose=FALSE)
  hv<-SeuratObject::VariableFeatures(obj)
  np<-min(N_PCS,length(hv)-1L,ncol(obj)-1L)
  if(np<5) stop("可用细胞/特征过少，无法稳定进行 PCA。")
  obj<-Seurat::ScaleData(obj,features=hv,verbose=FALSE)
  obj<-Seurat::RunPCA(obj,features=hv,npcs=np,seed.use=SEED,verbose=FALSE)
  red<-"pca"
  if(ANALYSIS_MODE=="multi_resolution" && !is.null(GRAPH_REDUCTION)) {
    if(!GRAPH_REDUCTION %in% names(original@reductions)) stop("GRAPH_REDUCTION 不在原对象中。")
    emb<-SeuratObject::Embeddings(original[[GRAPH_REDUCTION]])
    if(!all(colnames(obj) %in% rownames(emb))) stop("指定降维缺少部分分析细胞。")
    obj[["pi_graph_embedding"]]<-SeuratObject::CreateDimReducObject(
      embeddings=emb[colnames(obj),,drop=FALSE],key="PIG_",assay="RNA")
    red<-"pi_graph_embedding"
  }
  dims<-seq_len(min(N_PCS,ncol(SeuratObject::Embeddings(obj[[red]]))))
  groups<-analysis_groups();clcols<-groups$column
  if(ANALYSIS_MODE=="multi_resolution") {
    obj<-Seurat::FindNeighbors(obj,reduction=red,dims=dims,graph.name=c("pi_nn","pi_snn"),verbose=FALSE)
    for(r in RESOLUTIONS) {
      obj<-Seurat::FindClusters(obj,graph.name="pi_snn",resolution=r,algorithm=1,
                                random.seed=SEED,verbose=FALSE)
      col<-paste0("pi_",resolution_key(r));obj[[col]]<-as.character(SeuratObject::Idents(obj))
    }
  }
  if(!is.null(UMAP_REDUCTION) && UMAP_REDUCTION %in% names(original@reductions)) {
    em<-SeuratObject::Embeddings(original[[UMAP_REDUCTION]])
    if(!all(colnames(obj) %in% rownames(em))) stop("原 UMAP 缺少部分细胞，请设 UMAP_REDUCTION <- NULL 重算。")
    obj[["pi_umap"]]<-SeuratObject::CreateDimReducObject(embeddings=em[colnames(obj),1:2,drop=FALSE],key="PIUMAP_",assay="RNA")
  } else obj<-Seurat::RunUMAP(obj,reduction=red,dims=dims,reduction.name="pi_umap",seed.use=SEED,verbose=FALSE)
  md<-obj[[]];md$cell<-rownames(md)
  csv(md,file.path(OUT_DIR,"01_cell_metadata_clusters.csv"))
  csv(data.frame(groups,n_clusters=vapply(clcols,function(k)length(unique(md[[k]])),integer(1))),
      file.path(OUT_DIR,"01_cluster_counts.csv"))
  tab<-lapply(clcols,function(k){x<-as.data.frame(table(cluster=md[[k]],condition=md$pi_condition));x$resolution<-k;x})
  csv(bind(tab),file.path(OUT_DIR,"01_cluster_condition_cell_counts.csv"))
  csv(data.frame(Gene=hv),file.path(OUT_DIR,"01_PCA_variable_genes.csv"))
  rm(original);invisible(gc());obj
}

aggregate_cluster <- function(counts,md,col) {
  key<-paste(match(md$pi_sample,unique(md$pi_sample)),match(md[[col]],unique(md[[col]])),sep="/")
  lv<-unique(key);j<-match(key,lv)
  d<-Matrix::sparseMatrix(i=seq_along(j),j=j,x=1,dims=c(length(j),length(lv)))
  pb<-counts %*% d;m<-md[match(lv,key),,drop=FALSE]
  m$cluster<-as.character(m[[col]]);m$n_cells<-tabulate(j,nbins=length(lv))
  m$lib_size<-as.numeric(Matrix::colSums(pb));m$pb_id<-paste0("pb",seq_len(nrow(m)))
  colnames(pb)<-m$pb_id;rownames(m)<-m$pb_id
  m$pass<-m$n_cells>=MIN_CELLS_PER_SAMPLE_CLUSTER & m$lib_size>=MIN_PB_LIBRARY
  list(counts=pb,meta=m)
}
sample_design <- function(m,st) {
  m<-m[m$pi_condition %in% c("CK",st) & m$pass,,drop=FALSE]
  if(!is.null(PAIR_COL)) {
    if(anyDuplicated(paste(m$pi_pair,m$pi_condition,sep="/")))
      stop("同一 cluster 中一个配对单位/条件存在多个样本；请核对设计或合并技术重复。")
    pairs<-intersect(m$pi_pair[m$pi_condition=="CK"],m$pi_pair[m$pi_condition==st])
    m<-m[m$pi_pair %in% pairs,,drop=FALSE]
  }
  tb<-table(factor(m$pi_condition,levels=c("CK",st)))
  if(any(tb<MIN_REPLICATES_PER_CONDITION)) stop("群内 QC 后可用样本：CK=",tb[1],"；",st,"=",tb[2],
                                                "；每组需要 ",MIN_REPLICATES_PER_CONDITION,"。这不是全实验的重复数；查看 QC_sample_group_counts.csv。")
  m$pi_condition<-factor(as.character(m$pi_condition),levels=c("CK",st))
  if(!is.null(PAIR_COL)) m$pi_pair<-factor(m$pi_pair)
  cov<-if(length(COVARIATE_COLS))paste0("pi_cov",seq_along(COVARIATE_COLS)) else character()
  # 比较内常数协变量无须入模；变化但与处理共线时下面会明确报错。
  cov<-cov[vapply(cov,function(k)length(unique(m[[k]]))>1,logical(1))]
  for(k in cov) if(is.character(m[[k]])||is.factor(m[[k]])) m[[k]]<-factor(m[[k]])
  form<-stats::reformulate(c(if(!is.null(PAIR_COL))"pi_pair",cov,"pi_condition"))
  mm<-stats::model.matrix(form,m)
  if(qr(mm)$rank<ncol(mm) || nrow(mm)<=ncol(mm)) stop("设计不满秩或没有剩余自由度；检查批次与处理混杂。")
  list(meta=m,formula=form)
}
run_augur <- function(norm,md,col,st) {
  take<-md$pi_condition %in% c("CK",st);m<-droplevels(md[take,,drop=FALSE])
  tab<-table(as.character(m[[col]]),factor(m$pi_condition,levels=c("CK",st)))
  eligible<-rownames(tab)[apply(tab,1,min)>=AUGUR_SUBSAMPLE_SIZE]
  m<-m[as.character(m[[col]]) %in% eligible,,drop=FALSE]
  if(!nrow(m)) stop("Augur：所有群均未达到每处理最低细胞数。")
  x<-norm[,rownames(m),drop=FALSE];x<-x[Matrix::rowSums(x)>0,,drop=FALSE]
  message(sprintf("Augur %s: cells=%d, genes=%d, celltypes=%d, n_subsamples=%d, subsample_size=%d, folds=%d, threads=%d",
                  st, ncol(x), nrow(x), length(eligible),
                  AUGUR_N_SUBSAMPLES, AUGUR_SUBSAMPLE_SIZE, AUGUR_FOLDS, N_THREADS))
  # 关键修复：Augur::calculate_auc 基于 future。RStudio / RStudio Server 下
  # 如果全局 future plan 是 multiprocess（Linux 落成 multicore），fork 会不稳定，
  # 一旦 worker 崩溃，主 R 进程可能只打印 "Execution halted" 而无具体错误信息。
  # 这里统一改为 multisession（RStudio 安全），计算结束再恢复原 plan。
  if(!requireNamespace("future",quietly=TRUE)) stop("Augur 依赖 future 包，请安装。")
  old_plan<-future::plan()
  on.exit(future::plan(old_plan),add=TRUE)
  if(N_THREADS>1L) future::plan(future::multisession,workers=N_THREADS)
  else future::plan(future::sequential)
  # 参数按实际 Augur 版本自适应；老版本 Augur 没有 n_threads。
  auc_args<-names(formals(Augur::calculate_auc))
  call_args<-list(input=x,meta=m,cell_type_col=col,label_col="pi_condition",
                  n_subsamples=AUGUR_N_SUBSAMPLES,subsample_size=AUGUR_SUBSAMPLE_SIZE,
                  folds=AUGUR_FOLDS,min_cells=AUGUR_SUBSAMPLE_SIZE,show_progress=TRUE)
  if("n_threads" %in% auc_args) call_args$n_threads<-N_THREADS
  out<-tryCatch(
    do.call(Augur::calculate_auc,call_args),
    error=function(e){
      msg<-sprintf("Augur calculate_auc 失败（cells=%d, genes=%d, threads=%d）：%s",
                   ncol(x),nrow(x),N_THREADS,conditionMessage(e))
      message(msg)
      stop(msg)
    }
  )
  if(!all(c("cell_type","auc") %in% names(out$AUC))) stop("Augur 返回结构不符，需核对安装版本。")
  # 原始矩阵不重复存入 checkpoint。
  out$X<-NULL;out$feature_importance<-NULL
  rm(x);invisible(gc())
  list(value=data.frame(cluster=as.character(out$AUC$cell_type),value=as.numeric(out$AUC$auc)),raw=out)
}
scdist_metadata <- function(md,col,cluster,st) {
  m<-droplevels(md[as.character(md[[col]])==cluster & md$pi_condition %in% c("CK",st),,drop=FALSE])
  cc<-table(m$pi_sample);keep<-names(cc)[cc>=MIN_CELLS_PER_SAMPLE_CLUSTER]
  m<-m[m$pi_sample %in% keep,,drop=FALSE]
  if(!is.null(PAIR_COL)) {
    pairs<-intersect(m$pi_pair[m$pi_condition=="CK"],m$pi_pair[m$pi_condition==st])
    m<-m[m$pi_pair %in% pairs,,drop=FALSE]
  }
  u<-unique(m[,c("pi_sample","pi_condition")]);nr<-table(factor(u$pi_condition,levels=c("CK",st)))
  if(any(nr<MIN_REPLICATES_PER_CONDITION)) stop("scDist 群内细胞数过滤后可用样本：CK=",nr[1],"；",st,"=",nr[2],
                                                "；每组需要 ",MIN_REPLICATES_PER_CONDITION,"。查看 QC_sample_group_counts.csv。")
  m$pi_condition<-factor(as.character(m$pi_condition),levels=c("CK",st));m$pi_sample<-factor(m$pi_sample)
  cov<-if(length(COVARIATE_COLS))paste0("pi_cov",seq_along(COVARIATE_COLS)) else character()
  cov<-cov[vapply(cov,function(k)length(unique(m[[k]]))>1,logical(1))]
  for(k in cov) if(is.character(m[[k]])||is.factor(m[[k]]))m[[k]]<-factor(m[[k]])
  fx<-c("pi_condition",cov);re<-c("pi_sample",if(!is.null(PAIR_COL))"pi_pair")
  mm<-model.matrix(reformulate(fx),m)
  if(qr(mm)$rank<ncol(mm)) stop("scDist 固定效应设计不满秩。")
  list(meta=m,fixed=fx,random=re,nr=nr)
}
run_scdist <- function(norm,md,col,cluster,st) {
  design<-scdist_metadata(md,col,cluster,st)
  m<-design$meta;fx<-design$fixed;re<-design$random;nr<-design$nr
  x<-norm[,rownames(m),drop=FALSE];x<-x[Matrix::rowSums(x)>0,,drop=FALSE]
  if(min(dim(x))<=SCDIST_D) stop("scDist 维度不足，不自动改变 d。")
  memory_check(nrow(x),ncol(x))
  # 每个群单独调用，避免 scDist 一次把全部单细胞矩阵转 dense。
  m$pi_sccluster<-"cluster"
  out<-scDist::scDist(normalized_counts=as.matrix(x),meta.data=m,fixed.effects=fx,
                      random.effects=re,clusters="pi_sccluster",d=SCDIST_D,min.count.per.cell=SCDIST_MIN_CELLS)
  rr<-out$results
  if(nrow(rr)!=1L || !"Dist." %in% names(rr) || !is.finite(rr[["Dist."]][1]))
    stop("scDist 未返回有效距离，不能设为 0。")
  list(value=as.numeric(rr[["Dist."]][1]),results=rr,n_CK=nr[1],n_stress=nr[2])
}

# 每个群 x 每个真实样本的完整网格；零细胞仅用于 QC，不伪造零表达样本入模。
group_qc <- function(pb,md,col,grouping) {
  samples<-unique(md[,c("pi_sample","pi_condition"),drop=FALSE])
  samples$pi_sample<-as.character(samples$pi_sample);samples$pi_condition<-as.character(samples$pi_condition)
  cls<-sort(unique(as.character(md[[col]])))
  grid<-expand.grid(pi_sample=samples$pi_sample,cluster=cls,stringsAsFactors=FALSE)
  grid$pi_condition<-samples$pi_condition[match(grid$pi_sample,samples$pi_sample)]
  makekey<-function(s,c)paste(match(s,samples$pi_sample),match(c,cls),sep="/")
  ix<-match(makekey(grid$pi_sample,grid$cluster),makekey(pb$meta$pi_sample,pb$meta$cluster))
  grid$n_cells<-pb$meta$n_cells[ix];grid$lib_size<-pb$meta$lib_size[ix]
  grid$n_cells[is.na(ix)]<-0L;grid$lib_size[is.na(ix)]<-0
  grid$pass_cells<-grid$n_cells>=MIN_CELLS_PER_SAMPLE_CLUSTER
  grid$pass_library<-grid$lib_size>=MIN_PB_LIBRARY
  grid$pass_DESeq2<-grid$pass_cells & grid$pass_library
  grid$pass_scDist<-grid$pass_cells
  grid$reason<-vapply(seq_len(nrow(grid)),function(i) {
    if(grid$n_cells[i]==0L)return("no_cells_in_this_group")
    reasons<-c(if(!grid$pass_cells[i])"below_min_cells",if(!grid$pass_library[i])"below_min_pseudobulk_library")
    if(length(reasons))paste(reasons,collapse=";") else "passes_sample_QC"
  },character(1))
  grid$grouping<-grouping;grid$group_index<-match(grid$cluster,cls)
  grid$min_cells_threshold<-MIN_CELLS_PER_SAMPLE_CLUSTER;grid$min_library_threshold<-MIN_PB_LIBRARY
  checks<-list()
  for(cl in cls)for(st in c("Heat","Cold","PEG")) {
    a<-grid[grid$cluster==cl & grid$pi_condition %in% c("CK",st),,drop=FALSE]
    for(method in c("DESeq2","scDist","Augur")) {
      if(method=="DESeq2") {
        flag<-a$pass_DESeq2
        result<-tryCatch({sample_design(pb$meta[pb$meta$cluster==cl,,drop=FALSE],st);""},error=conditionMessage)
      } else if(method=="scDist") {
        flag<-a$pass_scDist
        result<-tryCatch({scdist_metadata(md,col,cl,st);""},error=conditionMessage)
      } else {
        flag<-a$n_cells>0
        nums<-vapply(c("CK",st),function(g)sum(a$n_cells[a$pi_condition==g]),numeric(1))
        result<-if(any(nums<AUGUR_SUBSAMPLE_SIZE))"Augur: insufficient total cells in one condition" else ""
      }
      checks[[length(checks)+1L]]<-data.frame(grouping=grouping,group_index=match(cl,cls),cluster=cl,
                                              condition=st,method=method,n_CK_global=sum(samples$pi_condition=="CK"),
                                              n_stress_global=sum(samples$pi_condition==st),
                                              n_CK_with_cells=sum(a$pi_condition=="CK" & a$n_cells>0),
                                              n_stress_with_cells=sum(a$pi_condition==st & a$n_cells>0),
                                              n_CK_pass_sample_QC=sum(a$pi_condition=="CK" & flag),
                                              n_stress_pass_sample_QC=sum(a$pi_condition==st & flag),
                                              cells_CK=sum(a$n_cells[a$pi_condition=="CK"]),cells_stress=sum(a$n_cells[a$pi_condition==st]),
                                              eligible=!nzchar(result),reason=if(nzchar(result))result else "eligible_before_model_fit",
                                              stringsAsFactors=FALSE)
    }
  }
  list(samples=grid,comparisons=bind(checks))
}
skipped_job <- function(key,reason) {
  message("跳过 ",key,": ",reason)
  list(ok=FALSE,value=NULL,error=reason,warnings=character(),key=key,status="skipped_QC")
}

# ------------------------- 4. 单细胞 Mixscape --------------------------------
# Seurat 在补算缺失基因时，按 Command 名称读取 Tools 中保存的 Neighbor。
# 某些版本二者名称不一致，得到 NULL 后 Indices(NULL) 报错。
# 只给本次 CalcPerturbSig 生成、且条形码/索引均已核对的原对象增加别名；不改邻居。
select_mixscape_neighbors <- function(tool_list,command,rep_cells,rep_controls,k,indices_fun) {
  expected<-make.names(paste0("PRTB_",names(rep_cells)))
  if(anyDuplicated(expected))stop("Mixscape 区组名称规范化后重复。")
  candidates<-names(tool_list)[grepl("CalcPerturbSig",names(tool_list),fixed=TRUE)]
  valid<-function(bundle) {
    if(!is.list(bundle)||!all(expected %in% names(bundle)))return(FALSE)
    all(vapply(seq_along(expected),function(i) {
      ix<-tryCatch(indices_fun(bundle[[expected[i]]]),error=function(e)NULL)
      !is.null(ix) && is.matrix(ix) && ncol(ix)==k &&
        !is.null(rownames(ix)) && !anyDuplicated(rownames(ix)) &&
        setequal(rownames(ix),rep_cells[[i]]) &&
        all(is.finite(ix)) && all(ix>=1 & ix<=length(rep_controls[[i]]) & ix==floor(ix))
    },logical(1)))
  }
  good<-candidates[vapply(tool_list[candidates],valid,logical(1))]
  if(!length(good))stop("没有找到与本次细胞、CK 池和邻居数一致的 Neighbor 缓存；不会调用 Indices(NULL)。")
  if(length(good)>1L && !all(vapply(good[-1L],function(x)identical(tool_list[[good[1L]]],tool_list[[x]]),logical(1))))
    stop("发现多个不一致的 Neighbor 缓存，不能猜选。")
  source<-if(command %in% good)command else good[1L]
  list(bundle=tool_list[[source]],source=source,command=command,aliased=!identical(source,command))
}
rebuild_mixscape_neighbors <- function(embedding,rep_cells,rep_controls,ndims,k,nn_helper) {
  # 缓存缺失时，用同一 PCA、CK 池、RANN 和 k 重建；与 CalcPerturbSig 的 NNHelper 调用一致。
  if(ncol(embedding)<ndims)stop("重建 Mixscape 邻居所需 PCA 维度不足。")
  names_out<-make.names(paste0("PRTB_",names(rep_cells)))
  if(anyDuplicated(names_out))stop("Mixscape 区组名称规范化后重复。")
  out<-lapply(seq_along(rep_cells),function(i) {
    if(length(rep_controls[[i]])<k)stop("重建邻居时 CK 细胞数不足。")
    if(!all(c(rep_cells[[i]],rep_controls[[i]]) %in% rownames(embedding)))stop("PCA 缺少所需细胞。")
    nn_helper(data=embedding[rep_controls[[i]],seq_len(ndims),drop=FALSE],
              query=embedding[rep_cells[[i]],seq_len(ndims),drop=FALSE],k=k,method="rann")
  })
  setNames(out,names_out)
}
repair_mixscape_neighbors <- function(sub,st,block) {
  commands<-names(sub@commands)
  candidates<-commands[grepl("CalcPerturbSig",commands,fixed=TRUE)]
  candidates<-candidates[vapply(candidates,function(k) {
    p<-sub@commands[[k]]@params
    identical(p$new.assay.name,"PRTB") && identical(p$assay,"RNA") && identical(p$reduction,"pca")
  },logical(1))]
  diagnostic<-list(command_names=commands,tool_names=names(sub@tools),
                   command_params=lapply(sub@commands,function(z)z@params),seurat_version=as.character(utils::packageVersion("Seurat")))
  saveRDS(diagnostic,file.path(OUT_DIR,paste0("Mixscape_",st,"_command_diagnostics.rds")))
  if(length(candidates)!=1L)stop("必须能唯一定位本次 CalcPerturbSig 命令，详见 Mixscape command diagnostics。")
  md<-sub[[]];cells<-colnames(sub)
  labs<-if(is.null(block))rep("rep1",length(cells)) else as.character(md[cells,block])
  reps<-unique(labs)
  rep_cells<-setNames(lapply(reps,function(x)cells[labs==x]),reps)
  rep_controls<-lapply(rep_cells,function(x)x[as.character(md[x,"pi_condition"])=="CK"])
  failure<-NULL
  selected<-tryCatch(select_mixscape_neighbors(sub@tools,candidates,rep_cells,rep_controls,MIX_NEIGHBORS,SeuratObject::Indices),
                     error=function(e){failure<<-conditionMessage(e);NULL})
  rebuilt<-FALSE
  if(is.null(selected)) {
    message("邻居缓存不可直接使用；用相同 PCA/CK 池重新查找邻居：",failure)
    helper<-getFromNamespace("NNHelper","Seurat")
    if(!all(c("data","query","k","method") %in% names(formals(helper))))
      stop("当前 Seurat NNHelper 接口不兼容，已保存命令诊断。")
    bundle<-rebuild_mixscape_neighbors(SeuratObject::Embeddings(sub[["pca"]]),rep_cells,rep_controls,MIX_NDIMS,MIX_NEIGHBORS,helper)
    selected<-select_mixscape_neighbors(setNames(list(bundle),candidates),candidates,rep_cells,rep_controls,MIX_NEIGHBORS,SeuratObject::Indices)
    rebuilt<-TRUE
  }
  # Tool<- 自动使用调用函数名，不能用它设置任意别名，因此在当前对象的 tools 槽中明确赋值。
  sub@tools[[selected$command]]<-selected$bundle
  csv(data.frame(condition=st,source_tool=selected$source,lookup_command=selected$command,
                 action=if(rebuilt)"rebuilt_same_PCA_CK_RANN" else if(selected$aliased)"added_validated_alias" else "already_aligned",
                 original_issue=if(is.null(failure))"" else failure,
                 n_blocks=length(reps),neighbors=MIX_NEIGHBORS),file.path(OUT_DIR,paste0("Mixscape_",st,"_neighbor_audit.csv")))
  message("Mixscape 邻居核验：",selected$source," -> ",selected$command)
  sub
}
run_mixscape <- function(obj,st) {
  cells<-rownames(obj[[]])[obj$pi_condition %in% c("CK",st)]
  sub<-subset(obj,cells=cells)
  sub$pi_condition<-factor(as.character(sub$pi_condition),levels=c("CK",st))
  if(sum(sub$pi_condition=="CK")<MIX_NEIGHBORS) stop("Mixscape 的 CK 数量不足。")
  block<-if(is.null(MIX_NEIGHBOR_BLOCK_COL)) NULL else "pi_mixblock"
  if(!is.null(block)) {
    tb<-table(sub[[block]][,1],factor(sub$pi_condition,levels=c("CK",st)))
    if(any(tb[,"CK"]<MIX_NEIGHBORS)) stop("Mixscape 每个邻居区组都必须有足够 CK；不能按处理特异 sample_id 分组。")
  }
  if(ncol(SeuratObject::Embeddings(sub[["pca"]]))<MIX_NDIMS) stop("Mixscape PCA 维度不足。")
  hv<-intersect(SeuratObject::VariableFeatures(sub),rownames(sub[["RNA"]]))
  if(!length(hv))stop("Mixscape 缺少 RNA 变量基因，不能自动改为所有基因。")
  # 初始 PRTB 只为 hv 计算，不能用全基因维度相乘造成溢出/错误高估。
  # 这是初始矩阵估算，不是后续补算基因及 Seurat 所有中间对象的峰值保证。
  memory_check(length(hv),ncol(sub),multiplier=4)
  sub<-Seurat::CalcPerturbSig(sub,assay="RNA",features=hv,slot="data",
                              gd.class="pi_condition",nt.cell.class="CK",reduction="pca",ndims=MIX_NDIMS,
                              num.neighbors=MIX_NEIGHBORS,split.by=block,new.assay.name="PRTB",verbose=TRUE)
  sub<-repair_mixscape_neighbors(sub,st,block)
  # 不按四种聚类重跑 Mixscape；每个 stress vs CK 只跑一次。
  # RunMixscape 内部会对需要的基因补算 PRTB 并缩放，所以不将负的 PRTB 再 LogNormalize。
  prior_tools<-mixscape_tool_names(sub@tools)
  if(length(prior_tools))sub@tools[prior_tools]<-NULL # 清理仅限当前子对象中的旧结果，防止误取历史pvec
  RunMixscape<-Seurat::RunMixscape
  sub<-RunMixscape(sub,assay="PRTB",slot="scale.data",labels="pi_condition",
                   nt.class.name="CK",new.class.name="pi_mixscape",min.de.genes=MIX_MIN_DE_GENES,
                   min.cells=5,de.assay="RNA",logfc.threshold=MIX_LOGFC_THRESHOLD,
                   iter.num=MIX_ITERATIONS,split.by=NULL,verbose=TRUE)
  sm<-sub[[]];cells<-colnames(sub)
  tools_found<-sub@tools[mixscape_tool_names(sub@tools)]
  # 先保存轻量原始结果，再提取；即使提取报错也不会丢失已算出的pvec。
  snapshot<-list(condition=st,cells=cells,metadata=sm[cells,"pi_condition",drop=FALSE],
                 tools=tools_found,all_tool_names=names(sub@tools),compute_signature=COMPUTE_SIGNATURE,
                 classification=sm[,grep("^pi_mixscape",names(sm)),drop=FALSE])
  saveRDS(snapshot,pi_path(paste0("Mixscape_",st,"_tool_snapshot.rds")))
  capture.output(str(tools_found,max.level=3),file=pi_path(paste0("Mixscape_",st,"_tool_structure.txt")))
  selected<-select_mixscape_pvec(tools_found,st,cells,as.character(sm[cells,"pi_condition"]))
  pi_csv(selected$audit,paste0("Mixscape_",st,"_tool_selection"))
  val<-selected$value
  message("Mixscape pvec来源：",selected$source,"；有效胁迫细胞：",
          sum(is.finite(val[as.character(sm[cells,"pi_condition"])==st])))
  list(value=val,tools=selected$tool,tool_source=selected$source,
       classification=sm[,grep("^pi_mixscape",names(sm)),drop=FALSE])
}


# ------------------------- 新版用户DEG输入与图形配置 --------------------------
PAPER_DIR <- file.path(Sys.getenv("SOYBEAN_PROJECT_ROOT", unset = "data"), "muscat_results_ck_ref_paper")
# 明确固定你的输入，不在找不到时换用其他DEG或旧DESeq2结果。
DEG_FILE <- file.path(PAPER_DIR,'All_DEG_filtered_paper_thresholds.csv')
# 如有意改用 combined.csv，请只修改上面的单一路径，不把两表相加。
DEG_COLUMNS <- list(gene=NULL,celltype=NULL,contrast=NULL,condition=NULL,state=NULL)
DEG_CONDITION_IS_VS_CK <- FALSE # 仅当表中单独Heat/Cold/PEG明确表示vs CK时才改TRUE
DEG_CELLTYPE_MAP <- character() # 可指定 c('DEG表标签'='Seurat实际细胞类型')；在细胞类型重命名之后再应用
ZERO_DEG_COMPARISONS <- NULL    # 可选已确认0 DEG的data.frame(Celltype,Condition)；默认未知不补零
EXPORT_TIFF <- FALSE
PLOT_POINT_SIZE <- 0.13
VIOLIN_POINT_MAX_PER_GROUP <- 500L # 只限制叠加小点；小提琴密度与统计使用全部完整PI细胞
OUTPUT_PREFIX <- 'PI_UserDEG_'
STRESSES <- c('Heat','Cold','PEG')
PI_COMPONENTS <- c('UserDEG_count','Augur_annotation','scDist_annotation','Mixscape')
USER_PI_VERSION <- 'user_deg_annotation_4_components_v3_CK2_input'
COMPUTE_SIGNATURE <- NULL
OLD_CACHE_VALID <- FALSE
PI_JOBS <- list()
CACHE_MANIFESTS <- list()
CURRENT_MANIFEST <- NULL
LEGACY_MANIFEST <- NULL

pi_path <- function(name) file.path(OUT_DIR,paste0(OUTPUT_PREFIX,name))
pi_csv <- function(d,name) csv(d,pi_path(paste0(name,'.csv')))
need <- function(d,cols,label) {
  if(!all(cols %in% names(d)))stop(label,'缺列：',paste(setdiff(cols,names(d)),collapse=', '))
}
unique_rows <- function(d,cols,label) {
  need(d,cols,label)
  if(anyNA(d[cols]) || anyDuplicated(d[cols]))stop(label,'存在NA或重复键：',paste(cols,collapse=' / '))
}
digest_object <- function(x) {
  f<-tempfile(fileext='.rds');on.exit(unlink(f));saveRDS(x,f,version=2)
  unname(tools::md5sum(f))
}
read_csv_pi <- function(path) {
  d<-read.csv(path,check.names=FALSE,stringsAsFactors=FALSE,fileEncoding='UTF-8-BOM')
  if(anyDuplicated(names(d)))stop('CSV存在重复列：',path)
  d
}
deg_column <- function(d,role,aliases,required=FALSE,multiple=FALSE) {
  explicit<-DEG_COLUMNS[[role]]
  if(!is.null(explicit)) {
    if(!all(explicit %in% names(d)))stop('DEG_COLUMNS$',role,'指定列不存在。')
    if(!multiple && length(explicit)!=1L)stop('该角色只能指定一列：',role)
    return(explicit)
  }
  hit<-names(d)[canon(names(d)) %in% canon(aliases)]
  if(!length(hit)) {
    if(required)stop('无法识别DEG的',role,'列；请填写DEG_COLUMNS。实际列：',paste(names(d),collapse=', '))
    return(NULL)
  }
  if(multiple)return(hit)
  if(length(hit)>1L) {
    v<-lapply(hit,function(n)trimws(as.character(d[[n]])))
    if(!all(vapply(v[-1],identical,logical(1),v[[1]])))
      stop('DEG列存在不同内容，请用DEG_COLUMNS$',role,'指定：',paste(hit,collapse=', '))
  }
  hit[1]
}
parse_ck_pair <- function(v) {
  # 接受明确的stress-vs-CK；反向比较也可用于计数，PI不使用logFC方向。
  x<-canon(v);out<-rep(NA_character_,length(v))
  for(st in STRESSES)for(a in CONDITION_MAP[[st]])for(b in CONDITION_MAP$CK) {
    labels<-c(paste(a,'vs',b),paste(b,'vs',a),paste(a,'versus',b),paste(b,'versus',a),
              paste(a,'minus',b),paste(b,'minus',a),paste(a,b,sep='-'),paste(b,a,sep='-'),
              paste(a,b,sep='/'),paste(b,a,sep='/'))
    out[!is.na(x) & x %in% canon(labels)]<-st
  }
  out
}
standardize_user_deg <- function(raw) {
  gc<-deg_column(raw,'gene',c('Gene','gene','gene_id','gene_name','geneid','feature','id'),TRUE)
  cc<-deg_column(raw,'celltype',c('Celltype','cell_type','celltype','cell_type_final','celltype_final','annotation','cluster_id'),TRUE)
  xc<-deg_column(raw,'contrast',c('contrast','comparison','contrast_name'),multiple=TRUE)
  tc<-deg_column(raw,'condition',c('Condition','treatment','stress','group'))
  sc<-deg_column(raw,'state',c('State','state','substate'))
  n<-nrow(raw)
  # 关键：DEG 的 Celltype 也应用与表达对象相同的 lateral root -> Procambium 重命名，
  # 否则 DEG 计数无法与表达对象的细胞类型对齐。
  d<-data.frame(Gene=trimws(as.character(raw[[gc]])),
                Celltype=apply_celltype_rename(trimws(as.character(raw[[cc]]))),
                Condition=NA_character_,State=if(is.null(sc))NA_character_ else trimws(as.character(raw[[sc]])),
                Status='Unresolved_comparison',Source_row=seq_len(n),stringsAsFactors=FALSE)
  if(anyNA(d$Gene) || any(!nzchar(d$Gene)))stop('DEG中存在空基因ID。')
  for(i in seq_len(n)) {
    labs<-if(length(xc))trimws(as.character(unlist(raw[i,xc,drop=FALSE]))) else character()
    labs<-labs[!is.na(labs) & nzchar(labs)]
    pars<-parse_ck_pair(labs)
    # 明确的状态比较/处理间比较不能被另一列的CK比较悄悄覆盖。
    unknown_pair<-is.na(pars) & grepl('vs|versus|minus|/| - ',labs,ignore.case=TRUE)
    valid<-unique(pars[!is.na(pars)])
    if(length(valid)>1L || any(unknown_pair)) {d$Status[i]<-'Conflicting_or_nonCK_comparison';next}
    if(length(valid)==1L) {d$Condition[i]<-valid;d$Status[i]<-'Explicit_CK_comparison';next}
    if(!length(labs) && !is.null(tc)) {
      lab<-as.character(raw[[tc]][i]);tr<-parse_ck_pair(lab)
      if(is.na(tr) && DEG_CONDITION_IS_VS_CK)tr<-map_condition(lab)
      if(!is.na(tr) && tr %in% STRESSES) {d$Condition[i]<-tr;d$Status[i]<-'Condition_column_CK_comparison'}
    }
  }
  d
}
load_user_deg <- function() {
  if(!file.exists(DEG_FILE))stop('指定DEG不存在：',DEG_FILE,'；不会自动更换输入。')
  # 即便使用已核验标准化缓存，也读取原表记录行数及实际列，不靠文件名猜测。
  raw<-read_csv_pi(DEG_FILE)
  pi_csv(data.frame(Column=names(raw)),'DEG_input_columns')
  md5<-unname(tools::md5sum(DEG_FILE));reader<-'Direct_user_DEG'
  saved<-file.path(OUT_DIR,'Nature_DE_filtered_annotated.csv')
  prov<-file.path(OUT_DIR,'Nature_DE_input_status.csv')
  d<-NULL
  # 用户指定列时直接解析，避免已有标准化文件遮蔽新配置。
  use_saved<-all(vapply(DEG_COLUMNS,is.null,logical(1))) && !DEG_CONDITION_IS_VS_CK
  if(use_saved && file.exists(saved) && file.exists(prov)) {
    audit<-read_csv_pi(prov)
    if(all(c('File','MD5') %in% names(audit))) {
      matchrow<-which(gsub('\\','/',audit$File,fixed=TRUE)==gsub('\\','/',DEG_FILE,fixed=TRUE) &
                        tolower(audit$MD5)==tolower(md5))
      if(length(matchrow)==1L) {
        q<-read_csv_pi(saved)
        required<-c('Gene','Celltype','State','Condition','Source_file','Source_row','Contrast_conflict')
        if(all(required %in% names(q)) && nrow(q)==nrow(raw) &&
           !anyNA(q$Source_file) && all(gsub('\\','/',q$Source_file,fixed=TRUE)==gsub('\\','/',DEG_FILE,fixed=TRUE))) {
          unique_rows(q,'Source_row','Nature标准化DEG')
          if(!setequal(q$Source_row,seq_len(nrow(raw))))stop('标准化DEG原始行号不完整。')
          d<-q[c('Gene','Celltype','State','Condition','Source_row')]
          # 标准化缓存里可能仍是旧标签（Lateral root），这里再套一次同样的重命名，保证与表达对象一致。
          d$Celltype<-apply_celltype_rename(d$Celltype)
          conflict<-as.character(q$Contrast_conflict)
          if(anyNA(conflict)||any(!conflict %in% c('TRUE','FALSE')))stop('标准化DEG冲突字段无效。')
          d$Status<-ifelse(conflict=='TRUE','Conflicting_or_nonCK_comparison',
                           ifelse(d$Condition %in% STRESSES,'Verified_standardized_user_DEG','Unresolved_comparison'))
          d$Condition[conflict=='TRUE']<-NA_character_
          reader<-'Verified_Nature_annotation_of_same_user_DEG'
        }
      }
    }
  }
  if(is.null(d))d<-standardize_user_deg(raw)
  # 状态级差异不能混作整细胞类型差异；不按名称推断状态与细胞类型的映射。
  state<-!is.na(d$State) & nzchar(trimws(d$State))
  d$Status[state]<-'State_level_comparison_excluded'
  d$Condition[state]<-NA_character_
  if(length(DEG_CELLTYPE_MAP)) {
    if(is.null(names(DEG_CELLTYPE_MAP)) || anyDuplicated(names(DEG_CELLTYPE_MAP)))stop('DEG_CELLTYPE_MAP需要唯一原标签。')
    k<-match(d$Celltype,names(DEG_CELLTYPE_MAP));yes<-!is.na(k)
    d$Celltype[yes]<-unname(DEG_CELLTYPE_MAP[k[yes]])
  }
  pi_csv(d,'DEG_row_audit')
  pi_csv(data.frame(File=DEG_FILE,MD5=md5,Rows=nrow(raw),Reader=reader,
                    Threshold='Already_filtered_by_user; no new FDR/logFC filter',
                    Excluded_state_rows=sum(state)),'DEG_source')
  if(!any(d$Condition %in% STRESSES))stop('你的DEG表没有可确认的整细胞类型胁迫vs CK比较；见DEG_row_audit。不会替换其他DE表。')
  list(data=d,md5=md5)
}
count_user_degs <- function(deg,celltypes) {
  # Gene ID原样去重；不移除大豆版本后缀、不混合转录本/基因ID。
  usable<-deg[deg$Condition %in% STRESSES & !is.na(deg$Celltype) & nzchar(deg$Celltype),,drop=FALSE]
  bad<-setdiff(unique(usable$Celltype),celltypes)
  pi_csv(usable[usable$Celltype %in% bad,,drop=FALSE],'DEG_unmatched_celltypes')
  grid<-expand.grid(Celltype=celltypes,Condition=STRESSES,stringsAsFactors=FALSE)
  grid$DEG_count<-NA_integer_;grid$Input_rows<-0L;grid$Status<-'Unknown_no_rows_in_filtered_file'
  for(i in seq_len(nrow(grid))) {
    z<-usable[usable$Celltype==grid$Celltype[i] & usable$Condition==grid$Condition[i],,drop=FALSE]
    if(nrow(z)) {
      grid$DEG_count[i]<-length(unique(z$Gene));grid$Input_rows[i]<-nrow(z)
      grid$Status[i]<-'Counted_unique_genes_in_user_DEG'
    }
  }
  if(!is.null(ZERO_DEG_COMPARISONS)) {
    unique_rows(ZERO_DEG_COMPARISONS,c('Celltype','Condition'),'显式零DEG比较')
    for(i in seq_len(nrow(ZERO_DEG_COMPARISONS))) {
      j<-which(grid$Celltype==ZERO_DEG_COMPARISONS$Celltype[i] & grid$Condition==ZERO_DEG_COMPARISONS$Condition[i])
      if(length(j)!=1L || is.finite(grid$DEG_count[j]))stop('显式零DEG记录不存在或与非零计数冲突。')
      grid$DEG_count[j]<-0L;grid$Status[j]<-'User_confirmed_zero_DEG'
    }
  }
  grid
}

# ------------------------- 可核验缓存：不改目录、不沿用旧DESeq2计数 ------------
compute_manifest <- function() {
  keys<-c('MERGED_RDS','ANALYSIS_MODE','EXPECTED_REPLICATES','CELLTYPE_COL','CONDITION_COL','SAMPLE_COL','PAIR_COL',
          'COVARIATE_COLS','MIX_NEIGHBOR_BLOCK_COL','CELLTYPE_KEEP','CELLTYPE_RENAME_MAP','CONDITION_MAP','RESOLUTIONS','ASSAY','N_HVG','N_PCS',
          'GRAPH_REDUCTION','UMAP_REDUCTION','SEED','N_THREADS','MIN_CELLS_PER_SAMPLE_CLUSTER','MIN_REPLICATES_PER_CONDITION',
          'MIN_PB_LIBRARY','AUGUR_N_SUBSAMPLES','AUGUR_SUBSAMPLE_SIZE','AUGUR_FOLDS','SCDIST_D','SCDIST_MIN_CELLS',
          'MIX_NEIGHBORS','MIX_NDIMS','MIX_MIN_DE_GENES','MIX_LOGFC_THRESHOLD','MIX_ITERATIONS','DENSE_MEMORY_BUDGET_GB')
  fns<-c('canon','apply_celltype_rename','map_condition','pick_col','same_partition','resolve_sample_column','get_metadata','raw_counts',
         'analysis_groups','resolution_key','prepare','get_expression','aggregate_cluster','group_qc','run_augur',
         'scdist_metadata','run_scdist','extract_mix_pvec','select_mixscape_neighbors',
         'rebuild_mixscape_neighbors','repair_mixscape_neighbors','run_mixscape',
         'mixscape_tool_names','select_mixscape_pvec')
  pkgs<-c('Seurat','SeuratObject','Matrix','Augur','scDist','mixtools')
  versions<-vapply(pkgs,function(x)if(requireNamespace(x,quietly=TRUE))as.character(utils::packageVersion(x)) else NA_character_,character(1))
  list(config=mget(keys,envir=environment(main)),input=file.info(MERGED_RDS)[,c('size','mtime'),drop=FALSE],
       packages=versions,functions=lapply(mget(fns,envir=environment(main)),function(f)paste(deparse(body(f)),collapse='\n')))
}
compatible_manifest <- function(old,new) {
  keys<-names(new$config);fn<-names(new$functions);pk<-names(new$packages)
  is.list(old) && all(keys %in% names(old$config)) && all(fn %in% names(old$functions)) &&
    all(pk %in% names(old$packages)) && !anyNA(new$packages) &&
    identical(old$config[keys],new$config) && identical(old$input,new$input) &&
    identical(old$functions[fn],new$functions) && identical(old$packages[pk],new$packages)
}
job_seed <- function(key) as.integer((SEED+sum(utf8ToInt(key)*seq_along(utf8ToInt(key))))%%.Machine$integer.max)
get_saved_components <- function() {
  paths<-file.path(OUT_DIR,c('04_cell_components_raw.csv','05_cell_PI.csv','01_UMAP_coordinates.csv'))
  if(!OLD_CACHE_VALID || !all(file.exists(paths)))return(NULL)
  raw<-read_csv_pi(paths[1]);md<-read_csv_pi(paths[2]);um<-read_csv_pi(paths[3])
  need(raw,c('cell',PI_COMPONENTS[-1]),'原非DEG分量')
  need(md,c('cell','pi_condition','pi_celltype','pi_sample'),'原细胞元数据')
  unique_rows(raw,'cell','原分量');unique_rows(md,'cell','原元数据');unique_rows(um,'cell','原UMAP')
  if(!setequal(raw$cell,md$cell) || !all(md$cell %in% um$cell))stop('原缓存细胞集合不一致。')
  raw<-raw[match(md$cell,raw$cell),c('cell',PI_COMPONENTS[-1]),drop=FALSE]
  if(!all(md$pi_condition %in% STRESSES))stop('原PI元数据含非胁迫条件。')
  # 旧05表中的PI本身完全不参与新计算。
  if(any(!is.finite(as.matrix(raw[PI_COMPONENTS[-1]]))))return(NULL)
  message('复用通过原manifest核验的Augur/scDist/Mixscape逐细胞数值；忽略旧DESeq2计数及旧PI。')
  list(cells=md,nondeg=raw,umap=um,source='Validated_existing_component_CSVs')
}
compute_other_components <- function() {
  deps<-c('Seurat','SeuratObject','Matrix','Augur','scDist','mixtools')
  miss<-deps[!vapply(deps,requireNamespace,logical(1),quietly=TRUE)]
  if(length(miss))stop('计算其他分量缺少包：',paste(miss,collapse=', '))
  dir.create(file.path(OUT_DIR,'checkpoints'),showWarnings=FALSE)
  prep<-user_checkpoint('prepared',prepare);register_pi(prep)
  if(!prep$ok)stop(prep$error)
  obj<-prep$value;rm(prep);invisible(gc())
  md<-obj[[]];norm<-get_expression(obj,'RNA','data');counts<-raw_counts(obj,rownames(md))
  cm<-md[md$pi_condition %in% STRESSES,,drop=FALSE];cm$cell<-rownames(cm)
  raw<-matrix(NA_real_,nrow(cm),3,dimnames=list(cm$cell,PI_COMPONENTS[-1]))
  pb<-aggregate_cluster(counts,md,'pi_celltype');qc<-group_qc(pb,md,'pi_celltype','annotation')
  pi_csv(qc$samples,'QC_sample_celltype');pi_csv(qc$comparisons,'QC_methods')
  clusters<-sort(unique(as.character(md$pi_celltype)))
  for(st in STRESSES) {
    check<-qc$comparisons[qc$comparisons$condition==st,,drop=FALSE]
    key<-paste('Augur','annotation',st,sep='_')
    z<-if(any(check$eligible[check$method=='Augur']))user_checkpoint(key,function()run_augur(norm,md,'pi_celltype',st)) else skipped_job(key,'No eligible celltypes')
    register_pi(z)
    take<-cm$pi_condition==st
    if(z$ok) {
      vals<-z$value$value;unique_rows(vals,'cluster','Augur')
      raw[take,'Augur_annotation']<-vals$value[match(cm$pi_celltype[take],vals$cluster)]
    }
    for(k in seq_along(clusters)) {
      cl<-clusters[k];take<-cm$pi_condition==st & cm$pi_celltype==cl
      if(!any(take))next
      key<-paste('scDist','annotation',st,paste0('cluster',k),sep='_')
      chk<-check[check$method=='scDist' & check$cluster==cl,,drop=FALSE]
      if(nrow(chk)!=1L)stop('scDist QC记录不唯一。')
      z<-if(chk$eligible)user_checkpoint(key,function()run_scdist(norm,md,'pi_celltype',cl,st)) else skipped_job(key,chk$reason)
      register_pi(z);if(z$ok)raw[take,'scDist_annotation']<-z$value$value
    }
    z<-user_checkpoint(paste0('Mixscape_',st),function()run_mixscape(obj,st));register_pi(z)
    take<-cm$pi_condition==st
    if(z$ok)raw[take,'Mixscape']<-z$value$value[cm$cell[take]]
  }
  em<-SeuratObject::Embeddings(obj[['pi_umap']])
  list(cells=cm,nondeg=data.frame(cell=cm$cell,raw,check.names=FALSE),
       umap=data.frame(cell=rownames(em),UMAP1=em[,1],UMAP2=em[,2]),source='Validated_checkpoints_or_recomputed_nonDEG_methods')
}

# ------------------------- 4分量合成与图4风格排版 -----------------------------
combine_user_pi <- function(raw) {
  raw<-as.matrix(raw)
  if(!identical(colnames(raw),PI_COMPONENTS))stop('PI分量名或顺序不正确。')
  if(any(raw[,'UserDEG_count']<0,na.rm=TRUE))stop('DEG计数不能为负。')
  if(any(raw[,'Augur_annotation']<0 | raw[,'Augur_annotation']>1,na.rm=TRUE))stop('Augur AUC应在[0,1]。')
  # Mixscape pvec是原始投影得分，不是posterior probability；可以小于0或大于1。
  # 不截断原始模型分量；在下面按各分量的实际范围统一min-max缩放。
  scaled<-raw;audits<-list()
  for(j in seq_len(ncol(raw))) {
    x<-raw[,j];ok<-is.finite(x);scaled[,j]<-NA_real_;lo<-hi<-NA_real_;status<-'No_values'
    if(any(ok)) {
      lo<-min(x[ok]);hi<-max(x[ok])
      if(lo==hi){scaled[ok,j]<-0;status<-'Constant_zero_convention_no_discriminating_information'}
      else {scaled[ok,j]<-(x[ok]-lo)/(hi-lo);status<-'Pooled_Heat_Cold_PEG_minmax'}
    }
    audits[[j]]<-data.frame(Component=colnames(raw)[j],Minimum=lo,Maximum=hi,
                            N_valid=sum(ok),N_missing=sum(!ok),Status=status)
  }
  n<-rowSums(is.finite(scaled));pi<-rep(NA_real_,nrow(raw));full<-n==4L
  pi[full]<-rowMeans(scaled[full,,drop=FALSE])
  list(PI=pi,scaled=scaled,n=n,audit=bind(audits))
}
sample_summary <- function(cells) {
  groups<-split(seq_len(nrow(cells)),interaction(cells$pi_sample,cells$pi_celltype,drop=TRUE))
  bind(lapply(groups,function(ii) {
    z<-cells[ii,,drop=FALSE];v<-z$PI[is.finite(z$PI)]
    data.frame(Sample=z$pi_sample[1],Condition=as.character(z$pi_condition[1]),Celltype=z$pi_celltype[1],
               N_cells=nrow(z),N_complete=length(v),Coverage=length(v)/nrow(z),
               Mean_PI=if(length(v))mean(v) else NA_real_)
  }))
}
pi_theme <- function(size=9) ggplot2::theme_classic(base_size=size)+ggplot2::theme(
  strip.background=ggplot2::element_rect(fill='#EFF2F3',color='#9AA2A8',linewidth=.3),
  strip.text=ggplot2::element_text(face='bold',size=size),plot.title=ggplot2::element_text(face='bold',size=size+1),
  plot.caption=ggplot2::element_text(size=7,hjust=0),legend.title=ggplot2::element_text(size=8),
  plot.margin=ggplot2::margin(5,5,5,5))
save_pi_plot <- function(p,name,w,h) {
  ggplot2::ggsave(pi_path(paste0(name,'.pdf')),p,device=if(capabilities('cairo'))grDevices::cairo_pdf else 'pdf',width=w,height=h,units='mm',limitsize=FALSE)
  ggplot2::ggsave(pi_path(paste0(name,'.png')),p,width=w,height=h,units='mm',dpi=300,limitsize=FALSE,bg='white')
  if(EXPORT_TIFF)ggplot2::ggsave(pi_path(paste0(name,'.tiff')),p,width=w,height=h,units='mm',dpi=600,compression='lzw',limitsize=FALSE)
}
plot_user_pi <- function(cm,scaled,um,samples) {
  need(um,c('cell','UMAP1','UMAP2'),'UMAP');unique_rows(um,'cell','UMAP')
  ix<-match(cm$cell,um$cell)
  if(anyNA(ix))stop('UMAP缺少PI细胞。')
  d<-cm;d$UMAP1<-um$UMAP1[ix];d$UMAP2<-um$UMAP2[ix]
  if(any(!is.finite(as.matrix(d[c('UMAP1','UMAP2')]))))stop('UMAP含非有限坐标。')
  d$pi_condition<-factor(d$pi_condition,levels=STRESSES)
  ct<-unique(as.character(d$pi_celltype));d$pi_celltype<-factor(d$pi_celltype,levels=ct)
  colors<-c(Heat='#D55E00',Cold='#4C87B8',PEG='#399B58')
  bg<-um[c('UMAP1','UMAP2')]
  labels<-c(UserDEG_count='User DEG count',Augur_annotation='Augur AUC',scDist_annotation='scDist distance',Mixscape='Mixscape')
  panels<-lapply(seq_along(PI_COMPONENTS),function(j) {
    x<-d[c('cell','UMAP1','UMAP2','pi_condition','pi_celltype')];x$Score<-scaled[,j]
    x$Method<-factor(labels[j],levels=unname(labels));x
  })
  a<-ggplot2::ggplot(bind(panels),ggplot2::aes(UMAP1,UMAP2))+
    ggplot2::geom_point(data=bg,color='#E6E8EA',size=PLOT_POINT_SIZE,inherit.aes=TRUE)+
    ggplot2::geom_point(ggplot2::aes(color=Score),size=PLOT_POINT_SIZE)+
    ggplot2::facet_grid(Method~pi_condition)+ggplot2::coord_equal()+
    ggplot2::scale_color_viridis_c(limits=c(0,1),option='D',na.value='#C8CCD0',name='Scaled\ncomponent')+
    pi_theme(8)+ggplot2::theme(axis.text=ggplot2::element_blank(),axis.ticks=ggplot2::element_blank(),
                               axis.line=ggplot2::element_blank(),legend.position='bottom')+
    ggplot2::labs(title='a  Components within annotated cell types',x=NULL,y=NULL)
  b<-ggplot2::ggplot(d,ggplot2::aes(UMAP1,UMAP2))+
    ggplot2::geom_point(data=bg,color='#E6E8EA',size=PLOT_POINT_SIZE)+
    ggplot2::geom_point(ggplot2::aes(color=PI),size=PLOT_POINT_SIZE)+ggplot2::coord_equal()+
    ggplot2::scale_color_viridis_c(limits=c(0,1),option='D',na.value='#C8CCD0',name='Perturbation\nindex')+
    pi_theme(8)+ggplot2::theme(axis.text=ggplot2::element_blank(),axis.ticks=ggplot2::element_blank(),
                               axis.line=ggplot2::element_blank(),legend.position='bottom')+
    ggplot2::labs(title='b  User-DEG ensemble PI',x=NULL,y=NULL)
  good<-d[is.finite(d$PI),,drop=FALSE]
  # 固定随机种子只用于小点抽稀；小提琴仍使用全部完整细胞。
  set.seed(SEED)
  groups<-split(seq_len(nrow(good)),interaction(good$pi_condition,good$pi_celltype,drop=TRUE))
  dots<-good[unlist(lapply(groups,function(ii)if(length(ii)>VIOLIN_POINT_MAX_PER_GROUP)sample(ii,VIOLIN_POINT_MAX_PER_GROUP) else ii)),,drop=FALSE]
  s<-samples[is.finite(samples$Mean_PI),,drop=FALSE]
  s$pi_condition<-factor(s$Condition,levels=STRESSES);s$pi_celltype<-factor(s$Celltype,levels=ct)
  cplot<-ggplot2::ggplot(good,ggplot2::aes(pi_condition,PI,fill=pi_condition))+
    ggplot2::geom_violin(scale='width',trim=TRUE,linewidth=.25,na.rm=TRUE)+
    ggplot2::geom_point(data=dots,position=ggplot2::position_jitter(width=.11,height=0,seed=SEED),size=.16,alpha=.12,color='black')+
    ggplot2::geom_point(data=s,ggplot2::aes(x=pi_condition,y=Mean_PI),inherit.aes=FALSE,
                        position=ggplot2::position_jitter(width=.09,height=0,seed=SEED),shape=21,fill='white',color='black',size=1.7,stroke=.4)+
    ggplot2::facet_wrap(~pi_celltype,ncol=4,drop=FALSE)+ggplot2::scale_fill_manual(values=colors,drop=FALSE)+
    ggplot2::scale_y_continuous(limits=c(0,1),breaks=c(0,.25,.5,.75,1),expand=ggplot2::expansion(mult=c(.01,.03)))+
    pi_theme(8)+ggplot2::theme(legend.position='none',panel.grid.major.y=ggplot2::element_line(color='#E6E8EA',linewidth=.2))+
    ggplot2::labs(title='c  Cell-type distributions and biological samples',x=NULL,y='Perturbation index',
                  caption='Violin: all complete cells; faint dots: display subsample; white dots: biological-sample means.\nFour components; missing components are not averaged. Distributions are descriptive, not independent biological replicates.')
  save_pi_plot(a,'A_components_UMAP',210,210)
  save_pi_plot(b+ggplot2::facet_wrap(~pi_condition,nrow=1),'B_PI_UMAP',210,85)
  ch<-max(90,ceiling(length(ct)/4)*43+25)
  save_pi_plot(cplot,'C_celltype_violin',230,ch)
  # grid原生组合，无patchwork/cowplot额外依赖。
  ga<-ggplot2::ggplotGrob(a);gb<-ggplot2::ggplotGrob(b+ggplot2::facet_wrap(~pi_condition,ncol=1));gc<-ggplot2::ggplotGrob(cplot)
  draw<-function() {
    grid::grid.newpage();grid::pushViewport(grid::viewport(layout=grid::grid.layout(2,2,
                                                                                    heights=grid::unit(c(150,ch),'mm'),widths=grid::unit(c(.66,.34),'null'))))
    for(z in list(list(ga,1,1),list(gb,1,2),list(gc,2,1:2))) {
      grid::pushViewport(grid::viewport(layout.pos.row=z[[2]],layout.pos.col=z[[3]]));grid::grid.draw(z[[1]]);grid::popViewport()
    }
    grid::popViewport()
  }
  h<-150+ch
  if(capabilities('cairo'))grDevices::cairo_pdf(pi_path('Figure4.pdf'),width=240/25.4,height=h/25.4)
  else grDevices::pdf(pi_path('Figure4.pdf'),width=240/25.4,height=h/25.4)
  tryCatch(draw(),finally=grDevices::dev.off())
  grDevices::png(pi_path('Figure4.png'),width=240,height=h,units='mm',res=300,bg='white')
  tryCatch(draw(),finally=grDevices::dev.off())
  pi_csv(d,'Figure_source_cells');pi_csv(samples,'Figure_source_samples')
}

# ------------------------- v3: Mixscape结果提取与可核验续算 ------------------
# SeuratObject::`Tool<-` 根据调用栈命名 tools；Seurat::RunMixscape 调用可能
# 保存为 Seurat::RunMixscape / Seurat..RunMixscape，而非固定 RunMixscape。
# 仅读取本次 RunMixscape 的 pvec，不用后验概率或其他分量代替。
mixscape_tool_names <- function(x) {
  n<-names(x)
  n[grepl('(^|[^[:alnum:]])RunMixscape($|[^[:alnum:]])',n)]
}
extract_mix_pvec <- function(tool_result,st,cells) {
  if(anyNA(cells) || anyDuplicated(cells))stop('Mixscape输入条形码不唯一。')
  if(!is.list(tool_result) || is.null(tool_result[[st]]))
    stop('Mixscape结果缺少处理键 ',st,'；实际键：',paste(names(tool_result),collapse=', '),
         '。请查看本次 tool_snapshot.rds；不能据此直接断言差异基因不足。')
  parts<-tool_result[[st]]
  if(is.data.frame(parts))parts<-list(parts)
  if(!is.list(parts) || !length(parts))stop('Mixscape处理结果为空：',st)
  frames<-lapply(parts,function(d) {
    if(!is.data.frame(d) || !'pvec' %in% names(d) || !is.numeric(d$pvec))
      stop('Mixscape结果不是含数值pvec的数据框；不会改用posterior probability。')
    bc<-rownames(d)
    if(is.null(bc) || anyNA(bc) || anyDuplicated(bc) || any(!nzchar(bc)))
      stop('Mixscape pvec缺少唯一细胞条形码。')
    if(!all(bc %in% cells))stop('Mixscape pvec含不属于本次比较的条形码；不自动删除前后缀。')
    data.frame(cell=bc,value=d$pvec,stringsAsFactors=FALSE)
  })
  z<-bind(unname(frames))
  if(anyDuplicated(z$cell))stop('Mixscape不同结果块有重复细胞；不能自动求平均。')
  val<-z$value[match(cells,z$cell)];names(val)<-cells
  if(!any(is.finite(val)))stop('Mixscape pvec没有有限数值。')
  val
}
select_mixscape_pvec <- function(tool_list,st,cells,condition) {
  if(length(condition)!=length(cells))stop('Mixscape细胞与处理标签长度不一致。')
  candidates<-mixscape_tool_names(tool_list)
  if(!length(candidates))stop('本次没有找到RunMixscape工具槽。实际槽：',paste(names(tool_list),collapse=', '))
  attempts<-lapply(candidates,function(k)tryCatch({
    v<-extract_mix_pvec(tool_list[[k]],st,cells)
    parts<-tool_list[[k]][[st]]
    if(is.data.frame(parts))parts<-list(parts)
    for(d in parts)if('pi_condition' %in% names(d)) {
      actual<-as.character(condition[match(rownames(d),cells)])
      if(anyNA(actual) || !identical(as.character(d$pi_condition),actual))
        stop('pvec数据框处理标签与本次元数据不一致。')
    }
    if(!any(is.finite(v[as.character(condition)==st])))stop('pvec没有胁迫细胞的有效值。')
    list(ok=TRUE,value=v,error='')
  },error=function(e)list(ok=FALSE,value=NULL,error=conditionMessage(e))))
  audit<-data.frame(Tool=candidates,OK=vapply(attempts,`[[`,logical(1),'ok'),
                    Error=vapply(attempts,`[[`,character(1),'error'),stringsAsFactors=FALSE)
  good<-which(audit$OK)
  if(!length(good))stop('RunMixscape槽读取失败：',paste(paste(audit$Tool,audit$Error,sep=': '),collapse=' | '))
  if(length(good)>1L && !all(vapply(good[-1L],function(i)
    identical(attempts[[good[1L]]]$value,attempts[[i]]$value),logical(1))))
    stop('多个RunMixscape槽含不同pvec，拒绝猜选。')
  chosen<-if('RunMixscape' %in% candidates[good])match('RunMixscape',candidates) else good[1L]
  list(value=attempts[[chosen]]$value,tool=tool_list[[candidates[chosen]]],source=candidates[chosen],audit=audit)
}

# 缓存迁移只允许非Mixscape任务忽略本次Mixscape函数的变化。
# 其余输入、全部配置、包版本和全部已有计算函数仍必须完全一致。
mixscape_function_names <- function() c('extract_mix_pvec','run_mixscape',
                                        'mixscape_tool_names','select_mixscape_pvec','select_mixscape_neighbors',
                                        'rebuild_mixscape_neighbors','repair_mixscape_neighbors')
reusable_manifest <- function(old,new,key) {
  if(grepl('^Mixscape_',key))return(compatible_manifest(old,new))
  if(!is.list(old) || !is.list(old$functions))return(FALSE)
  remove<-mixscape_function_names()
  old$functions<-old$functions[setdiff(names(old$functions),remove)]
  new$functions<-new$functions[setdiff(names(new$functions),remove)]
  compatible_manifest(old,new)
}
initialize_pi_cache <- function(manifest) {
  checkpoint_dir<-file.path(OUT_DIR,'checkpoints')
  dir.create(checkpoint_dir,showWarnings=FALSE,recursive=TRUE)
  paths<-unique(c(pi_path('compute_manifest.rds'),
                  list.files(checkpoint_dir,pattern='^UserPI_manifest_.*[.]rds$',full.names=TRUE)))
  found<-list()
  for(p in paths[file.exists(paths)]) {
    m<-tryCatch(readRDS(p),error=function(e)NULL)
    if(is.list(m) && !is.null(m$functions))found[[digest_object(m)]]<-m
  }
  # 先保留上一轮manifest，再更新入口文件，否则下一次无法验证旧成功缓存。
  for(sig in names(found)) {
    dest<-file.path(checkpoint_dir,paste0('UserPI_manifest_',sig,'.rds'))
    if(!file.exists(dest))saveRDS(found[[sig]],dest)
  }
  found[[COMPUTE_SIGNATURE]]<-manifest
  CACHE_MANIFESTS<<-found;CURRENT_MANIFEST<<-manifest
  oldfile<-file.path(OUT_DIR,'run_manifest.rds')
  LEGACY_MANIFEST<<-if(file.exists(oldfile))tryCatch(readRDS(oldfile),error=function(e)NULL) else NULL
  OLD_CACHE_VALID<<-compatible_manifest(LEGACY_MANIFEST,manifest)
  saveRDS(manifest,file.path(checkpoint_dir,paste0('UserPI_manifest_',COMPUTE_SIGNATURE,'.rds')))
  saveRDS(manifest,pi_path('compute_manifest.rds'))
  pi_csv(data.frame(Signature=names(found),
                    Reuse_nonMixscape=vapply(found,reusable_manifest,logical(1),new=manifest,key='prepared'),
                    Reuse_Mixscape=vapply(found,reusable_manifest,logical(1),new=manifest,key='Mixscape_Heat')),'cache_validation')
}
validate_pi_job <- function(key,value) {
  if(key=='prepared')return(invisible(TRUE))
  if(!is.list(value))stop('任务没有返回预期的list：',key)
  if(grepl('^Augur_',key)) {
    v<-value$value;unique_rows(v,'cluster','Augur缓存')
    if(!is.numeric(v$value) || !any(is.finite(v$value)) ||
       any(v$value[is.finite(v$value)]<0 | v$value[is.finite(v$value)]>1))
      stop('Augur未返回有效AUC。')
  } else if(grepl('^scDist_',key)) {
    if(!is.numeric(value$value) || length(value$value)!=1L || !is.finite(value$value))
      stop('scDist未返回单个有效距离。')
  } else if(grepl('^Mixscape_',key)) {
    v<-value$value
    if(!is.numeric(v) || is.null(names(v)) || anyDuplicated(names(v)) || !any(is.finite(v)))
      stop('Mixscape未返回含唯一条形码的有效pvec。')
  }
  invisible(TRUE)
}
user_checkpoint <- function(key,fun) {
  fresh<-file.path(OUT_DIR,'checkpoints',paste0('UserPI_',COMPUTE_SIGNATURE,'_',key,'.rds'))
  candidates<-character()
  if(RESUME) {
    candidates<-fresh
    for(sig in names(CACHE_MANIFESTS))if(reusable_manifest(CACHE_MANIFESTS[[sig]],CURRENT_MANIFEST,key))
      candidates<-c(candidates,file.path(OUT_DIR,'checkpoints',paste0('UserPI_',sig,'_',key,'.rds')))
    if(reusable_manifest(LEGACY_MANIFEST,CURRENT_MANIFEST,key))
      candidates<-c(candidates,file.path(OUT_DIR,'checkpoints',paste0(key,'.rds')))
  }
  for(path in unique(candidates))if(file.exists(path)) {
    z<-tryCatch(readRDS(path),error=function(e)NULL)
    if(is.list(z) && isTRUE(z$ok) && identical(z$key,key)) {
      issue<-tryCatch({validate_pi_job(key,z$value);NULL},error=function(e)conditionMessage(e))
      if(is.null(issue)) {z$cache_source<-path;message('复用成功缓存：',key);return(z)}
      message('缓存结果无效，将重新计算 ',key,'：',issue)
    }
  }
  set.seed(job_seed(key));message('计算：',key);warns<-character();error_call<-''
  z<-tryCatch(withCallingHandlers({
    v<-fun();validate_pi_job(key,v);list(ok=TRUE,value=v,error='')
  },warning=function(w) {warns<<-c(warns,conditionMessage(w));invokeRestart('muffleWarning')}),
  error=function(e) {
    error_call<<-paste(deparse(conditionCall(e)),collapse=' ')
    list(ok=FALSE,value=NULL,error=conditionMessage(e))
  })
  z$key<-key;z$warnings<-unique(warns);z$cache_source<-'';z$error_call<-error_call
  saveRDS(z,fresh)
  if(!z$ok)message('失败 [',key,']：',z$error)
  z
}
register_pi <- function(z) {
  PI_JOBS[[length(PI_JOBS)+1L]]<<-data.frame(Job=z$key,OK=z$ok,Error=z$error,
                                             Error_call=if(is.null(z$error_call))'' else z$error_call,
                                             Warnings=paste(z$warnings,collapse=' | '),Cache=if(is.null(z$cache_source))'' else z$cache_source)
  pi_csv(bind(PI_JOBS),'jobs')
}
write_pi_component_diagnostics <- function(cm,raw) {
  grouped<-split(seq_len(nrow(cm)),interaction(cm$pi_condition,cm$pi_celltype,drop=TRUE))
  tab<-bind(lapply(grouped,function(ii)bind(lapply(seq_len(ncol(raw)),function(j)
    data.frame(Condition=as.character(cm$pi_condition[ii[1]]),Celltype=as.character(cm$pi_celltype[ii[1]]),
               Component=colnames(raw)[j],N_cells=length(ii),N_valid=sum(is.finite(raw[ii,j])),
               N_missing=sum(!is.finite(raw[ii,j])))))))
  pi_csv(tab,'component_coverage')
  missing<-apply(!is.finite(raw),1,function(x)paste(colnames(raw)[x],collapse=';'))
  pi_csv(data.frame(cell=cm$cell,Condition=cm$pi_condition,Celltype=cm$pi_celltype,
                    Missing_components=missing),'missing_components')
  for(j in seq_len(ncol(raw)))message('有效分量 ',colnames(raw)[j],'：',sum(is.finite(raw[,j])),'/',nrow(raw))
  invisible(tab)
}


main <- function() {
  if(ANALYSIS_MODE!='annotated')stop('用户细胞类型DEG只能用于annotated四分量。没有四种分辨率对应DEG时，不能生成论文13分量。')
  if(!dir.exists(OUT_DIR))dir.create(OUT_DIR,recursive=TRUE,showWarnings=FALSE)
  if(!file.exists(MERGED_RDS))stop('原表达RDS不存在：',MERGED_RDS)
  if(!requireNamespace('ggplot2',quietly=TRUE))stop('缺少ggplot2。')
  status<-'Running';stage<-'DEG_input'
  pi_csv(data.frame(Status=status,Stage=stage,Version=USER_PI_VERSION,Time=as.character(Sys.time())),'run_status')
  on.exit({
    if(status=='Running')status<-'Error_or_interrupted'
    pi_csv(data.frame(Status=status,Stage=stage,Version=USER_PI_VERSION,Time=as.character(Sys.time())),'run_status')
    capture.output(sessionInfo(),file=pi_path('sessionInfo.txt'))
  },add=TRUE)
  message('R 版本：', R.version.string)
  message('表达输入（含补回CK2）：',MERGED_RDS)
  message('输出目录不变：',OUT_DIR)
  message('DEG计数唯一来源：',DEG_FILE)
  message('细胞类型重命名：Lateral root 相关标签 -> Procambium')
  message('Augur 并行：', N_THREADS, ' threads (future multisession)')
  dg<-load_user_deg()
  stage<-'Other_components'
  manifest<-compute_manifest();COMPUTE_SIGNATURE<<-digest_object(manifest)
  initialize_pi_cache(manifest)
  PI_JOBS<<-list()
  dat<-get_saved_components()
  if(is.null(dat))dat<-compute_other_components()
  cm<-dat$cells;unique_rows(cm,'cell','本次细胞元数据')
  unique_rows(unique(cm[c('pi_sample','pi_condition')]),'pi_sample','样本处理映射')
  counts<-count_user_degs(dg$data,sort(unique(as.character(cm$pi_celltype))))
  pi_csv(counts,'DEG_counts_by_celltype_condition')
  if(!any(is.finite(counts$DEG_count)))stop('DEG细胞类型与表达对象无法匹配；见DEG_unmatched_celltypes。')
  key<-function(ct,tr)paste(ct,tr,sep='\034')
  j<-match(key(cm$pi_celltype,cm$pi_condition),key(counts$Celltype,counts$Condition))
  nr<-dat$nondeg;unique_rows(nr,'cell','非DEG分量')
  ix<-match(cm$cell,nr$cell);if(anyNA(ix))stop('非DEG分量缺细胞。')
  raw<-cbind(UserDEG_count=counts$DEG_count[j],as.matrix(nr[ix,PI_COMPONENTS[-1],drop=FALSE]))
  rownames(raw)<-cm$cell
  write_pi_component_diagnostics(cm,raw)
  stage<-'Combine_PI';sc<-combine_user_pi(raw)
  cm$PI<-sc$PI;cm$n_components<-sc$n;cm$expected_components<-4L
  cm$analysis_mode<-'User_DEG_annotation_4_components'
  cm$PI_status<-ifelse(sc$n==4,'complete_4','missing_component')
  ss<-sample_summary(cm)
  pi_csv(data.frame(cell=cm$cell,raw,check.names=FALSE),'cell_components_raw')
  pi_csv(data.frame(cell=cm$cell,sc$scaled,check.names=FALSE),'cell_components_scaled')
  pi_csv(sc$audit,'scaling_ranges');pi_csv(cm,'cell_PI');pi_csv(ss,'sample_celltype_PI')
  coverage<-aggregate(cbind(N_cells,N_complete)~Condition+Celltype,ss,sum)
  coverage$Coverage<-coverage$N_complete/coverage$N_cells;pi_csv(coverage,'coverage')
  saveRDS(list(version=USER_PI_VERSION,DEG_file=DEG_FILE,DEG_md5=dg$md5,other_source=dat$source,
               compute_signature=COMPUTE_SIGNATURE,cells=cm,raw=raw,scaled=sc$scaled,scaling=sc$audit,samples=ss,umap=dat$umap),pi_path('results.rds'))
  writeLines(c('User-DEG四分量PI：New Phytologist Fig.4结构的明确适配版。',
               paste('唯一DEG：',DEG_FILE),paste('原输出目录：',OUT_DIR),
               '细胞类型重命名：Lateral root 及其变体在表达对象和DEG表中统一改为 Procambium。',
               'DEG数量：用户已筛选表，唯一Gene×Celltype×Treatment去重；不新增FDR/logFC过滤、不重跑DESeq2。',
               '该限制指PI的DEG计数分量；Mixscape仍按自身算法在表达对象内寻找特征基因，不读其他DEG文件。',
               'v3修复Mixscape工具槽名称读取，核验pvec与处理/条形码；成功的其他方法缓存经manifest核验后继续使用。',
               'v3.1显式设置 future::plan(multisession)，避免 RStudio 下 multicore fork 崩溃导致 "Execution halted"。',
               '相同阈值应在各比较间一致；用户筛选阈值可能不同于原文PI的FDR<0.05。',
               'PI = (minmax(UserDEG_count) + minmax(Augur AUC) + minmax(scDist) + minmax(Mixscape))/4。',
               '每个分量在Heat/Cold/PEG所有可用细胞间统一缩放；不分处理分别缩放。',
               '常数分量按原适配约定置0，表示无区分信息；缩放审计会标记。',
               '缺任一分量即无完整PI；没有DEG记录不等于已检验且0 DEG；Endodermis状态比较不纳入。',
               '4分量不等于论文4分辨率13分量；未把同一细胞类型DEG复制到四种分辨率。',
               '同一细胞类型和处理的前三项在细胞间共享，细胞内分布主要受Mixscape影响。',
               'PI是数据集内相对表达扰动，不是损伤概率、耐受性、差异显著性或因果证据。',
               '小提琴为描述性分布，白点为生物学样本均值；不把每个细胞当独立生物学重复。',
               '灰色表示参考背景或缺失分量；CK是参照，不人为指定CK PI=0。',
               '原05_cell_PI.csv和ensemble_PI_results.rds保留；新的结果都带PI_UserDEG_前缀。',
               '若旧Nature脚本仍读取原PI文件，它不会自动切换为本次用户DEG版；应明确指定新的results.rds/cell_PI.csv。',
               '论文：Robertson et al. New Phytologist (2026), DOI 10.1111/nph.71378, Methods 13 / Fig.4。'
  ),pi_path('README.txt'),useBytes=TRUE)
  if(!any(is.finite(cm$PI))) {
    status<-'No_complete_PI'
    failed<-bind(PI_JOBS);failed<-failed[!failed$OK,,drop=FALSE]
    if(nrow(failed))message(paste(paste(failed$Job,failed$Error,sep=': '),collapse='\n'))
    absent<-colnames(raw)[colSums(is.finite(raw))==0L]
    stop('没有完整四分量PI。完全缺失的分量：',if(length(absent))paste(absent,collapse=', ') else '无（各分量有效细胞未重叠）',
         '。详见 PI_UserDEG_jobs.csv、component_coverage.csv、missing_components.csv 和 Mixscape_*_tool_structure.txt。')
  }
  stage<-'Publication_plots';plot_user_pi(cm,sc$scaled,dat$umap,ss)
  status<-if(all(is.finite(cm$PI)))'Completed' else 'Completed_with_missing_PI'
  message('完成：',sum(is.finite(cm$PI)),'/',nrow(cm),' 个细胞有完整PI。图：',pi_path('Figure4.pdf'))
  invisible(list(cells=cm,samples=ss,output_dir=OUT_DIR))
}
if(!isTRUE(getOption('soy.userpi.skip_main',FALSE)))main()

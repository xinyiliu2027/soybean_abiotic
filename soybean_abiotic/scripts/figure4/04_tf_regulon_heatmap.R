# Source: 图4.docx, extracted text lines 2174-2687
# 读取01输出；GENIE3推断候选regulon
# See README.md and docs/KNOWN_LIMITATIONS.md before running.

##TF热图
# v1.2 修复GENIE3缺少doRNG：并行依赖预检、实际运行预检、缺依赖时单核继续；保留TF下载修复。
# 完整独立脚本：大豆 RNA 推断候选 regulon 热图（并非 ATAC/motif 活性）
# 首次安装：
# install.packages(c('SeuratObject','Matrix','pheatmap','digest'))
# 四核GENIE3所需依赖（缺失时本脚本自动单核运行）：
# install.packages(c('doRNG','doParallel','foreach','dplyr'))
# if(!requireNamespace('BiocManager',quietly=TRUE)) install.packages('BiocManager')
# BiocManager::install(c('GENIE3','UCell','BiocParallel','limma'))
OUT_DIR <- file.path(Sys.getenv("SOYBEAN_PROJECT_ROOT", unset = "data"), "CK2_added_v3_Results/Ensemble_PI_Heat_Cold_PEG_with_CK2")
INPUT_RDS <- file.path(OUT_DIR,'Endo3_v5_all_cells_modules_scored.rds')
MODULE_CSV <- file.path(OUT_DIR,'Endo3_v5_network_modules.csv')
TF_LIST_FILE <- file.path(OUT_DIR,'PlantTFDB_Gma_TF_list.txt')
# 可选核实后的映射表 original_id,gene_id；不要把不同基因组版本直接强行对应。
ID_MAP_FILE <- NULL
# 可选同一大豆基因版本的先验表 TF,target,source；设定后只保留该表支持的边。
# 默认没有使用 Cell 文章 ATAC 数据，也没有用根瘤网络代替根内皮层网络。
PRIOR_EDGES_FILE <- NULL
AUTO_DOWNLOAD_TF <- TRUE
PREFIX <- 'TFreg_v1_'
SEED <- 20261027L
N_CORES <- 4L # 请求核数；启动时根据依赖和微型运行检查决定实际核数
GENIE3_CORES <- 1L # main()中自动设置，不需要手动修改
N_TREES <- 500L
CELLS_PER_POOL <- 10L
MAX_POOLS_PER_SAMPLE_TYPE <- 12L
MAX_TARGET_GENES <- 5000L
MAX_TF_GENES <- 1500L
MIN_GENE_POOL_FRACTION <- 0.02 # 保留少数细胞类型/阶段特异基因
TARGET_BATCH <- 200L
TOP_TFS_PER_TARGET <- 5L
MAX_TARGETS_PER_TF <- 100L
MIN_TARGETS_PER_TF <- 10L
MIN_POOLED_R <- 0.15
MIN_SAMPLE_R <- 0.05
MIN_POSITIVE_SAMPLES <- 2L # 同时要求至少一个条件的两个重复均正相关
MIN_SAMPLE_TYPE_CELLS <- 20L
UCELL_MAX_RANK <- 1500L
SCORE_CHUNK <- 1500L
TOP_HEATMAP_TFS <- 60L
REUSE_CACHE <- TRUE
RUN_EXPLORATORY_TESTS <- TRUE
VERSION <- 'TFreg_v1.2_parallel_preflight'
design <- data.frame(Sample=c('CK1','CK2','highTemp1','highTemp2','lowTemp1','lowTemp2','PEG1','PEG2'),
                     Condition=rep(c('CK','Heat','Cold','PEG'),each=2),stringsAsFactors=FALSE)
conditions <- c('CK','Heat','Cold','PEG')
stages <- paste0('Endodermis_P',1:3)
type_order <- c('Root Cap','Epidermis','Cortex','Dividing Cell','Meristem','QC',stages,
                'Pericycle','Procambium','Phloem','Xylem')
out <- function(s) file.path(OUT_DIR,paste0(PREFIX,s))
tabout <- function(x,s) write.csv(x,out(s),row.names=FALSE,na='')
msg <- function(...) message(format(Sys.time(),'%H:%M:%S'),' ',paste0(...,collapse=''))
canon_gene <- function(x) {
  x<-as.character(x);y<-toupper(trimws(x))
  hit<-grepl('^GLYMA[._-]?[0-9]{2}G[0-9]{6}(\\.[0-9]+)*$',y)
  y[hit]<-sub('^GLYMA[._-]?([0-9]{2}G[0-9]{6}).*$','Glyma.\\1',y[hit])
  y[!hit]<-trimws(x[!hit]);y
}
# v1.1：不再依赖物种网页；优先本地，随后官方静态文件，支持内容识别和校验。
TF_DOWNLOAD_URLS <- c(
  'https://planttfdb.gao-lab.org/download/TF_list/Gma_TF_list.txt.gz',
  'http://planttfdb.gao-lab.org/download/TF_list/Gma_TF_list.txt.gz',
  'https://planttfdb.cbi.pku.edu.cn/download/TF_list/Gma_TF_list.txt.gz',
  'http://planttfdb.cbi.pku.edu.cn/download/TF_list/Gma_TF_list.txt.gz'
)
TF_DOWNLOAD_RETRIES <- 2L
TF_DOWNLOAD_TIMEOUT <- 45L
# source() 和 Rscript 两种运行方式都记录脚本所在目录，便于离线时把名单放在脚本旁边。
SCRIPT_DIR <- local({
  ff<-lapply(sys.frames(),function(z)z$ofile)
  ff<-Filter(function(z)is.character(z)&&length(z)==1L,ff)
  if(length(ff)) dirname(normalizePath(ff[[length(ff)]],mustWork=FALSE)) else {
    a<-grep('^--file=',commandArgs(FALSE),value=TRUE)
    if(length(a)) dirname(normalizePath(sub('^--file=','',a[1]),mustWork=FALSE)) else getwd()
  }
})
read_tf_text <- function(path) {
  if(!file.exists(path)||is.na(file.info(path)$size)||file.info(path)$size==0) stop('文件缺失或为空。')
  con<-file(path,'rb');magic<-tryCatch(readBin(con,'raw',n=2),finally=close(con))
  # 不根据文件扩展名猜压缩格式：改名为 .txt 的 gzip 也能读取。
  gz<-identical(as.integer(magic),c(31L,139L))
  con<-if(gz) gzfile(path,'rt') else file(path,'rt')
  lines<-tryCatch(withCallingHandlers(readLines(con,warn=FALSE),
                                      warning=function(w)stop(conditionMessage(w))),finally=close(con))
  if(!length(lines)) stop('没有文本内容。')
  if(any(grepl('<!doctype|<html|<body|<script',head(lines,30),ignore.case=TRUE)))
    stop('返回的是HTML网页/错误页，不是TF名单。')
  lines
}
parse_tf_text <- function(lines) {
  families<-strsplit('AP2 ARF ARR-B B3 BBR-BPC BES1 C2H2 C3H CAMTA CO-like CPP DBB Dof E2F/DP EIL ERF FAR1 G2-like GATA GeBP GRAS GRF HB-PHD HB-other HD-ZIP HRT-like HSF LBD LFY LSD M-type_MADS MIKC_MADS MYB MYB_related NAC NF-X1 NF-YA NF-YB NF-YC Nin-like RAV S1Fa-like SAP SBP SRS STAT TALE TCP Trihelix VOZ WOX WRKY Whirly YABBY ZF-HD bHLH bZIP',' ')[[1]]
  d<-do.call(rbind,lapply(lines,function(z) {
    ids<-regmatches(z,gregexpr('Glyma[._-]?[0-9]{2}G[0-9]{6}(\\.[0-9]+)*',z,ignore.case=TRUE))[[1]]
    if(!length(ids)) return(NULL)
    fields<-trimws(strsplit(gsub('"','',z),'[\t,; ]+')[[1]])
    fam<-families[tolower(families) %in% tolower(fields)]
    data.frame(TF=unique(canon_gene(ids)),Family=if(length(fam))paste(fam,collapse='/') else 'Unknown',stringsAsFactors=FALSE)
  }))
  if(is.null(d)||!nrow(d)) stop('找不到标准Glyma基因ID。')
  d<-aggregate(Family~TF,unique(d),function(z)paste(sort(unique(z)),collapse='/'))
  if(nrow(d)<100L) stop('少于100个不同的大豆TF基因；不能用重复行通过校验。')
  if(mean(d$Family=='Unknown')>0.5) warning('超过一半TF缺少家族标签；允许作为TF名单使用，但请核对来源。')
  d
}
fetch_tf_resource <- function(url,dest) {
  old<-options(timeout=TF_DOWNLOAD_TIMEOUT);on.exit(options(old),add=TRUE)
  method<-if(isTRUE(capabilities('libcurl'))) 'libcurl' else 'auto'
  status<-suppressWarnings(download.file(url,dest,method=method,mode='wb',quiet=TRUE))
  if(!identical(as.integer(status),0L)) stop('download.file 返回非零状态。')
  invisible(dest)
}
load_tf_list <- function() {
  if(!dir.exists(OUT_DIR)) dir.create(OUT_DIR,recursive=TRUE)
  cached<-file.path(OUT_DIR,'PlantTFDB_Gma_TF_list_validated.tsv')
  filenames<-c('PlantTFDB_Gma_TF_list.txt','PlantTFDB_Gma_TF_list.txt.gz',
               'Gma_TF_list.txt','Gma_TF_list.txt.gz','PlantTFDB_Gma_TF_list_validated.tsv')
  local_paths<-unique(c(TF_LIST_FILE,cached,paste0(TF_LIST_FILE,'.gz'),
                        unlist(lapply(unique(c(dirname(TF_LIST_FILE),OUT_DIR,SCRIPT_DIR,getwd())),function(p)file.path(p,filenames)))))
  audit<-data.frame(Source=character(),Attempt=integer(),Status=character(),Detail=character())
  record<-function(src,attempt,status,detail) {
    audit<<-rbind(audit,data.frame(Source=src,Attempt=attempt,Status=status,Detail=detail))
    tabout(audit,'TF_resource_attempts.csv')
  }
  accept<-function(d,src,rawpath) {
    checksum<-unname(tools::md5sum(rawpath))
    # 缓存为普通TSV；先写临时文件，校验后再替换缓存；不覆盖用户提供的原名单。
    candidate<-tempfile(tmpdir=OUT_DIR,fileext='.tsv')
    on.exit(unlink(candidate),add=TRUE)
    write.table(data.frame(Gene_ID=d$TF,Family=d$Family),candidate,sep='\t',quote=FALSE,row.names=FALSE)
    checked<-parse_tf_text(read_tf_text(candidate))
    if(!identical(d$TF,checked$TF)) stop('缓存写入后的TF校验失败。')
    if(!file.copy(candidate,cached,overwrite=TRUE)) stop('无法写入验证后的TF缓存。')
    tabout(d,'TF_catalog.csv')
    writeLines(c(paste('Source:',src),paste('Source_file_MD5:',checksum),
                 paste('Unique_TF_genes:',nrow(d)),paste('Validated_cache:',cached),
                 'PlantTFDB reference: Glycine max Wm82.a2.v1; ID syntax checks do not establish genome-version equivalence.',
                 paste('Loaded_at:',Sys.time())),out('TF_download_source.txt'))
    msg('TF名单已就绪：',nrow(d),'个不同基因；来源：',src)
    d
  }
  for(p in local_paths[file.exists(local_paths)]) {
    d<-tryCatch(parse_tf_text(read_tf_text(p)),error=function(e) {
      record(p,0L,'Invalid_local',conditionMessage(e));NULL
    })
    if(!is.null(d)) {record(p,0L,'Validated',paste(nrow(d),'unique genes'));return(accept(d,p,p))}
  }
  if(AUTO_DOWNLOAD_TF) {
    for(attempt in seq_len(TF_DOWNLOAD_RETRIES)) {
      for(u in TF_DOWNLOAD_URLS) {
        msg('尝试TF直接下载 [',attempt,'/',TF_DOWNLOAD_RETRIES,']：',u)
        tmp<-tempfile(fileext='.bin')
        d<-tryCatch({fetch_tf_resource(u,tmp);parse_tf_text(read_tf_text(tmp))},error=function(e) {
          record(u,attempt,'Failed',conditionMessage(e));NULL
        })
        if(!is.null(d)) {
          record(u,attempt,'Validated',paste(nrow(d),'unique genes'))
          result<-tryCatch(accept(d,u,tmp),finally=unlink(tmp));return(result)
        }
        unlink(tmp)
      }
      if(attempt<TF_DOWNLOAD_RETRIES) Sys.sleep(2)
    }
  }
  writeLines(c('TF resource unavailable: analysis has NOT started.',
               'Download the official Gma_TF_list.txt.gz on another network if necessary.',
               'Put the gz file in OUT_DIR or beside this script; decompression is optional.',
               TF_DOWNLOAD_URLS,paste('Output directory:',OUT_DIR),
               paste('Configured local path:',TF_LIST_FILE),
               paste('Attempt audit:',out('TF_resource_attempts.csv'))),out('TF_RESOURCE_REQUIRED.txt'))
  stop('未获取到有效大豆TF名单，网络分析尚未开始。\n',
       '请下载官方 Gma_TF_list.txt.gz，保持.gz格式即可，放到输出目录或脚本旁边。\n',
       '下载地址：',TF_DOWNLOAD_URLS[1],'\n',
       '逐个地址的失败原因见：',out('TF_resource_attempts.csv'),call.=FALSE)
}

load_data <- function() {
  if(!file.exists(INPUT_RDS)) stop('不存在：',INPUT_RDS,
                                   '\n需要含 Endodermis_P1/P2/P3 的 v5 对象；先完成 v5，或修改 INPUT_RDS。')
  obj<-readRDS(INPUT_RDS)
  if(!inherits(obj,'Seurat')) stop('输入必须为Seurat对象。')
  md<-obj@meta.data
  if(!all(c('Sample','Celltype') %in% names(md))) stop('需要metadata列Sample、Celltype：请使用v5输出。')
  if(is.null(rownames(md))||anyDuplicated(rownames(md))) stop('metadata细胞名缺失/重复。')
  md$Sample<-as.character(md$Sample);md$Celltype<-as.character(md$Celltype)
  md$Celltype[tolower(md$Celltype) %in% c('lateral root','procambium')]<-'Procambium'
  md<-md[md$Sample %in% design$Sample,,drop=FALSE]
  cond<-design$Condition[match(md$Sample,design$Sample)]
  if('Condition' %in% names(md)&&any(!is.na(md$Condition)&as.character(md$Condition)!=cond)) stop('样本与Condition冲突。')
  md$Condition<-cond;md$cell<-rownames(md)
  if(anyNA(md$Celltype)||any(!nzchar(md$Celltype))) stop('细胞类型缺失。')
  if(!all(stages %in% md$Celltype)) stop('缺少Endodermis_P1/P2/P3；不会把旧分支标签自动混为三阶段。')
  if(any(grepl('^Endodermis',md$Celltype)&!md$Celltype %in% stages)) stop('存在不属于三阶段的Endodermis标签。')
  eligible<-rep(TRUE,nrow(md))
  if('Analysis_eligible' %in% names(md)) {
    z<-as.logical(md$Analysis_eligible);eligible<-!is.na(z)&z
  }
  tabout(data.frame(cell=md$cell,Sample=md$Sample,Celltype=md$Celltype,Included=eligible),'cell_selection.csv')
  md<-md[eligible,,drop=FALSE]
  if(!setequal(md$Sample,design$Sample)) stop('未包含8个独立样本。')
  a<-obj@assays[['RNA']]
  if(is.null(a)) stop('缺少RNA assay。')
  if(inherits(a,'Assay5')) {
    a<-SeuratObject::JoinLayers(a,layers='counts',new='counts')
    counts<-SeuratObject::LayerData(a,layer='counts')
  } else if(utils::packageVersion('SeuratObject') >= '5.0.0') {
    counts<-SeuratObject::LayerData(a,layer='counts')
  } else counts<-SeuratObject::GetAssayData(a,slot='counts')
  if(!all(md$cell %in% colnames(counts))) stop('counts未覆盖metadata。')
  counts<-methods::as(counts[,md$cell,drop=FALSE],'dgCMatrix')
  if(any(!is.finite(counts@x))||any(counts@x<0)||any(abs(counts@x-round(counts@x))>1e-6)) stop('需要原始非负整数RNA counts。')
  original<-rownames(counts);ids<-canon_gene(original)
  if(!is.null(ID_MAP_FILE)) {
    mp<-read.csv(ID_MAP_FILE,stringsAsFactors=FALSE)
    stopifnot(all(c('original_id','gene_id') %in% names(mp)),!anyDuplicated(mp$original_id))
    ii<-match(original,mp$original_id);hit<-!is.na(ii);ids[hit]<-canon_gene(mp$gene_id[ii[hit]])
  }
  if(anyNA(ids)||any(!nzchar(ids))||anyDuplicated(ids)) stop('基因ID缺失或规范化后重复；请核对gene/transcript层级。')
  tabout(data.frame(original_id=original,gene_id=ids),'gene_id_mapping.csv');rownames(counts)<-ids
  if(any(Matrix::colSums(counts)==0)) stop('存在零counts细胞。')
  tabout(as.data.frame(table(Sample=md$Sample,Celltype=md$Celltype)),'sample_celltype_counts.csv')
  list(counts=counts,md=md)
}
make_pools <- function(counts,md) {
  set.seed(SEED);groups<-split(seq_len(nrow(md)),paste(md$Sample,md$Celltype,sep='|'))
  pools<-list();pm<-list();n<-0L
  for(ix in groups) {
    k<-min(MAX_POOLS_PER_SAMPLE_TYPE,length(ix)%/%CELLS_PER_POOL)
    if(k<2L) next
    ix<-sample(ix,k*CELLS_PER_POOL)
    for(j in seq_len(k)) {
      n<-n+1L;z<-ix[((j-1L)*CELLS_PER_POOL+1L):(j*CELLS_PER_POOL)]
      pools[[n]]<-z
      pm[[n]]<-data.frame(Pool=paste0('pool',n),Sample=md$Sample[z[1]],Condition=md$Condition[z[1]],
                          Celltype=md$Celltype[z[1]],Cells=length(z))
    }
  }
  if(n<40) stop('可用表达池少于40个。')
  pm<-do.call(rbind,pm)
  if(!setequal(pm$Sample,design$Sample)) stop('构建表达池后丢失样本。')
  tabout(pm,'pool_metadata.csv')
  tabout(data.frame(cell=md$cell[unlist(pools)],Pool=rep(pm$Pool,lengths(pools))),'pool_membership.csv')
  membership<-Matrix::sparseMatrix(i=unlist(pools),j=rep(seq_len(n),lengths(pools)),x=1,dims=c(ncol(counts),n))
  mat<-counts %*% membership;mat<-mat %*% Matrix::Diagonal(x=1e4/Matrix::colSums(mat))
  mat@x<-log1p(mat@x);colnames(mat)<-pm$Pool
  # 同一样本同一类型/阶段内随机无重复聚合。表达池用于网络学习，不作为生物学重复。
  list(expr=mat,meta=pm)
}
learn_regulons <- function(pools,tf,key) {
  e<-pools$expr;mu<-Matrix::rowMeans(e);sq<-e;sq@x<-sq@x^2
  vv<-Matrix::rowMeans(sq)-mu^2;rm(sq)
  keep<-Matrix::rowMeans(e>0)>=MIN_GENE_POOL_FRACTION&vv>1e-8
  ranked<-rownames(e)[order(vv,decreasing=TRUE)];genes<-ranked[ranked %in% rownames(e)[keep]]
  regulators<-head(genes[genes %in% tf$TF],MAX_TF_GENES)
  tabout(data.frame(TF=tf$TF,Detected=tf$TF %in% rownames(e),Selected=tf$TF %in% regulators),'TF_detection_audit.csv')
  if(length(regulators)<50) stop('有效TF不足50个，检查基因ID和基因组版本。')
  if(file.exists(MODULE_CSV)) {
    mod<-read.csv(MODULE_CSV,stringsAsFactors=FALSE)
    if(!all(c('gene_name','module') %in% names(mod))) stop('模块表需要gene_name,module列。')
    mg<-canon_gene(mod$gene_name[!is.na(mod$module)&!tolower(mod$module) %in% c('grey','gray','0')])
    targets<-head(genes[genes %in% mg],MAX_TARGET_GENES);src<-'assigned WGCNA genes ranked by RNA variance'
  } else {targets<-head(genes,MAX_TARGET_GENES);src<-'variable expressed RNA genes; module file absent'}
  if(length(targets)<100) stop('有效靶基因不足100个。')
  tabout(data.frame(gene=targets,Selection=src),'target_gene_universe.csv')
  x<-as.matrix(e[union(targets,regulators),,drop=FALSE]);key<-digest::digest(list(key,targets,regulators,x))
  msg('GENIE3：',length(regulators),' TF，',length(targets),'靶基因，',ncol(x),'表达池。')
  ed<-list();batches<-split(targets,ceiling(seq_along(targets)/TARGET_BATCH))
  for(b in seq_along(batches)) {
    f<-out(sprintf('cache_%s_GENIE3_%03d.rds',key,b))
    if(REUSE_CACHE&&file.exists(f)) z<-readRDS(f) else {
      set.seed(SEED+b)
      w<-GENIE3::GENIE3(x,regulators=regulators,targets=batches[[b]],treeMethod='RF',
                        K='sqrt',nTrees=N_TREES,nCores=GENIE3_CORES,verbose=TRUE)
      if(is.null(rownames(w))||is.null(colnames(w))) stop('GENIE3输出缺少基因名。')
      z<-do.call(rbind,lapply(batches[[b]],function(g) {
        ww<-w[regulators,g];names(ww)<-regulators
        ww<-head(sort(ww[names(ww)!=g&is.finite(ww)&ww>0],decreasing=TRUE),TOP_TFS_PER_TARGET)
        if(!length(ww)) return(NULL)
        data.frame(TF=names(ww),target=g,weight=as.numeric(ww))
      }));saveRDS(z,f)
    }
    ed[[b]]<-z;msg('完成网络批次 ',b,'/',length(batches))
  }
  ed<-do.call(rbind,ed)
  if(is.null(ed)||!nrow(ed)) stop('GENIE3没有有效边。')
  ed$r<-vapply(seq_len(nrow(ed)),function(i)cor(x[ed$TF[i],],x[ed$target[i],]),numeric(1))
  bysample<-split(seq_len(ncol(x)),pools$meta$Sample)
  rr<-vapply(bysample,function(ix)vapply(seq_len(nrow(ed)),function(i) {
    a<-x[ed$TF[i],ix];b<-x[ed$target[i],ix]
    if(length(ix)<10||sd(a)<1e-8||sd(b)<1e-8) return(NA_real_)
    cor(a,b)
  },numeric(1)),numeric(nrow(ed)))
  colnames(rr)<-paste0('r_',names(bysample));ed<-cbind(ed,rr)
  ed$positive_samples<-rowSums(rr>=MIN_SAMPLE_R,na.rm=TRUE)
  ed$replicate_pair_supported<-Reduce(`|`,lapply(conditions,function(co) {
    cn<-paste0('r_',design$Sample[design$Condition==co]);z<-rr[,cn,drop=FALSE]
    rowSums(is.finite(z)&z>=MIN_SAMPLE_R)==2L
  }))
  ed$source<-'RNA_GENIE3_positive_association'
  ed$keep<-is.finite(ed$r)&ed$r>=MIN_POOLED_R&ed$positive_samples>=MIN_POSITIVE_SAMPLES&ed$replicate_pair_supported
  if(!is.null(PRIOR_EDGES_FILE)) {
    pr<-read.csv(PRIOR_EDGES_FILE,stringsAsFactors=FALSE)
    if(!all(c('TF','target','source') %in% names(pr))) stop('先验表需要TF,target,source列。')
    if(anyNA(pr$source)||any(!nzchar(pr$source))) stop('先验来源不能为空。')
    pr$TF<-canon_gene(pr$TF);pr$target<-canon_gene(pr$target)
    m<-match(paste(ed$TF,ed$target),paste(pr$TF,pr$target));ed$keep<-ed$keep&!is.na(m)
    ed$source[!is.na(m)]<-paste0('RNA_plus_prior:',pr$source[m[!is.na(m)]])
  }
  tabout(ed,'all_candidate_edges_with_QC.csv');d<-ed[ed$keep,,drop=FALSE]
  if(!nrow(d)) stop('无通过筛选的边；检查审计表，不要为出图随意降低阈值。')
  d<-do.call(rbind,lapply(split(d,d$TF),function(z)head(z[order(-z$weight),],MAX_TARGETS_PER_TF)))
  n<-table(d$TF);d<-d[d$TF %in% names(n)[n>=MIN_TARGETS_PER_TF],,drop=FALSE]
  if(!nrow(d)) stop('无至少10个靶基因的候选regulon。')
  tabout(d,'retained_TF_target_edges.csv');reg<-lapply(split(d$target,d$TF),unique)
  saveRDS(reg,out('candidate_regulons.rds'));reg
}
score_cells <- function(counts,regulons) {
  # 用全部RNA基因作为UCell排名背景，不只保留TF/模块基因。
  safe<-sprintf('R%04d',seq_along(regulons));orig<-names(regulons);names(regulons)<-safe
  sc<-matrix(NA_real_,ncol(counts),length(regulons),dimnames=list(colnames(counts),orig))
  chunks<-split(seq_len(ncol(counts)),ceiling(seq_len(ncol(counts))/SCORE_CHUNK))
  for(k in seq_along(chunks)) {
    ix<-chunks[[k]]
    z<-as.matrix(UCell::ScoreSignatures_UCell(counts[,ix,drop=FALSE],features=regulons,
                                              maxRank=min(UCELL_MAX_RANK,nrow(counts)),BPPARAM=BiocParallel::SerialParam()))
    want<-paste0(safe,'_UCell')
    if(!all(colnames(counts)[ix] %in% rownames(z))||!all(want %in% colnames(z))) stop('UCell细胞/签名名称不匹配。')
    sc[ix,]<-z[colnames(counts)[ix],want,drop=FALSE]
    if(k%%10L==0L) msg('评分批次 ',k,'/',length(chunks))
  }
  if(any(!is.finite(sc))) stop('UCell存在非有限值。');sc
}
sample_average <- function(scores,md) {
  stopifnot(identical(rownames(scores),md$cell))
  groups<-split(seq_len(nrow(md)),paste(md$Sample,md$Celltype,sep='|'))
  groups<-groups[lengths(groups)>=MIN_SAMPLE_TYPE_CELLS]
  if(!length(groups)) stop('没有达到最小细胞数的样本×类型组合。')
  meta<-do.call(rbind,lapply(groups,function(ix)data.frame(Sample=md$Sample[ix[1]],
                                                           Condition=md$Condition[ix[1]],Celltype=md$Celltype[ix[1]],n_cells=length(ix))))
  mat<-do.call(rbind,lapply(groups,function(ix)colMeans(scores[ix,,drop=FALSE])))
  rownames(meta)<-rownames(mat)<-names(groups);list(mat=mat,meta=meta)
}
condition_average <- function(av) {
  keys<-unlist(lapply(conditions,function(co)paste(co,type_order,sep='|')))
  keys<-keys[keys %in% paste(av$meta$Condition,av$meta$Celltype,sep='|')]
  m<-matrix(NA_real_,length(keys),ncol(av$mat),dimnames=list(keys,colnames(av$mat)))
  for(k in keys) {
    ix<-which(paste(av$meta$Condition,av$meta$Celltype,sep='|')==k)
    # 两个生物学重复各占一半；缺重复则保留NA，不补零。
    if(length(ix)==2&&length(unique(av$meta$Sample[ix]))==2) m[k,]<-colMeans(av$mat[ix,,drop=FALSE])
  };m
}
z_columns <- function(m) {
  mu<-colMeans(m,na.rm=TRUE);s<-apply(m,2,sd,na.rm=TRUE);s[!is.finite(s)|s<1e-8]<-1
  sweep(sweep(m,2,mu,'-'),2,s,'/')
}
draw_heatmap <- function(m,tf,file,title) {
  tabout(data.frame(Group=rownames(m),m,check.names=FALSE),paste0(file,'_raw.csv'))
  good<-rowSums(is.finite(m))==ncol(m)
  if(any(!good)) tabout(data.frame(Group=rownames(m)[!good],Reason='missing biological replicate'),paste0(file,'_omitted_rows.csv'))
  m<-m[good,,drop=FALSE]
  if(nrow(m)<2||ncol(m)<2) {msg(file,'：行列不足，保留原始表。');return(FALSE)}
  z<-z_columns(m)
  tabout(data.frame(Group=rownames(z),z,check.names=FALSE),paste0(file,'_Z.csv'))
  ann<-data.frame(Family=tf$Family[match(colnames(z),tf$TF)],row.names=colnames(z));ann$Family[is.na(ann$Family)]<-'Unknown'
  width<-max(10,min(24,4+ncol(z)*0.18));height<-max(5,2+nrow(z)*0.25)
  z[z>3]<-3;z[z< -3]<- -3
  ph<-pheatmap::pheatmap(z,color=colorRampPalette(c('#2166AC','#67A9CF','#FFFFBF','#FDAE61','#B2182B'))(100),
                         breaks=seq(-3,3,length.out=101),scale='none',cluster_rows=FALSE,cluster_cols=TRUE,
                         annotation_col=ann,border_color=NA,fontsize=10,fontsize_col=7,angle_col=90,
                         main=paste0(title,'\nColumn Z-score; sample-balanced means'),silent=TRUE)
  pdf(out(paste0(file,'.pdf')),width=width,height=height,useDingbats=FALSE)
  grid::grid.newpage();grid::grid.draw(ph$gtable);dev.off()
  png(out(paste0(file,'.png')),width=width,height=height,units='in',res=220)
  grid::grid.newpage();grid::grid.draw(ph$gtable);dev.off();all(good)
}
test_treatments <- function(av) {
  ans<-list();status<-list()
  for(ct in unique(av$meta$Celltype)) {
    ix<-which(av$meta$Celltype==ct);ix<-ix[match(design$Sample,av$meta$Sample[ix])]
    status[[ct]]<-data.frame(Celltype=ct,Status=if(anyNA(ix))'Skipped: incomplete eight samples' else 'Tested')
    if(anyNA(ix)) next
    cond<-factor(av$meta$Condition[ix],levels=conditions);mm<-model.matrix(~0+cond);colnames(mm)<-conditions
    fit<-limma::lmFit(t(av$mat[ix,,drop=FALSE]),mm)
    contrast<-limma::makeContrasts(Heat-CK,Cold-CK,PEG-CK,levels=mm)
    fit<-limma::eBayes(limma::contrasts.fit(fit,contrast),trend=TRUE)
    for(j in 1:3) {
      d<-limma::topTable(fit,coef=j,number=Inf,sort.by='none')
      # limma的logFC列在这里是UCell分数差，不是RNA log2FC。
      ans[[length(ans)+1L]]<-data.frame(TF=rownames(d),Celltype=ct,Treatment=conditions[j+1L],
                                        Delta_UCell=d$logFC,P_value=d$P.Value,n_CK=2,n_Treatment=2,Interpretation='Exploratory_same_data_network_and_test')
    }
  }
  tabout(do.call(rbind,status),'comparison_status.csv')
  if(!length(ans)) return(NULL)
  d<-do.call(rbind,ans);d$FDR_global<-p.adjust(d$P_value,'BH')
  tabout(d,'exploratory_treatment_tests.csv');d
}
# v1.2：在读取大RDS和构网前检查GENIE3并行依赖，并做微型运行检查。
probe_genie3 <- function(cores) {
  # 保留调用者随机数状态；真实构网仍在各批次使用SEED+b。
  had_seed<-exists('.Random.seed',envir=.GlobalEnv,inherits=FALSE)
  if(had_seed) old_seed<-get('.Random.seed',envir=.GlobalEnv)
  on.exit({
    if(had_seed) assign('.Random.seed',old_seed,envir=.GlobalEnv)
    else if(exists('.Random.seed',envir=.GlobalEnv,inherits=FALSE)) rm('.Random.seed',envir=.GlobalEnv)
    if(cores>1L) {
      try(doParallel::stopImplicitCluster(),silent=TRUE)
      try(foreach::registerDoSEQ(),silent=TRUE)
    }
  },add=TRUE)
  set.seed(SEED)
  x<-matrix(rnorm(8L*20L),nrow=8L,dimnames=list(paste0('check_gene_',1:8),NULL))
  w<-GENIE3::GENIE3(x,regulators=rownames(x)[1:4],targets=rownames(x)[5:8],
                    treeMethod='RF',K='sqrt',nTrees=5L,nCores=cores,verbose=FALSE)
  if(!is.matrix(w)||any(!is.finite(w))||!all(rownames(x)[5:8] %in% colnames(w)))
    stop('GENIE3微型检查返回了无效权重矩阵。')
  invisible(TRUE)
}
configure_genie3 <- function(requested,
                             has_namespace=function(p)requireNamespace(p,quietly=TRUE),probe=probe_genie3) {
  if(length(requested)!=1L||!is.finite(requested)||requested<1||requested!=as.integer(requested))
    stop('N_CORES必须为正整数。')
  requested<-as.integer(requested)
  deps<-c('doRNG','doParallel','foreach','dplyr')
  present<-vapply(deps,function(p)isTRUE(tryCatch(has_namespace(p),error=function(e)FALSE)),logical(1))
  tabout(data.frame(Package=deps,Loadable=present),'GENIE3_parallel_dependencies.csv')
  effective<-requested;reason<-'Requested mode available'
  if(requested>1L&&any(!present)) {
    effective<-1L
    reason<-paste('Missing/unloadable parallel packages:',paste(deps[!present],collapse=', '))
    msg(reason,'；本次自动使用单核，保留500棵树等原分析参数。')
    msg("如需四核，在运行前安装：install.packages(c('doRNG','doParallel','foreach','dplyr'))")
  }
  error<-tryCatch({probe(effective);NULL},error=function(e)conditionMessage(e))
  if(!is.null(error)&&effective>1L) {
    reason<-paste('Parallel preflight failed:',error)
    effective<-1L;msg(reason,'；检查单核模式。')
    error<-tryCatch({probe(1L);NULL},error=function(e)conditionMessage(e))
  }
  state<-data.frame(Requested_cores=requested,Effective_cores=effective,
                    Preflight=if(is.null(error))'Passed' else 'Failed',Reason=reason,
                    Error=if(is.null(error))'' else error)
  tabout(state,'GENIE3_runtime_mode.csv')
  if(!is.null(error)) stop('GENIE3单核检查也失败：',error,
                           '\n请检查GENIE3安装；尚未读取分析RDS。',call.=FALSE)
  msg('GENIE3运行检查通过；请求 ',requested,' 核，实际使用 ',effective,' 核。')
  effective
}

main <- function() {
  options(timeout=max(300,getOption('timeout')))
  packages<-c('SeuratObject','Matrix','pheatmap','digest','GENIE3','UCell','BiocParallel')
  if(RUN_EXPLORATORY_TESTS) packages<-c(packages,'limma')
  missing<-packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
  if(length(missing)) stop('缺少包：',paste(missing,collapse=', '),'；安装命令见脚本开头。')
  if(!dir.exists(OUT_DIR)) dir.create(OUT_DIR,recursive=TRUE)
  # 重跑只更新 TFreg_v1_ 文件及下载的TF名单，不写回输入RDS/旧WGCNA结果。
  for(status_file in c('SUCCESS.txt','PARTIAL.txt','TF_RESOURCE_REQUIRED.txt')) if(file.exists(out(status_file))) unlink(out(status_file))
  GENIE3_CORES <<- configure_genie3(N_CORES)
  tf<-load_tf_list();dat<-load_data();type_order<<-unique(c(type_order,sort(unique(dat$md$Celltype))))
  msg('读入 ',ncol(dat$counts),' 个合格细胞。');pools<-make_pools(dat$counts,dat$md)
  key<-digest::digest(list(VERSION,file.info(INPUT_RDS)[c('size','mtime')],tf,
                           if(file.exists(MODULE_CSV))tools::md5sum(MODULE_CSV) else NULL,
                           if(!is.null(PRIOR_EDGES_FILE))tools::md5sum(PRIOR_EDGES_FILE) else NULL,
                           if(!is.null(ID_MAP_FILE))tools::md5sum(ID_MAP_FILE) else NULL,
                           SEED,N_TREES,GENIE3_CORES,TOP_TFS_PER_TARGET,MAX_TARGET_GENES,MAX_TF_GENES,TARGET_BATCH,
                           MIN_GENE_POOL_FRACTION,MIN_POOLED_R,MIN_SAMPLE_R,MIN_POSITIVE_SAMPLES,packageVersion('GENIE3')))
  regulons<-learn_regulons(pools,tf,key)
  skey<-digest::digest(list(key,regulons,tools::md5sum(INPUT_RDS),UCELL_MAX_RANK,packageVersion('UCell')))
  sf<-out(paste0('cache_',skey,'_UCell.rds'))
  if(REUSE_CACHE&&file.exists(sf)) scores<-readRDS(sf) else {scores<-score_cells(dat$counts,regulons);saveRDS(scores,sf)}
  av<-sample_average(scores,dat$md)
  tabout(cbind(av$meta,as.data.frame(av$mat,check.names=FALSE)),'sample_regulon_scores.csv')
  means<-condition_average(av);ck<-means[grepl('^CK[|]',rownames(means)),,drop=FALSE];rownames(ck)<-sub('^CK[|]','',rownames(ck))
  v<-apply(means,2,var,na.rm=TRUE);v[!is.finite(v)]<- -Inf
  selected<-names(head(sort(v[v>1e-12],decreasing=TRUE),TOP_HEATMAP_TFS))
  if(length(selected)<2) stop('有变异的regulon不足两个。')
  tabout(data.frame(TF=names(regulons),Targets=lengths(regulons),Family=tf$Family[match(names(regulons),tf$TF)],
                    Displayed=names(regulons) %in% selected),'regulon_summary.csv')
  tabout(data.frame(Group=rownames(means),means,check.names=FALSE),'all_regulons_condition_means.csv')
  saveRDS(list(scores=scores,metadata=dat$md,regulons=regulons,VERSION=VERSION),out('scores_and_regulons.rds'))
  plots<-list()
  plots$CK<-draw_heatmap(ck[,selected,drop=FALSE],tf,'Figure01_CK_celltype_regulon','RNA-inferred candidate regulons | CK')
  plots$All<-draw_heatmap(means[,selected,drop=FALSE],tf,'Figure02_all_conditions_regulon','RNA-inferred candidate regulons | all conditions')
  ek<-unlist(lapply(conditions,function(co)paste(co,stages,sep='|')))
  em<-matrix(NA_real_,length(ek),ncol(means),dimnames=list(ek,colnames(means)))
  present<-intersect(ek,rownames(means));em[present,]<-means[present,,drop=FALSE]
  ev<-apply(em,2,var,na.rm=TRUE);ev[!is.finite(ev)]<- -Inf
  es<-names(head(sort(ev[ev>1e-12],decreasing=TRUE),TOP_HEATMAP_TFS))
  tabout(data.frame(TF=es,Selection=rep('Endodermis variance, not significance',length(es))),'Endodermis_displayed_TFs.csv')
  plots$Endo<-draw_heatmap(em[,es,drop=FALSE],tf,'Figure03_Endodermis_P1P2P3_regulon','Endodermis candidate regulons | P1 to P3')
  eix<-which(av$meta$Celltype %in% stages)
  eix<-eix[order(match(av$meta$Sample[eix],design$Sample),match(av$meta$Celltype[eix],stages))]
  plots$Replicates<-draw_heatmap(av$mat[eix,es,drop=FALSE],tf,'Figure04_Endodermis_individual_samples','Endodermis | individual biological samples')
  # TF自身表达的对照图；归一化分母仍是全部RNA counts。
  ex<-dat$counts[selected,,drop=FALSE] %*% Matrix::Diagonal(x=1e4/Matrix::colSums(dat$counts))
  ex@x<-log1p(ex@x);colnames(ex)<-dat$md$cell
  ea<-sample_average(t(as.matrix(ex)),dat$md);exprmean<-condition_average(ea)
  plots$Expression<-draw_heatmap(exprmean,tf,'Figure05_TF_expression_control','TF expression | mean log1p(CP10K), not regulon activity')
  if(RUN_EXPLORATORY_TESTS) test_treatments(av)
  tabout(data.frame(Figure=names(plots),Generated=unlist(plots)),'figure_status.csv')
  capture.output(sessionInfo(),file=out('sessionInfo.txt'))
  info<-c('RNA-inferred candidate regulons; NOT experimentally verified direct binding.',
          'Positive coexpression is not proof of transcriptional activation.',
          'Pools are not biological replicates. Treatment tests use the eight original samples.',
          'Learning and testing use the same dataset; treatment statistics are exploratory.',
          'P1/P2/P3 are inherited from the input. No reclustering or new trajectory was run.',
          paste('Prior table:',if(is.null(PRIOR_EDGES_FILE))'None (RNA-only inference)' else PRIOR_EDGES_FILE),
          paste('GENIE3 actual cores:',GENIE3_CORES),
          paste('Output:',OUT_DIR),paste('Input:',INPUT_RDS))
  if(all(unlist(plots))) {writeLines(info,out('SUCCESS.txt'));msg('全部完成。输出：',OUT_DIR,'；前缀：',PREFIX)}
  else {writeLines(info,out('PARTIAL.txt'));msg('评分完成但部分热图未生成，见 figure_status.csv 和 omitted_rows.csv。')}
}
if(sys.nframe()==0L) main()

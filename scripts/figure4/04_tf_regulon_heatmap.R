# Source: Figure 4 document, extracted text lines 2174-2687
# Reads script 01 outputs and infers candidate regulons using GENIE3.
# See README.md and docs/KNOWN_LIMITATIONS.md before running.

# TF heatmap.
# v1.2: precheck GENIE3 parallel dependencies and a small runtime example; fall back to one core if needed, retaining TF download fixes.
# Standalone heatmap of soybean RNA-inferred candidate regulons, not ATAC/motif activity.
# Initial installation commands.
# install.packages(c('SeuratObject','Matrix','pheatmap','digest'))
# Dependencies for four-core GENIE3; use one core automatically if unavailable.
# install.packages(c('doRNG','doParallel','foreach','dplyr'))
# if(!requireNamespace('BiocManager',quietly=TRUE)) install.packages('BiocManager')
# BiocManager::install(c('GENIE3','UCell','BiocParallel','limma'))
OUT_DIR <- file.path(Sys.getenv("SOYBEAN_PROJECT_ROOT", unset = "data"), "CK2_added_v3_Results/Ensemble_PI_Heat_Cold_PEG_with_CK2")
INPUT_RDS <- file.path(OUT_DIR,'Endo3_v5_all_cells_modules_scored.rds')
MODULE_CSV <- file.path(OUT_DIR,'Endo3_v5_network_modules.csv')
TF_LIST_FILE <- file.path(OUT_DIR,'PlantTFDB_Gma_TF_list.txt')
# Optional verified original_id,gene_id mapping; do not force direct matches across genome versions.
ID_MAP_FILE <- NULL
# Optional TF,target,source prior-edge table using the same soybean gene version; when provided, retain only supported edges.
# By default, neither ATAC data from the cited Cell paper nor a nodule network is used as a root-Endodermis network.
PRIOR_EDGES_FILE <- NULL
AUTO_DOWNLOAD_TF <- TRUE
PREFIX <- 'TFreg_v1_'
SEED <- 20261027L
N_CORES <- 4L # Requested core count; dependency and small runtime checks determine the effective count.
GENIE3_CORES <- 1L # Set automatically in main(); no manual change is needed.
N_TREES <- 500L
CELLS_PER_POOL <- 10L
MAX_POOLS_PER_SAMPLE_TYPE <- 12L
MAX_TARGET_GENES <- 5000L
MAX_TF_GENES <- 1500L
MIN_GENE_POOL_FRACTION <- 0.02 # Retain genes specific to rare cell types or stages.
TARGET_BATCH <- 200L
TOP_TFS_PER_TARGET <- 5L
MAX_TARGETS_PER_TF <- 100L
MIN_TARGETS_PER_TF <- 10L
MIN_POOLED_R <- 0.15
MIN_SAMPLE_R <- 0.05
MIN_POSITIVE_SAMPLES <- 2L # Also require positive correlations in both replicates of at least one condition.
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
# v1.1: use local TF files first, then official static files, with content detection and validation rather than a species webpage.
TF_DOWNLOAD_URLS <- c(
  'https://planttfdb.gao-lab.org/download/TF_list/Gma_TF_list.txt.gz',
  'http://planttfdb.gao-lab.org/download/TF_list/Gma_TF_list.txt.gz',
  'https://planttfdb.cbi.pku.edu.cn/download/TF_list/Gma_TF_list.txt.gz',
  'http://planttfdb.cbi.pku.edu.cn/download/TF_list/Gma_TF_list.txt.gz'
)
TF_DOWNLOAD_RETRIES <- 2L
TF_DOWNLOAD_TIMEOUT <- 45L
# Record the script directory under both source() and Rscript so an offline TF list can be placed beside the script.
SCRIPT_DIR <- local({
  ff<-lapply(sys.frames(),function(z)z$ofile)
  ff<-Filter(function(z)is.character(z)&&length(z)==1L,ff)
  if(length(ff)) dirname(normalizePath(ff[[length(ff)]],mustWork=FALSE)) else {
    a<-grep('^--file=',commandArgs(FALSE),value=TRUE)
    if(length(a)) dirname(normalizePath(sub('^--file=','',a[1]),mustWork=FALSE)) else getwd()
  }
})
read_tf_text <- function(path) {
  if(!file.exists(path)||is.na(file.info(path)$size)||file.info(path)$size==0) stop("File is missing or empty.")
  con<-file(path,'rb');magic<-tryCatch(readBin(con,'raw',n=2),finally=close(con))
  # Detect compression from content rather than extension; gzip files renamed .txt remain readable.
  gz<-identical(as.integer(magic),c(31L,139L))
  con<-if(gz) gzfile(path,'rt') else file(path,'rt')
  lines<-tryCatch(withCallingHandlers(readLines(con,warn=FALSE),
                                      warning=function(w)stop(conditionMessage(w))),finally=close(con))
  if(!length(lines)) stop("No text content.")
  if(any(grepl('<!doctype|<html|<body|<script',head(lines,30),ignore.case=TRUE)))
    stop("Received an HTML/error page rather than a TF list.")
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
  if(is.null(d)||!nrow(d)) stop("No standard Glyma gene IDs found.")
  d<-aggregate(Family~TF,unique(d),function(z)paste(sort(unique(z)),collapse='/'))
  if(nrow(d)<100L) stop("Fewer than 100 distinct soybean TF genes; duplicated rows cannot satisfy validation.")
  if(mean(d$Family=='Unknown')>0.5) warning("More than half the TFs lack family labels; the list can be used, but verify its source.")
  d
}
fetch_tf_resource <- function(url,dest) {
  old<-options(timeout=TF_DOWNLOAD_TIMEOUT);on.exit(options(old),add=TRUE)
  method<-if(isTRUE(capabilities('libcurl'))) 'libcurl' else 'auto'
  status<-suppressWarnings(download.file(url,dest,method=method,mode='wb',quiet=TRUE))
  if(!identical(as.integer(status),0L)) stop("download.file returned a nonzero status.")
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
    # Cache as plain TSV using a validated temporary file before replacement; do not overwrite a user-supplied original list.
    candidate<-tempfile(tmpdir=OUT_DIR,fileext='.tsv')
    on.exit(unlink(candidate),add=TRUE)
    write.table(data.frame(Gene_ID=d$TF,Family=d$Family),candidate,sep='\t',quote=FALSE,row.names=FALSE)
    checked<-parse_tf_text(read_tf_text(candidate))
    if(!identical(d$TF,checked$TF)) stop("TF validation failed after writing the cache.")
    if(!file.copy(candidate,cached,overwrite=TRUE)) stop("Cannot write the validated TF cache.")
    tabout(d,'TF_catalog.csv')
    writeLines(c(paste('Source:',src),paste('Source_file_MD5:',checksum),
                 paste('Unique_TF_genes:',nrow(d)),paste('Validated_cache:',cached),
                 'PlantTFDB reference: Glycine max Wm82.a2.v1; ID syntax checks do not establish genome-version equivalence.',
                 paste('Loaded_at:',Sys.time())),out('TF_download_source.txt'))
    msg("TF list ready: ",nrow(d)," distinct genes; source: ",src)
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
        msg("Trying direct TF download [",attempt,'/',TF_DOWNLOAD_RETRIES,']：',u)
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
  stop("No valid soybean TF list obtained; network analysis has not started.\n",
       "Download the official Gma_TF_list.txt.gz and place it, still compressed, in the output or script directory.\n",
       "Download URL: ",TF_DOWNLOAD_URLS[1],'\n',
       "Per-URL failure details: ",out('TF_resource_attempts.csv'),call.=FALSE)
}

load_data <- function() {
  if(!file.exists(INPUT_RDS)) stop("Does not exist: ",INPUT_RDS,
                                   "\nA v5 object with Endodermis_P1/P2/P3 is required; complete v5 first or change INPUT_RDS.")
  obj<-readRDS(INPUT_RDS)
  if(!inherits(obj,'Seurat')) stop("Input must be a Seurat object.")
  md<-obj@meta.data
  if(!all(c('Sample','Celltype') %in% names(md))) stop("Metadata requires Sample and Celltype; use the v5 output.")
  if(is.null(rownames(md))||anyDuplicated(rownames(md))) stop("Metadata cell names are missing/duplicated.")
  md$Sample<-as.character(md$Sample);md$Celltype<-as.character(md$Celltype)
  md$Celltype[tolower(md$Celltype) %in% c('lateral root','procambium')]<-'Procambium'
  md<-md[md$Sample %in% design$Sample,,drop=FALSE]
  cond<-design$Condition[match(md$Sample,design$Sample)]
  if('Condition' %in% names(md)&&any(!is.na(md$Condition)&as.character(md$Condition)!=cond)) stop("Sample identity conflicts with Condition.")
  md$Condition<-cond;md$cell<-rownames(md)
  if(anyNA(md$Celltype)||any(!nzchar(md$Celltype))) stop("Missing cell types.")
  if(!all(stages %in% md$Celltype)) stop("Endodermis_P1/P2/P3 missing; legacy branch labels will not be treated as three stages automatically.")
  if(any(grepl('^Endodermis',md$Celltype)&!md$Celltype %in% stages)) stop("Endodermis labels outside the three stages detected.")
  eligible<-rep(TRUE,nrow(md))
  if('Analysis_eligible' %in% names(md)) {
    z<-as.logical(md$Analysis_eligible);eligible<-!is.na(z)&z
  }
  tabout(data.frame(cell=md$cell,Sample=md$Sample,Celltype=md$Celltype,Included=eligible),'cell_selection.csv')
  md<-md[eligible,,drop=FALSE]
  if(!setequal(md$Sample,design$Sample)) stop("Not all eight independent samples are present.")
  a<-obj@assays[['RNA']]
  if(is.null(a)) stop("RNA assay missing.")
  if(inherits(a,'Assay5')) {
    a<-SeuratObject::JoinLayers(a,layers='counts',new='counts')
    counts<-SeuratObject::LayerData(a,layer='counts')
  } else if(utils::packageVersion('SeuratObject') >= '5.0.0') {
    counts<-SeuratObject::LayerData(a,layer='counts')
  } else counts<-SeuratObject::GetAssayData(a,slot='counts')
  if(!all(md$cell %in% colnames(counts))) stop("Counts do not cover metadata cells.")
  counts<-methods::as(counts[,md$cell,drop=FALSE],'dgCMatrix')
  if(any(!is.finite(counts@x))||any(counts@x<0)||any(abs(counts@x-round(counts@x))>1e-6)) stop("Raw nonnegative integer RNA counts are required.")
  original<-rownames(counts);ids<-canon_gene(original)
  if(!is.null(ID_MAP_FILE)) {
    mp<-read.csv(ID_MAP_FILE,stringsAsFactors=FALSE)
    stopifnot(all(c('original_id','gene_id') %in% names(mp)),!anyDuplicated(mp$original_id))
    ii<-match(original,mp$original_id);hit<-!is.na(ii);ids[hit]<-canon_gene(mp$gene_id[ii[hit]])
  }
  if(anyNA(ids)||any(!nzchar(ids))||anyDuplicated(ids)) stop("Gene IDs are missing or duplicated after normalization; check gene versus transcript identifiers.")
  tabout(data.frame(original_id=original,gene_id=ids),'gene_id_mapping.csv');rownames(counts)<-ids
  if(any(Matrix::colSums(counts)==0)) stop("Zero-count cells detected.")
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
  if(n<40) stop("Fewer than 40 usable expression pools.")
  pm<-do.call(rbind,pm)
  if(!setequal(pm$Sample,design$Sample)) stop("Samples were lost during expression-pool construction.")
  tabout(pm,'pool_metadata.csv')
  tabout(data.frame(cell=md$cell[unlist(pools)],Pool=rep(pm$Pool,lengths(pools))),'pool_membership.csv')
  membership<-Matrix::sparseMatrix(i=unlist(pools),j=rep(seq_len(n),lengths(pools)),x=1,dims=c(ncol(counts),n))
  mat<-counts %*% membership;mat<-mat %*% Matrix::Diagonal(x=1e4/Matrix::colSums(mat))
  mat@x<-log1p(mat@x);colnames(mat)<-pm$Pool
  # Aggregate disjoint random cells within each sample and cell type/stage; expression pools support network learning, not biological replication.
  list(expr=mat,meta=pm)
}
learn_regulons <- function(pools,tf,key) {
  e<-pools$expr;mu<-Matrix::rowMeans(e);sq<-e;sq@x<-sq@x^2
  vv<-Matrix::rowMeans(sq)-mu^2;rm(sq)
  keep<-Matrix::rowMeans(e>0)>=MIN_GENE_POOL_FRACTION&vv>1e-8
  ranked<-rownames(e)[order(vv,decreasing=TRUE)];genes<-ranked[ranked %in% rownames(e)[keep]]
  regulators<-head(genes[genes %in% tf$TF],MAX_TF_GENES)
  tabout(data.frame(TF=tf$TF,Detected=tf$TF %in% rownames(e),Selected=tf$TF %in% regulators),'TF_detection_audit.csv')
  if(length(regulators)<50) stop("Fewer than 50 valid TFs; check gene IDs and genome version.")
  if(file.exists(MODULE_CSV)) {
    mod<-read.csv(MODULE_CSV,stringsAsFactors=FALSE)
    if(!all(c('gene_name','module') %in% names(mod))) stop("Module table requires gene_name and module columns.")
    mg<-canon_gene(mod$gene_name[!is.na(mod$module)&!tolower(mod$module) %in% c('grey','gray','0')])
    targets<-head(genes[genes %in% mg],MAX_TARGET_GENES);src<-'assigned WGCNA genes ranked by RNA variance'
  } else {targets<-head(genes,MAX_TARGET_GENES);src<-'variable expressed RNA genes; module file absent'}
  if(length(targets)<100) stop("Fewer than 100 valid target genes.")
  tabout(data.frame(gene=targets,Selection=src),'target_gene_universe.csv')
  x<-as.matrix(e[union(targets,regulators),,drop=FALSE]);key<-digest::digest(list(key,targets,regulators,x))
  msg('GENIE3：',length(regulators),' TF，',length(targets),"target genes, ",ncol(x),"expression pools.")
  ed<-list();batches<-split(targets,ceiling(seq_along(targets)/TARGET_BATCH))
  for(b in seq_along(batches)) {
    f<-out(sprintf('cache_%s_GENIE3_%03d.rds',key,b))
    if(REUSE_CACHE&&file.exists(f)) z<-readRDS(f) else {
      set.seed(SEED+b)
      w<-GENIE3::GENIE3(x,regulators=regulators,targets=batches[[b]],treeMethod='RF',
                        K='sqrt',nTrees=N_TREES,nCores=GENIE3_CORES,verbose=TRUE)
      if(is.null(rownames(w))||is.null(colnames(w))) stop("GENIE3 output lacks gene names.")
      z<-do.call(rbind,lapply(batches[[b]],function(g) {
        ww<-w[regulators,g];names(ww)<-regulators
        ww<-head(sort(ww[names(ww)!=g&is.finite(ww)&ww>0],decreasing=TRUE),TOP_TFS_PER_TARGET)
        if(!length(ww)) return(NULL)
        data.frame(TF=names(ww),target=g,weight=as.numeric(ww))
      }));saveRDS(z,f)
    }
    ed[[b]]<-z;msg("Completed network batch ",b,'/',length(batches))
  }
  ed<-do.call(rbind,ed)
  if(is.null(ed)||!nrow(ed)) stop("GENIE3 produced no valid edges.")
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
    if(!all(c('TF','target','source') %in% names(pr))) stop("Prior table requires TF, target and source columns.")
    if(anyNA(pr$source)||any(!nzchar(pr$source))) stop("Prior sources cannot be empty.")
    pr$TF<-canon_gene(pr$TF);pr$target<-canon_gene(pr$target)
    m<-match(paste(ed$TF,ed$target),paste(pr$TF,pr$target));ed$keep<-ed$keep&!is.na(m)
    ed$source[!is.na(m)]<-paste0('RNA_plus_prior:',pr$source[m[!is.na(m)]])
  }
  tabout(ed,'all_candidate_edges_with_QC.csv');d<-ed[ed$keep,,drop=FALSE]
  if(!nrow(d)) stop("No edges passed filtering; inspect audit tables rather than lowering thresholds just to obtain a plot.")
  d<-do.call(rbind,lapply(split(d,d$TF),function(z)head(z[order(-z$weight),],MAX_TARGETS_PER_TF)))
  n<-table(d$TF);d<-d[d$TF %in% names(n)[n>=MIN_TARGETS_PER_TF],,drop=FALSE]
  if(!nrow(d)) stop("No candidate regulon contains at least 10 target genes.")
  tabout(d,'retained_TF_target_edges.csv');reg<-lapply(split(d$target,d$TF),unique)
  saveRDS(reg,out('candidate_regulons.rds'));reg
}
score_cells <- function(counts,regulons) {
  # Use all RNA genes as the UCell ranking background, not only TF or module genes.
  safe<-sprintf('R%04d',seq_along(regulons));orig<-names(regulons);names(regulons)<-safe
  sc<-matrix(NA_real_,ncol(counts),length(regulons),dimnames=list(colnames(counts),orig))
  chunks<-split(seq_len(ncol(counts)),ceiling(seq_len(ncol(counts))/SCORE_CHUNK))
  for(k in seq_along(chunks)) {
    ix<-chunks[[k]]
    z<-as.matrix(UCell::ScoreSignatures_UCell(counts[,ix,drop=FALSE],features=regulons,
                                              maxRank=min(UCELL_MAX_RANK,nrow(counts)),BPPARAM=BiocParallel::SerialParam()))
    want<-paste0(safe,'_UCell')
    if(!all(colnames(counts)[ix] %in% rownames(z))||!all(want %in% colnames(z))) stop("UCell cell/signature names do not match.")
    sc[ix,]<-z[colnames(counts)[ix],want,drop=FALSE]
    if(k%%10L==0L) msg("Scoring batch ",k,'/',length(chunks))
  }
  if(any(!is.finite(sc))) stop("UCell contains non-finite values.");sc
}
sample_average <- function(scores,md) {
  stopifnot(identical(rownames(scores),md$cell))
  groups<-split(seq_len(nrow(md)),paste(md$Sample,md$Celltype,sep='|'))
  groups<-groups[lengths(groups)>=MIN_SAMPLE_TYPE_CELLS]
  if(!length(groups)) stop("No sample-by-type combination meets the minimum cell count.")
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
    # Give equal weight to the two biological replicates; retain NA when a replicate is missing rather than filling zero.
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
  if(nrow(m)<2||ncol(m)<2) {msg(file,": insufficient rows/columns; retaining the raw table.");return(FALSE)}
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
      # Here limma's logFC column is a UCell score difference, not RNA log2FC.
      ans[[length(ans)+1L]]<-data.frame(TF=rownames(d),Celltype=ct,Treatment=conditions[j+1L],
                                        Delta_UCell=d$logFC,P_value=d$P.Value,n_CK=2,n_Treatment=2,Interpretation='Exploratory_same_data_network_and_test')
    }
  }
  tabout(do.call(rbind,status),'comparison_status.csv')
  if(!length(ans)) return(NULL)
  d<-do.call(rbind,ans);d$FDR_global<-p.adjust(d$P_value,'BH')
  tabout(d,'exploratory_treatment_tests.csv');d
}
# v1.2: check GENIE3 parallel dependencies and a small runtime example before loading the large RDS or building the network.
probe_genie3 <- function(cores) {
  # Preserve the caller's random state; actual network batches still use SEED + b.
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
    stop("GENIE3 runtime precheck returned an invalid weight matrix.")
  invisible(TRUE)
}
configure_genie3 <- function(requested,
                             has_namespace=function(p)requireNamespace(p,quietly=TRUE),probe=probe_genie3) {
  if(length(requested)!=1L||!is.finite(requested)||requested<1||requested!=as.integer(requested))
    stop("N_CORES must be a positive integer.")
  requested<-as.integer(requested)
  deps<-c('doRNG','doParallel','foreach','dplyr')
  present<-vapply(deps,function(p)isTRUE(tryCatch(has_namespace(p),error=function(e)FALSE)),logical(1))
  tabout(data.frame(Package=deps,Loadable=present),'GENIE3_parallel_dependencies.csv')
  effective<-requested;reason<-'Requested mode available'
  if(requested>1L&&any(!present)) {
    effective<-1L
    reason<-paste('Missing/unloadable parallel packages:',paste(deps[!present],collapse=', '))
    msg(reason,"; automatically using one core while retaining 500 trees and other analysis parameters.")
    msg("For four cores, first run: install.packages(c('doRNG','doParallel','foreach','dplyr'))")
  }
  error<-tryCatch({probe(effective);NULL},error=function(e)conditionMessage(e))
  if(!is.null(error)&&effective>1L) {
    reason<-paste('Parallel preflight failed:',error)
    effective<-1L;msg(reason,"; checking single-core mode.")
    error<-tryCatch({probe(1L);NULL},error=function(e)conditionMessage(e))
  }
  state<-data.frame(Requested_cores=requested,Effective_cores=effective,
                    Preflight=if(is.null(error))'Passed' else 'Failed',Reason=reason,
                    Error=if(is.null(error))'' else error)
  tabout(state,'GENIE3_runtime_mode.csv')
  if(!is.null(error)) stop("GENIE3 single-core precheck also failed: ",error,
                           "\nCheck the GENIE3 installation; the analysis RDS has not been read.",call.=FALSE)
  msg("GENIE3 precheck passed; requested ",requested," cores, using ",effective," cores.")
  effective
}

main <- function() {
  options(timeout=max(300,getOption('timeout')))
  packages<-c('SeuratObject','Matrix','pheatmap','digest','GENIE3','UCell','BiocParallel')
  if(RUN_EXPLORATORY_TESTS) packages<-c(packages,'limma')
  missing<-packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
  if(length(missing)) stop("Missing packages: ",paste(missing,collapse=', '),"; installation commands are at the top of the script.")
  if(!dir.exists(OUT_DIR)) dir.create(OUT_DIR,recursive=TRUE)
  # Reruns update only TFreg_v1_ outputs and the downloaded TF list, preserving the input RDS and existing WGCNA results.
  for(status_file in c('SUCCESS.txt','PARTIAL.txt','TF_RESOURCE_REQUIRED.txt')) if(file.exists(out(status_file))) unlink(out(status_file))
  GENIE3_CORES <<- configure_genie3(N_CORES)
  tf<-load_tf_list();dat<-load_data();type_order<<-unique(c(type_order,sort(unique(dat$md$Celltype))))
  msg("Loaded ",ncol(dat$counts)," eligible cells.");pools<-make_pools(dat$counts,dat$md)
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
  if(length(selected)<2) stop("Fewer than two variable regulons.")
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
  # Control plot of TF expression; the normalization denominator still uses all RNA counts.
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
  if(all(unlist(plots))) {writeLines(info,out('SUCCESS.txt'));msg("Completed. Output: ",OUT_DIR,"; prefix: ",PREFIX)}
  else {writeLines(info,out('PARTIAL.txt'));msg("Scoring completed but some heatmaps were not generated; see figure_status.csv and omitted_rows.csv.")}
}
if(sys.nframe()==0L) main()

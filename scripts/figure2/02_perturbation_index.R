# Source: Figure 2 document, extracted text lines 475-1718
# Standalone script requiring a CK2-complete RDS and an external combined DEG table.
# See README.md and docs/KNOWN_LIMITATIONS.md before running.

# CK2-complete input variant: retain the four-component PI workflow and use the all-cell RDS with CK2 restored.
# Compare only Heat/Cold/PEG versus CK, retaining all annotated cell types; no TF or Endodermis-stage analysis here.
# Keep the external DEG file unchanged; changing the expression RDS does not recompute those DEGs.
# Do not overwrite the input RDS; reruns update results with the same PI_UserDEG_ prefix.
# Validate caches using input paths, file metadata, configuration and function signatures.
# Sourcing this script runs main once; do not call main again afterwards.
# R Console：source(file.path(Sys.getenv("SOYBEAN_PROJECT_ROOT", unset = "data"), "Soybean_UserDEG_PI_Figure4_v3_CK2.R"), encoding="UTF-8")
# Terminal：Rscript /data01/liuxy/Project/Soybean_UserDEG_PI_Figure4_v3_CK2.R

# The DEG component uses the specified filtered DEG table; other components are Augur, scDist and Mixscape.
# Retain the configured OUT_DIR; do not run DESeq2 or substitute Endodermis state contrasts for stress contrasts.
# The cited paper's 13 components differ from these four cell-type components; do not duplicate DEG counts to mimic four resolutions.
# New outputs use PI_UserDEG_ in the configured directory; existing PI caches are not overwritten.
# Usage from the repository root: source("scripts/figure2/02_perturbation_index.R", encoding = "UTF-8")
# =============================================================================
# 1. Configuration.
MERGED_RDS <- file.path(Sys.getenv("SOYBEAN_PROJECT_ROOT", unset = "data"), "CK2_added_v3_Results/seurat_obj_with_CK2_v3.rds")
ANALYSIS_MODE <- "annotated" # This variant uses existing cell types and four components; do not replicate a DEG count to fabricate 13 components.
# Keep the configured output directory independent of the new MERGED_RDS location.
OUT_DIR <- file.path(Sys.getenv("SOYBEAN_PROJECT_ROOT", unset = "data"), "CK2_added_v3_Results")
# Check sufficient cells in both replicates of each group, including zero-cell combinations in QC tables.
# Use content-signature caches when computation functions change; the output directory need not change.
RESUME <- TRUE
CELLTYPE_COL <- "cell_type"   # Annotation column in the CK2-complete RDS; retain all cell types instead of substituting P1/P2/P3.
CONDITION_COL <- "treatment"  # Select the treatment column explicitly to avoid conflicts with Group.
SAMPLE_COL <- "sample_id"     # The source notes this column passed its checks; sample counts are still validated and exported.
EXPECTED_REPLICATES <- 2L     # The source assumes two replicates per treatment; stop and report metadata inconsistencies.
PAIR_COL <- NULL              # Specify only for a genuine paired design; paired observations share this ID across treatments.
COVARIATE_COLS <- character() # Optional sample-level covariates, such as c("batch"); the design must remain full rank.
MIX_NEIGHBOR_BLOCK_COL <- NULL # Optional shared batch/paired blocks; each block must include CK.
CELLTYPE_KEEP <- NULL         # Retain all annotated cell types by default.
# Unify lateral-root label variants as Procambium.
# Normalize with canon() and match prefixes, including "Lateral root" and "Lateral_root".
# Also covers variants such as "Lateral root cells" and "Lateral roots".
# Apply the same mapping to expression metadata and the DEG table's Celltype column.
CELLTYPE_RENAME_MAP <- list(
  Procambium = c("Lateral root","Lateral Root","lateral root","Lateral_root","lateral_root",
                 "Lateral root cells","Lateral roots","lateral roots")
)
CONDITION_MAP <- list(
  CK=c("CK","Control","CTRL","REF","Reference","\u5bf9\u7167"),
  Heat=c("Heat","High_Temp","HighTemp","HighTemperature","HT","\u9ad8\u6e29"),
  Cold=c("Cold","Low_Temp","LowTemp","LowTemperature","LT","\u4f4e\u6e29"),
  PEG=c("PEG","PEG6000","PEG8000")
)
RESOLUTIONS <- c(0.1,0.2,0.4,0.8) # Retain the source resolution without forcing the number of clusters reported for rice.
ASSAY <- "RNA"
N_HVG <- 3000L
N_PCS <- 30L
GRAPH_REDUCTION <- NULL       # NULL computes RNA PCA; otherwise specify an existing reduction, such as harmony.
UMAP_REDUCTION <- "umap"      # Prefer the original embedding; compute UMAP only if unavailable.
SEED <- 20261021L
N_THREADS <- 4L               # Used by Augur; reduce to 2L if memory is limited.
MIN_CELLS_PER_SAMPLE_CLUSTER <- 20L
MIN_REPLICATES_PER_CONDITION <- 2L # At least two replicates; three or more are preferable.
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
DENSE_MEMORY_BUDGET_GB <- 32 # Estimated budget for major dense matrices, not automatic detection of available RAM.
SAVE_SEURAT_WITH_PI <- FALSE # TRUE also saves the reconstructed minimal object without overwriting the input RDS.
PIPELINE_VERSION <- "ensemble_pi_adapted_1.5_augur_future_fix"

# Retain legacy settings for cache validation; the new main does not run DESeq2 or access DESeq2_tables.
# Augur uses cell-level cross-validation by default, not leave-one-biological-sample-out validation; PI is neither significance nor damage probability.

# 2. Utilities and validation.
csv <- function(x,p) utils::write.csv(x,p,row.names=FALSE,na="",fileEncoding="UTF-8")
bind <- function(x) if(length(x)) do.call(rbind,x) else data.frame()
canon <- function(x) gsub("[[:space:]_.-]","",tolower(trimws(as.character(x))))
# Unify lateral-root-related cell-type labels as Procambium.
# Normalize with canon() before prefix matching to cover variants such as "Lateral root cells".
# Edit CELLTYPE_RENAME_MAP to change mappings; the function itself need not change.
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
    stop("Conflicting aliases in CONDITION_MAP.")
  unname(v[match(canon(x),k)])
}
pick_col <- function(md,explicit,candidates,label,required=TRUE,transform=as.character) {
  if(!is.null(explicit)) {
    if(length(explicit)!=1L || !explicit %in% names(md)) stop(label," does not exist: ",explicit)
    return(explicit)
  }
  hit<-names(md)[canon(names(md)) %in% canon(candidates)]
  if(length(hit)>1L) {
    v<-lapply(hit,function(h)transform(md[[h]]))
    if(!all(vapply(v[-1],identical,logical(1),v[[1]])))
      stop(label," conflicting candidate columns: ",paste(hit,collapse=", "),". Specify the column in the configuration.")
  }
  if(!length(hit)) {
    if(required) stop("Cannot identify ",label,"; inspect 00_metadata_inventory.csv and configure explicitly.")
    return(NULL)
  }
  hit[1]
}
same_partition <- function(a,b) {
  a<-trimws(as.character(a));b<-trimws(as.character(b))
  identical(match(a,unique(a)),match(b,unique(b)))
}
resolve_sample_column <- function(md,condition,explicit=NULL) {
  # Restrict md to analyzed cells so salt-treatment or unannotated cells do not affect sample-column selection.
  if(!nrow(md)) stop("No cells remain after treatment/cell-type filtering; cannot validate samples.")
  candidates<-c("sample_id","sampleID","biological_sample","sample","sample_name")
  hits<-unique(unlist(lapply(candidates,function(k)names(md)[canon(names(md))==canon(k)])))
  if(!is.null(explicit)) {
    if(length(explicit)!=1L || !explicit %in% names(md)) stop("Configured SAMPLE_COL does not exist: ",explicit)
    hits<-unique(c(explicit,hits))
  }
  # Use orig.ident only when no standard sample column exists; never silently override a conflicting sample column.
  if(!length(hits) && "orig.ident" %in% names(md)) hits<-"orig.ident"
  if(!length(hits)) stop("No candidate sample columns; inspect 00_metadata_inventory.csv.")
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
    message("Sample column: ",chosen,"；",why)
    chosen
  }
  if(!is.null(explicit)) {
    if(!explicit %in% eligible) stop("Configured sample column failed structural checks; see 00_sample_column_audit.csv.")
    return(decision(explicit,"user_explicit_selection"))
  }
  if(!length(eligible)) stop("All sample columns failed structural checks; see the sample audit and crosswalk CSVs. Biological replicates cannot be fabricated.")
  if(length(eligible)==1L) return(decision(eligible,"only_structurally_eligible_candidate"))
  equiv<-all(vapply(eligible[-1L],function(k)same_partition(md[[eligible[1L]]],md[[k]]),logical(1)))
  if(equiv) return(decision(eligible[1L],paste0("equivalent_cell_partition: ",paste(eligible,collapse=" = "))))
  csv(data.frame(selected_sample_column=NA_character_,reason="different_eligible_partitions_require_sample_design",
                 caveat="The metadata columns encode different sample groupings; choose using experimental sample records."),
      file.path(OUT_DIR,"00_sample_column_decision.csv"))
  stop("Candidate sample columns partition cells differently. Inspect 00_sample_column_crosswalk.csv and 00_sample_column_audit.csv against sample records, then set SAMPLE_COL.")
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
  # Rename lateral-root variants as Procambium before scope filtering so all subsequent counts use the same label.
  types<-apply_celltype_rename(trimws(as.character(old[[ct]])))
  scope<-!is.na(mapped)&!is.na(types)&nzchar(types)&
    !tolower(types) %in% c("undefined","unknown","unassigned","unclassified","na")
  if(!is.null(CELLTYPE_KEEP))scope<-scope & types %in% CELLTYPE_KEEP
  if(!all(c("CK","Heat","Cold","PEG") %in% mapped[scope]))
    stop("CK, Heat, Cold and PEG must be identified; inspect 00_condition_mapping.csv. Other treatments are excluded.")
  sa<-resolve_sample_column(old[scope,,drop=FALSE],mapped[scope],SAMPLE_COL)
  z<-data.frame(cell=rownames(old),pi_celltype=types,
                pi_condition=map_condition(old[[co]]),pi_sample=trimws(as.character(old[[sa]])),
                row.names=rownames(old),stringsAsFactors=FALSE)
  extra<-unique(c(PAIR_COL,COVARIATE_COLS,MIX_NEIGHBOR_BLOCK_COL))
  if(any(!extra %in% names(old))) stop("Covariate/pair/block column missing: ",paste(setdiff(extra,names(old)),collapse=", "))
  if(!is.null(PAIR_COL)) z$pi_pair<-as.character(old[[PAIR_COL]])
  if(length(COVARIATE_COLS)) for(i in seq_along(COVARIATE_COLS))
    z[[paste0("pi_cov",i)]]<-old[[COVARIATE_COLS[i]]]
  if(!is.null(MIX_NEIGHBOR_BLOCK_COL)) z$pi_mixblock<-as.character(old[[MIX_NEIGHBOR_BLOCK_COL]])
  keep<-!is.na(z$pi_condition)&!is.na(z$pi_celltype)&nzchar(z$pi_celltype)&
    !tolower(z$pi_celltype) %in% c("undefined","unknown","unassigned","unclassified","na")
  if(!is.null(CELLTYPE_KEEP)) keep<-keep & z$pi_celltype %in% CELLTYPE_KEEP
  z<-z[keep,,drop=FALSE]
  if(!nrow(z)||anyDuplicated(z$cell)) stop("No cells or duplicated barcodes.")
  if(anyNA(z)||any(!nzchar(z$pi_sample))) stop("Retained cells have missing metadata; complete the relevant columns first.")
  if(!all(c("CK","Heat","Cold","PEG") %in% z$pi_condition)) stop("CK, Heat, Cold and PEG must be identified; check the mapping table.")
  cols<-c("pi_sample","pi_condition",if(!is.null(PAIR_COL))"pi_pair",paste0("pi_cov",seq_along(COVARIATE_COLS)))
  # Handle paste0 with zero-length inputs explicitly.
  cols<-intersect(cols,names(z)); sm<-unique(z[,cols,drop=FALSE])
  csv(sm,file.path(OUT_DIR,"00_sample_mapping.csv"))
  if(anyDuplicated(sm$pi_sample)) stop("A sample ID maps to multiple treatment/covariate/pair values; correct the sample information.")
  if(nrow(sm)>nrow(z)/MIN_CELLS_PER_SAMPLE_CLUSTER) stop("Too many sample IDs; check for accidental use of cell barcodes.")
  nr<-table(factor(sm$pi_condition,levels=c("CK","Heat","Cold","PEG")))
  csv(data.frame(condition=names(nr),n_samples=as.integer(nr)),file.path(OUT_DIR,"00_sample_counts.csv"))
  control_counts<-as.data.frame(table(Sample=z$pi_sample[z$pi_condition=="CK"]),stringsAsFactors=FALSE)
  names(control_counts)<-c("Sample","Cells")
  csv(control_counts,file.path(OUT_DIR,"00_CK_control_sample_counts.csv"))
  message("CK samples and cell counts: ",paste(paste(control_counts$Sample,control_counts$Cells,sep="="),collapse=", "))
  message("Global sample counts: ",paste(paste(names(nr),as.integer(nr),sep="="),collapse=", "))
  if(!is.null(EXPECTED_REPLICATES) && any(nr!=EXPECTED_REPLICATES))
    stop("sample_id does not match two replicates per condition; see 00_sample_mapping.csv. Samples will not be merged or split automatically.")
  if(any(nr<MIN_REPLICATES_PER_CONDITION))
    stop("At least one condition lacks biological replicates for a complete PI; random cell splitting will not be used to fabricate replicates.")
  csv(data.frame(role=c("celltype","condition","sample"),column=c(ct,co,sa)),
      file.path(OUT_DIR,"00_resolved_columns.csv"))
  z$pi_condition<-factor(z$pi_condition,levels=c("CK","Heat","Cold","PEG"))
  z
}

raw_counts <- function(obj,cells) {
  if(!ASSAY %in% names(obj@assays)) stop("Raw RNA assay missing.")
  a<-obj[[ASSAY]]
  if(inherits(a,"Assay5")) {
    ls<-SeuratObject::Layers(a,search=NA);cl<-ls[grepl("^counts($|[.])",ls)]
    if(!length(cl)) stop("RNA counts layer missing.")
    if(length(cl)>1) {
      cc<-unlist(lapply(cl,function(k)colnames(SeuratObject::LayerData(a,layer=k))))
      if(anyDuplicated(cc)) stop("Counts layers contain duplicated cells and cannot be joined directly.")
      a<-SeuratObject::JoinLayers(a,layers="counts",new="counts");cl<-"counts"
    }
    x<-SeuratObject::LayerData(a,layer=cl)
  } else x<-SeuratObject::GetAssayData(a,slot="counts")
  if(!all(cells %in% colnames(x))) stop("Some cells lack RNA counts.")
  x<-x[,cells,drop=FALSE]
  if(!inherits(x,"sparseMatrix")) x<-methods::as(x,"dgCMatrix")
  if(is.null(rownames(x))||anyDuplicated(rownames(x))) stop("Gene IDs are missing or duplicated.")
  if(any(!is.finite(x@x))||any(x@x<0)||any(abs(x@x-round(x@x))>1e-6))
    stop("Raw nonnegative integer counts are required, not data/SCT/integrated values.")
  if(any(Matrix::colSums(x)==0)) stop("Some cells have zero total counts.")
  x[Matrix::rowSums(x)>0,,drop=FALSE]
}
resolution_key <- function(r) paste0("r",gsub("[.]","p",format(r,trim=TRUE,scientific=FALSE)))
analysis_groups <- function() {
  if(identical(ANALYSIS_MODE,"annotated"))
    return(data.frame(grouping="annotation",column="pi_celltype",resolution=NA_real_))
  if(!identical(ANALYSIS_MODE,"multi_resolution"))stop("Unknown ANALYSIS_MODE.")
  rk<-vapply(RESOLUTIONS,resolution_key,character(1))
  data.frame(grouping=rk,column=paste0("pi_",rk),resolution=RESOLUTIONS)
}
memory_check <- function(ng,nc,multiplier=3) {
  # Convert dimensions to double before multiplication to avoid integer overflow above 2^31 - 1.
  dims<-as.double(c(ng,nc,multiplier))
  if(length(dims)!=3L || any(!is.finite(dims)) || any(dims<0))stop("Invalid dimensions for memory estimation.")
  gb<-prod(dims)*8/1024^3
  if(!is.finite(gb))stop("Memory estimate exceeds the representable range.")
  if(gb>DENSE_MEMORY_BUDGET_GB) stop(sprintf("Estimated dense-matrix workload %.1f GB exceeds budget %.1f GB; adjust the budget or server. Genes will not be silently removed.",gb,DENSE_MEMORY_BUDGET_GB))
  invisible(gb)
}
get_expression <- function(obj,assay="RNA",layer="data") {
  if(utils::packageVersion("SeuratObject") >= "5.0.0")
    SeuratObject::GetAssayData(obj,assay=assay,layer=layer)
  else SeuratObject::GetAssayData(obj,assay=assay,slot=layer)
}
# 3. Preprocessing and shared clustering.
prepare <- function() {
  original<-readRDS(MERGED_RDS)
  if(!inherits(original,"Seurat")) stop("RDS does not contain a Seurat object.")
  md<-get_metadata(original);counts<-raw_counts(original,md$cell)
  # Reconstruct a minimal RNA object to avoid using unidentified integrated/data layers as expression input.
  obj<-SeuratObject::CreateSeuratObject(counts=counts,meta.data=md,assay="RNA")
  obj<-Seurat::NormalizeData(obj,normalization.method="LogNormalize",scale.factor=10000,verbose=FALSE)
  obj<-Seurat::FindVariableFeatures(obj,selection.method="vst",nfeatures=N_HVG,verbose=FALSE)
  hv<-SeuratObject::VariableFeatures(obj)
  np<-min(N_PCS,length(hv)-1L,ncol(obj)-1L)
  if(np<5) stop("Too few cells/features for stable PCA.")
  obj<-Seurat::ScaleData(obj,features=hv,verbose=FALSE)
  obj<-Seurat::RunPCA(obj,features=hv,npcs=np,seed.use=SEED,verbose=FALSE)
  red<-"pca"
  if(ANALYSIS_MODE=="multi_resolution" && !is.null(GRAPH_REDUCTION)) {
    if(!GRAPH_REDUCTION %in% names(original@reductions)) stop("GRAPH_REDUCTION is absent from the original object.")
    emb<-SeuratObject::Embeddings(original[[GRAPH_REDUCTION]])
    if(!all(colnames(obj) %in% rownames(emb))) stop("Specified reduction lacks some analyzed cells.")
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
    if(!all(colnames(obj) %in% rownames(em))) stop("Original UMAP lacks cells; set UMAP_REDUCTION <- NULL to recompute.")
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
      stop("Multiple samples share a pair/condition within a cluster; verify the design or combine technical replicates.")
    pairs<-intersect(m$pi_pair[m$pi_condition=="CK"],m$pi_pair[m$pi_condition==st])
    m<-m[m$pi_pair %in% pairs,,drop=FALSE]
  }
  tb<-table(factor(m$pi_condition,levels=c("CK",st)))
  if(any(tb<MIN_REPLICATES_PER_CONDITION)) stop("Samples passing within-group QC: CK=",tb[1],"；",st,"=",tb[2],
                                                "; required per group: ",MIN_REPLICATES_PER_CONDITION,". This is not the full experiment replicate count; see QC_sample_group_counts.csv.")
  m$pi_condition<-factor(as.character(m$pi_condition),levels=c("CK",st))
  if(!is.null(PAIR_COL)) m$pi_pair<-factor(m$pi_pair)
  cov<-if(length(COVARIATE_COLS))paste0("pi_cov",seq_along(COVARIATE_COLS)) else character()
  # Omit covariates constant within a comparison; report covariates collinear with treatment.
  cov<-cov[vapply(cov,function(k)length(unique(m[[k]]))>1,logical(1))]
  for(k in cov) if(is.character(m[[k]])||is.factor(m[[k]])) m[[k]]<-factor(m[[k]])
  form<-stats::reformulate(c(if(!is.null(PAIR_COL))"pi_pair",cov,"pi_condition"))
  mm<-stats::model.matrix(form,m)
  if(qr(mm)$rank<ncol(mm) || nrow(mm)<=ncol(mm)) stop("Rank-deficient design or no residual degrees of freedom; check batch/treatment confounding.")
  list(meta=m,formula=form)
}
run_augur <- function(norm,md,col,st) {
  take<-md$pi_condition %in% c("CK",st);m<-droplevels(md[take,,drop=FALSE])
  tab<-table(as.character(m[[col]]),factor(m$pi_condition,levels=c("CK",st)))
  eligible<-rownames(tab)[apply(tab,1,min)>=AUGUR_SUBSAMPLE_SIZE]
  m<-m[as.character(m[[col]]) %in% eligible,,drop=FALSE]
  if(!nrow(m)) stop("Augur: no group meets the minimum cells per treatment.")
  x<-norm[,rownames(m),drop=FALSE];x<-x[Matrix::rowSums(x)>0,,drop=FALSE]
  message(sprintf("Augur %s: cells=%d, genes=%d, celltypes=%d, n_subsamples=%d, subsample_size=%d, folds=%d, threads=%d",
                  st, ncol(x), nrow(x), length(eligible),
                  AUGUR_N_SUBSAMPLES, AUGUR_SUBSAMPLE_SIZE, AUGUR_FOLDS, N_THREADS))
  # Augur::calculate_auc uses future. Under RStudio or RStudio Server,
  # a global multiprocess plan may use multicore on Linux, where forked workers can be unstable.
  # A worker crash may produce only "Execution halted" in the main R process.
  # Use multisession for this computation and restore the previous plan afterwards.
  if(!requireNamespace("future",quietly=TRUE)) stop("Augur requires future; install it first.")
  old_plan<-future::plan()
  on.exit(future::plan(old_plan),add=TRUE)
  if(N_THREADS>1L) future::plan(future::multisession,workers=N_THREADS)
  else future::plan(future::sequential)
  # Adapt arguments to the installed Augur version; older versions lack n_threads.
  auc_args<-names(formals(Augur::calculate_auc))
  call_args<-list(input=x,meta=m,cell_type_col=col,label_col="pi_condition",
                  n_subsamples=AUGUR_N_SUBSAMPLES,subsample_size=AUGUR_SUBSAMPLE_SIZE,
                  folds=AUGUR_FOLDS,min_cells=AUGUR_SUBSAMPLE_SIZE,show_progress=TRUE)
  if("n_threads" %in% auc_args) call_args$n_threads<-N_THREADS
  out<-tryCatch(
    do.call(Augur::calculate_auc,call_args),
    error=function(e){
      msg<-sprintf("Augur calculate_auc failed (cells=%d, genes=%d, threads=%d): %s",
                   ncol(x),nrow(x),N_THREADS,conditionMessage(e))
      message(msg)
      stop(msg)
    }
  )
  if(!all(c("cell_type","auc") %in% names(out$AUC))) stop("Unexpected Augur output structure; check the installed version.")
  # Do not duplicate the original matrix in checkpoints.
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
  if(any(nr<MIN_REPLICATES_PER_CONDITION)) stop("scDist samples after within-group cell filtering: CK=",nr[1],"；",st,"=",nr[2],
                                                "; required per group: ",MIN_REPLICATES_PER_CONDITION,". See QC_sample_group_counts.csv.")
  m$pi_condition<-factor(as.character(m$pi_condition),levels=c("CK",st));m$pi_sample<-factor(m$pi_sample)
  cov<-if(length(COVARIATE_COLS))paste0("pi_cov",seq_along(COVARIATE_COLS)) else character()
  cov<-cov[vapply(cov,function(k)length(unique(m[[k]]))>1,logical(1))]
  for(k in cov) if(is.character(m[[k]])||is.factor(m[[k]]))m[[k]]<-factor(m[[k]])
  fx<-c("pi_condition",cov);re<-c("pi_sample",if(!is.null(PAIR_COL))"pi_pair")
  mm<-model.matrix(reformulate(fx),m)
  if(qr(mm)$rank<ncol(mm)) stop("scDist fixed-effects design is rank deficient.")
  list(meta=m,fixed=fx,random=re,nr=nr)
}
run_scdist <- function(norm,md,col,cluster,st) {
  design<-scdist_metadata(md,col,cluster,st)
  m<-design$meta;fx<-design$fixed;re<-design$random;nr<-design$nr
  x<-norm[,rownames(m),drop=FALSE];x<-x[Matrix::rowSums(x)>0,,drop=FALSE]
  if(min(dim(x))<=SCDIST_D) stop("Insufficient scDist dimensions; d will not be changed automatically.")
  memory_check(nrow(x),ncol(x))
  # Call scDist separately for each group to avoid densifying the entire single-cell matrix at once.
  m$pi_sccluster<-"cluster"
  out<-scDist::scDist(normalized_counts=as.matrix(x),meta.data=m,fixed.effects=fx,
                      random.effects=re,clusters="pi_sccluster",d=SCDIST_D,min.count.per.cell=SCDIST_MIN_CELLS)
  rr<-out$results
  if(nrow(rr)!=1L || !"Dist." %in% names(rr) || !is.finite(rr[["Dist."]][1]))
    stop("scDist returned no valid distance; it cannot be replaced with zero.")
  list(value=as.numeric(rr[["Dist."]][1]),results=rr,n_CK=nr[1],n_stress=nr[2])
}

# Build the full group-by-real-sample grid; zero-cell entries are for QC, not fabricated zero-expression samples.
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
  message("Skipping ",key,": ",reason)
  list(ok=FALSE,value=NULL,error=reason,warnings=character(),key=key,status="skipped_QC")
}

# 4. Single-cell Mixscape.
# When filling missing genes, Seurat retrieves a stored Neighbor from Tools using the Command name.
# Some versions use inconsistent names, leading to Indices(NULL) errors.
# Add aliases only to the current CalcPerturbSig object after checking barcodes and indices; do not alter neighbors.
select_mixscape_neighbors <- function(tool_list,command,rep_cells,rep_controls,k,indices_fun) {
  expected<-make.names(paste0("PRTB_",names(rep_cells)))
  if(anyDuplicated(expected))stop("Mixscape block names collide after normalization.")
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
  if(!length(good))stop("No Neighbor cache matches these cells, CK pool and neighbor count; Indices(NULL) will not be called.")
  if(length(good)>1L && !all(vapply(good[-1L],function(x)identical(tool_list[[good[1L]]],tool_list[[x]]),logical(1))))
    stop("Multiple inconsistent Neighbor caches found; cannot select automatically.")
  source<-if(command %in% good)command else good[1L]
  list(bundle=tool_list[[source]],source=source,command=command,aliased=!identical(source,command))
}
rebuild_mixscape_neighbors <- function(embedding,rep_cells,rep_controls,ndims,k,nn_helper) {
  # If the cache is absent, reconstruct using the same PCA, CK pool, RANN and k as CalcPerturbSig's NNHelper call.
  if(ncol(embedding)<ndims)stop("Insufficient PCA dimensions to reconstruct Mixscape neighbors.")
  names_out<-make.names(paste0("PRTB_",names(rep_cells)))
  if(anyDuplicated(names_out))stop("Mixscape block names collide after normalization.")
  out<-lapply(seq_along(rep_cells),function(i) {
    if(length(rep_controls[[i]])<k)stop("Too few CK cells to reconstruct neighbors.")
    if(!all(c(rep_cells[[i]],rep_controls[[i]]) %in% rownames(embedding)))stop("PCA lacks required cells.")
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
  if(length(candidates)!=1L)stop("The current CalcPerturbSig command must be uniquely identified; see Mixscape command diagnostics.")
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
    message("Neighbor cache unusable; reconstructing with the same PCA/CK pool: ",failure)
    helper<-getFromNamespace("NNHelper","Seurat")
    if(!all(c("data","query","k","method") %in% names(formals(helper))))
      stop("Incompatible Seurat NNHelper interface; command diagnostics saved.")
    bundle<-rebuild_mixscape_neighbors(SeuratObject::Embeddings(sub[["pca"]]),rep_cells,rep_controls,MIX_NDIMS,MIX_NEIGHBORS,helper)
    selected<-select_mixscape_neighbors(setNames(list(bundle),candidates),candidates,rep_cells,rep_controls,MIX_NEIGHBORS,SeuratObject::Indices)
    rebuilt<-TRUE
  }
  # Tool<- derives its name from the calling function; assign an explicit alias in the current tools slot instead.
  sub@tools[[selected$command]]<-selected$bundle
  csv(data.frame(condition=st,source_tool=selected$source,lookup_command=selected$command,
                 action=if(rebuilt)"rebuilt_same_PCA_CK_RANN" else if(selected$aliased)"added_validated_alias" else "already_aligned",
                 original_issue=if(is.null(failure))"" else failure,
                 n_blocks=length(reps),neighbors=MIX_NEIGHBORS),file.path(OUT_DIR,paste0("Mixscape_",st,"_neighbor_audit.csv")))
  message("Mixscape neighbor validation: ",selected$source," -> ",selected$command)
  sub
}
run_mixscape <- function(obj,st) {
  cells<-rownames(obj[[]])[obj$pi_condition %in% c("CK",st)]
  sub<-subset(obj,cells=cells)
  sub$pi_condition<-factor(as.character(sub$pi_condition),levels=c("CK",st))
  if(sum(sub$pi_condition=="CK")<MIX_NEIGHBORS) stop("Too few CK cells for Mixscape.")
  block<-if(is.null(MIX_NEIGHBOR_BLOCK_COL)) NULL else "pi_mixblock"
  if(!is.null(block)) {
    tb<-table(sub[[block]][,1],factor(sub$pi_condition,levels=c("CK",st)))
    if(any(tb[,"CK"]<MIX_NEIGHBORS)) stop("Each Mixscape neighbor block requires sufficient CK cells; do not group by treatment-specific sample_id.")
  }
  if(ncol(SeuratObject::Embeddings(sub[["pca"]]))<MIX_NDIMS) stop("Insufficient Mixscape PCA dimensions.")
  hv<-intersect(SeuratObject::VariableFeatures(sub),rownames(sub[["RNA"]]))
  if(!length(hv))stop("Mixscape lacks RNA variable genes; cannot automatically use all genes.")
  # Initial PRTB is computed only for hv; using all-gene dimensions would overestimate memory and risk overflow.
  # This estimates the initial matrix, not peak memory for additional genes or all Seurat intermediates.
  memory_check(length(hv),ncol(sub),multiplier=4)
  sub<-Seurat::CalcPerturbSig(sub,assay="RNA",features=hv,slot="data",
                              gd.class="pi_condition",nt.cell.class="CK",reduction="pca",ndims=MIX_NDIMS,
                              num.neighbors=MIX_NEIGHBORS,split.by=block,new.assay.name="PRTB",verbose=TRUE)
  sub<-repair_mixscape_neighbors(sub,st,block)
  # Run Mixscape once per stress-versus-CK comparison, not once per clustering resolution.
  # RunMixscape fills and scales the required PRTB genes internally; do not LogNormalize negative PRTB values.
  prior_tools<-mixscape_tool_names(sub@tools)
  if(length(prior_tools))sub@tools[prior_tools]<-NULL # Remove only previous results within the current subset to avoid retrieving historical pvec values.
  RunMixscape<-Seurat::RunMixscape
  sub<-RunMixscape(sub,assay="PRTB",slot="scale.data",labels="pi_condition",
                   nt.class.name="CK",new.class.name="pi_mixscape",min.de.genes=MIX_MIN_DE_GENES,
                   min.cells=5,de.assay="RNA",logfc.threshold=MIX_LOGFC_THRESHOLD,
                   iter.num=MIX_ITERATIONS,split.by=NULL,verbose=TRUE)
  sm<-sub[[]];cells<-colnames(sub)
  tools_found<-sub@tools[mixscape_tool_names(sub@tools)]
  # Save the lightweight raw output before extraction so pvec survives an extraction error.
  snapshot<-list(condition=st,cells=cells,metadata=sm[cells,"pi_condition",drop=FALSE],
                 tools=tools_found,all_tool_names=names(sub@tools),compute_signature=COMPUTE_SIGNATURE,
                 classification=sm[,grep("^pi_mixscape",names(sm)),drop=FALSE])
  saveRDS(snapshot,pi_path(paste0("Mixscape_",st,"_tool_snapshot.rds")))
  capture.output(str(tools_found,max.level=3),file=pi_path(paste0("Mixscape_",st,"_tool_structure.txt")))
  selected<-select_mixscape_pvec(tools_found,st,cells,as.character(sm[cells,"pi_condition"]))
  pi_csv(selected$audit,paste0("Mixscape_",st,"_tool_selection"))
  val<-selected$value
  message("Mixscape pvec source: ",selected$source,"; valid stress cells: ",
          sum(is.finite(val[as.character(sm[cells,"pi_condition"])==st])))
  list(value=val,tools=selected$tool,tool_source=selected$source,
       classification=sm[,grep("^pi_mixscape",names(sm)),drop=FALSE])
}


# User DEG inputs and plotting configuration.
PAPER_DIR <- file.path(Sys.getenv("SOYBEAN_PROJECT_ROOT", unset = "data"), "muscat_results_ck_ref_paper")
# Use the explicitly configured input; do not substitute another DEG table or old DESeq2 results if absent.
DEG_FILE <- file.path(PAPER_DIR,'All_DEG_filtered_paper_thresholds.csv')
# To use combined.csv instead, change the single input path; do not add the two tables together.
DEG_COLUMNS <- list(gene=NULL,celltype=NULL,contrast=NULL,condition=NULL,state=NULL)
DEG_CONDITION_IS_VS_CK <- FALSE # Set TRUE only if Heat/Cold/PEG labels in the table explicitly represent contrasts against CK.
DEG_CELLTYPE_MAP <- character() # Optional mapping c('DEG_label' = 'Seurat_cell_type'), applied after the standard cell-type renaming.
ZERO_DEG_COMPARISONS <- NULL    # Optional data.frame(Celltype, Condition) of verified zero-DEG comparisons; unknown comparisons are not zero-filled.
EXPORT_TIFF <- FALSE
PLOT_POINT_SIZE <- 0.13
VIOLIN_POINT_MAX_PER_GROUP <- 500L # Limit only overlaid points; violin densities and statistics use all cells with complete PI values.
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
  if(!all(cols %in% names(d)))stop(label,"Missing columns: ",paste(setdiff(cols,names(d)),collapse=', '))
}
unique_rows <- function(d,cols,label) {
  need(d,cols,label)
  if(anyNA(d[cols]) || anyDuplicated(d[cols]))stop(label,"Missing or duplicated keys: ",paste(cols,collapse=' / '))
}
digest_object <- function(x) {
  f<-tempfile(fileext='.rds');on.exit(unlink(f));saveRDS(x,f,version=2)
  unname(tools::md5sum(f))
}
read_csv_pi <- function(path) {
  d<-read.csv(path,check.names=FALSE,stringsAsFactors=FALSE,fileEncoding='UTF-8-BOM')
  if(anyDuplicated(names(d)))stop("Duplicated CSV columns: ",path)
  d
}
deg_column <- function(d,role,aliases,required=FALSE,multiple=FALSE) {
  explicit<-DEG_COLUMNS[[role]]
  if(!is.null(explicit)) {
    if(!all(explicit %in% names(d)))stop('DEG_COLUMNS$',role,"Specified column does not exist.")
    if(!multiple && length(explicit)!=1L)stop("Only one column may be specified for this role: ",role)
    return(explicit)
  }
  hit<-names(d)[canon(names(d)) %in% canon(aliases)]
  if(!length(hit)) {
    if(required)stop("Cannot identify DEG ",role," column; set DEG_COLUMNS. Available columns: ",paste(names(d),collapse=', '))
    return(NULL)
  }
  if(multiple)return(hit)
  if(length(hit)>1L) {
    v<-lapply(hit,function(n)trimws(as.character(d[[n]])))
    if(!all(vapply(v[-1],identical,logical(1),v[[1]])))
      stop("DEG columns contain conflicting values; set DEG_COLUMNS$",role," explicitly: ",paste(hit,collapse=', '))
  }
  hit[1]
}
parse_ck_pair <- function(v) {
  # Accept explicit stress-versus-CK contrasts; reversed contrasts may be counted because PI does not use logFC direction.
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
  # Apply the same lateral-root-to-Procambium renaming to the DEG Celltype column,
  # so DEG counts align with expression-object cell types.
  d<-data.frame(Gene=trimws(as.character(raw[[gc]])),
                Celltype=apply_celltype_rename(trimws(as.character(raw[[cc]]))),
                Condition=NA_character_,State=if(is.null(sc))NA_character_ else trimws(as.character(raw[[sc]])),
                Status='Unresolved_comparison',Source_row=seq_len(n),stringsAsFactors=FALSE)
  if(anyNA(d$Gene) || any(!nzchar(d$Gene)))stop("DEG table contains empty gene IDs.")
  for(i in seq_len(n)) {
    labs<-if(length(xc))trimws(as.character(unlist(raw[i,xc,drop=FALSE]))) else character()
    labs<-labs[!is.na(labs) & nzchar(labs)]
    pars<-parse_ck_pair(labs)
    # Do not let another column silently override explicit stage-level or treatment-versus-treatment contrasts.
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
  if(!file.exists(DEG_FILE))stop("Specified DEG file missing: ",DEG_FILE,"; the input will not be replaced automatically.")
  # Even with a validated standardized cache, read the original table's row count and columns rather than inferring them from its filename.
  raw<-read_csv_pi(DEG_FILE)
  pi_csv(data.frame(Column=names(raw)),'DEG_input_columns')
  md5<-unname(tools::md5sum(DEG_FILE));reader<-'Direct_user_DEG'
  saved<-file.path(OUT_DIR,'Nature_DE_filtered_annotated.csv')
  prov<-file.path(OUT_DIR,'Nature_DE_input_status.csv')
  d<-NULL
  # Parse explicitly configured columns directly so an existing standardized cache cannot mask new settings.
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
          unique_rows(q,'Source_row',"Nature-standardized DEG table")
          if(!setequal(q$Source_row,seq_len(nrow(raw))))stop("Standardized DEG source row numbers are incomplete.")
          d<-q[c('Gene','Celltype','State','Condition','Source_row')]
          # Reapply the same renaming to legacy cached labels such as Lateral root.
          d$Celltype<-apply_celltype_rename(d$Celltype)
          conflict<-as.character(q$Contrast_conflict)
          if(anyNA(conflict)||any(!conflict %in% c('TRUE','FALSE')))stop("Invalid standardized DEG conflict field.")
          d$Status<-ifelse(conflict=='TRUE','Conflicting_or_nonCK_comparison',
                           ifelse(d$Condition %in% STRESSES,'Verified_standardized_user_DEG','Unresolved_comparison'))
          d$Condition[conflict=='TRUE']<-NA_character_
          reader<-'Verified_Nature_annotation_of_same_user_DEG'
        }
      }
    }
  }
  if(is.null(d))d<-standardize_user_deg(raw)
  # Do not mix stage-level and whole-cell-type DEGs or infer stage-to-cell-type mappings from names.
  state<-!is.na(d$State) & nzchar(trimws(d$State))
  d$Status[state]<-'State_level_comparison_excluded'
  d$Condition[state]<-NA_character_
  if(length(DEG_CELLTYPE_MAP)) {
    if(is.null(names(DEG_CELLTYPE_MAP)) || anyDuplicated(names(DEG_CELLTYPE_MAP)))stop("DEG_CELLTYPE_MAP requires unique original labels.")
    k<-match(d$Celltype,names(DEG_CELLTYPE_MAP));yes<-!is.na(k)
    d$Celltype[yes]<-unname(DEG_CELLTYPE_MAP[k[yes]])
  }
  pi_csv(d,'DEG_row_audit')
  pi_csv(data.frame(File=DEG_FILE,MD5=md5,Rows=nrow(raw),Reader=reader,
                    Threshold='Already_filtered_by_user; no new FDR/logFC filter',
                    Excluded_state_rows=sum(state)),'DEG_source')
  if(!any(d$Condition %in% STRESSES))stop("No verified whole-cell-type stress-versus-CK comparisons in the DEG table; see DEG_row_audit. No other DE table will be substituted.")
  list(data=d,md5=md5)
}
count_user_degs <- function(deg,celltypes) {
  # Deduplicate exact gene IDs without removing soybean version suffixes or mixing transcript and gene identifiers.
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
    unique_rows(ZERO_DEG_COMPARISONS,c('Celltype','Condition'),"Explicit zero-DEG comparisons")
    for(i in seq_len(nrow(ZERO_DEG_COMPARISONS))) {
      j<-which(grid$Celltype==ZERO_DEG_COMPARISONS$Celltype[i] & grid$Condition==ZERO_DEG_COMPARISONS$Condition[i])
      if(length(j)!=1L || is.finite(grid$DEG_count[j]))stop("Explicit zero-DEG records are absent or conflict with nonzero counts.")
      grid$DEG_count[j]<-0L;grid$Status[j]<-'User_confirmed_zero_DEG'
    }
  }
  grid
}

# Auditable caches without changing directories or reusing old DESeq2 counts.
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
  need(raw,c('cell',PI_COMPONENTS[-1]),"Original non-DEG components")
  need(md,c('cell','pi_condition','pi_celltype','pi_sample'),"Original cell metadata")
  unique_rows(raw,'cell',"Original components");unique_rows(md,'cell',"Original metadata");unique_rows(um,'cell',"Original UMAP")
  if(!setequal(raw$cell,md$cell) || !all(md$cell %in% um$cell))stop("Original cache cell sets do not match.")
  raw<-raw[match(md$cell,raw$cell),c('cell',PI_COMPONENTS[-1]),drop=FALSE]
  if(!all(md$pi_condition %in% STRESSES))stop("Original PI metadata includes non-stress conditions.")
  # PI values in the legacy 05 table are not used in the new calculation.
  if(any(!is.finite(as.matrix(raw[PI_COMPONENTS[-1]]))))return(NULL)
  message("Reusing manifest-validated per-cell Augur/scDist/Mixscape values; ignoring old DESeq2 counts and PI values.")
  list(cells=md,nondeg=raw,umap=um,source='Validated_existing_component_CSVs')
}
compute_other_components <- function() {
  deps<-c('Seurat','SeuratObject','Matrix','Augur','scDist','mixtools')
  miss<-deps[!vapply(deps,requireNamespace,logical(1),quietly=TRUE)]
  if(length(miss))stop("Missing packages for other components: ",paste(miss,collapse=', '))
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
      if(nrow(chk)!=1L)stop("scDist QC records are not unique.")
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

# Combine four components and arrange the Figure 4-style panels.
combine_user_pi <- function(raw) {
  raw<-as.matrix(raw)
  if(!identical(colnames(raw),PI_COMPONENTS))stop("Incorrect PI component names or order.")
  if(any(raw[,'UserDEG_count']<0,na.rm=TRUE))stop("DEG counts cannot be negative.")
  if(any(raw[,'Augur_annotation']<0 | raw[,'Augur_annotation']>1,na.rm=TRUE))stop("Augur AUC must lie in [0, 1].")
  # Mixscape pvec contains raw projection scores, not posterior probabilities; values may fall outside [0, 1].
  # Do not clip raw model components; min-max scale each component using its observed pooled range below.
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
  if(anyNA(ix))stop("UMAP lacks PI cells.")
  d<-cm;d$UMAP1<-um$UMAP1[ix];d$UMAP2<-um$UMAP2[ix]
  if(any(!is.finite(as.matrix(d[c('UMAP1','UMAP2')]))))stop("UMAP contains non-finite coordinates.")
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
  # The fixed random seed controls point subsampling only; violins still use all complete cells.
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
  # Combine panels with base grid without adding patchwork/cowplot dependencies.
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

# v3: Mixscape result extraction and auditable resumption.
# SeuratObject::Tool<- names tools using the call stack; Seurat::RunMixscape may save results
# under Seurat::RunMixscape or Seurat..RunMixscape rather than RunMixscape.
# Read only pvec from the current RunMixscape result; do not substitute posterior probabilities or other components.
mixscape_tool_names <- function(x) {
  n<-names(x)
  n[grepl('(^|[^[:alnum:]])RunMixscape($|[^[:alnum:]])',n)]
}
extract_mix_pvec <- function(tool_result,st,cells) {
  if(anyNA(cells) || anyDuplicated(cells))stop("Mixscape input barcodes are not unique.")
  if(!is.list(tool_result) || is.null(tool_result[[st]]))
    stop("Mixscape output lacks treatment key ",st,"; actual keys: ",paste(names(tool_result),collapse=', '),
         ". Inspect the current tool_snapshot.rds; this alone does not establish insufficient DEGs.")
  parts<-tool_result[[st]]
  if(is.data.frame(parts))parts<-list(parts)
  if(!is.list(parts) || !length(parts))stop("Empty Mixscape treatment output: ",st)
  frames<-lapply(parts,function(d) {
    if(!is.data.frame(d) || !'pvec' %in% names(d) || !is.numeric(d$pvec))
      stop("Mixscape output is not a data frame with numeric pvec; posterior probability will not be substituted.")
    bc<-rownames(d)
    if(is.null(bc) || anyNA(bc) || anyDuplicated(bc) || any(!nzchar(bc)))
      stop("Mixscape pvec lacks unique cell barcodes.")
    if(!all(bc %in% cells))stop("Mixscape pvec contains barcodes outside the current comparison; prefixes/suffixes will not be removed automatically.")
    data.frame(cell=bc,value=d$pvec,stringsAsFactors=FALSE)
  })
  z<-bind(unname(frames))
  if(anyDuplicated(z$cell))stop("Mixscape result blocks contain duplicated cells; automatic averaging is not allowed.")
  val<-z$value[match(cells,z$cell)];names(val)<-cells
  if(!any(is.finite(val)))stop("Mixscape pvec has no finite values.")
  val
}
select_mixscape_pvec <- function(tool_list,st,cells,condition) {
  if(length(condition)!=length(cells))stop("Mixscape cell and treatment-label lengths differ.")
  candidates<-mixscape_tool_names(tool_list)
  if(!length(candidates))stop("No current RunMixscape tool slot found. Actual slots: ",paste(names(tool_list),collapse=', '))
  attempts<-lapply(candidates,function(k)tryCatch({
    v<-extract_mix_pvec(tool_list[[k]],st,cells)
    parts<-tool_list[[k]][[st]]
    if(is.data.frame(parts))parts<-list(parts)
    for(d in parts)if('pi_condition' %in% names(d)) {
      actual<-as.character(condition[match(rownames(d),cells)])
      if(anyNA(actual) || !identical(as.character(d$pi_condition),actual))
        stop("pvec treatment labels disagree with current metadata.")
    }
    if(!any(is.finite(v[as.character(condition)==st])))stop("pvec contains no valid stress-cell values.")
    list(ok=TRUE,value=v,error='')
  },error=function(e)list(ok=FALSE,value=NULL,error=conditionMessage(e))))
  audit<-data.frame(Tool=candidates,OK=vapply(attempts,`[[`,logical(1),'ok'),
                    Error=vapply(attempts,`[[`,character(1),'error'),stringsAsFactors=FALSE)
  good<-which(audit$OK)
  if(!length(good))stop("Failed to read RunMixscape slot: ",paste(paste(audit$Tool,audit$Error,sep=': '),collapse=' | '))
  if(length(good)>1L && !all(vapply(good[-1L],function(i)
    identical(attempts[[good[1L]]]$value,attempts[[i]]$value),logical(1))))
    stop("Multiple RunMixscape slots contain different pvec values; refusing to guess.")
  chosen<-if('RunMixscape' %in% candidates[good])match('RunMixscape',candidates) else good[1L]
  list(value=attempts[[chosen]]$value,tool=tool_list[[candidates[chosen]]],source=candidates[chosen],audit=audit)
}

# Cache migration may ignore these Mixscape-function changes only for non-Mixscape tasks.
# All other inputs, configuration, package versions and computation functions must still match.
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
  # Preserve the previous manifest before updating the entry file so prior successful caches remain verifiable.
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
  if(!is.list(value))stop("Task did not return the expected list: ",key)
  if(grepl('^Augur_',key)) {
    v<-value$value;unique_rows(v,'cluster',"Augur cache")
    if(!is.numeric(v$value) || !any(is.finite(v$value)) ||
       any(v$value[is.finite(v$value)]<0 | v$value[is.finite(v$value)]>1))
      stop("Augur returned no valid AUC.")
  } else if(grepl('^scDist_',key)) {
    if(!is.numeric(value$value) || length(value$value)!=1L || !is.finite(value$value))
      stop("scDist did not return one valid distance.")
  } else if(grepl('^Mixscape_',key)) {
    v<-value$value
    if(!is.numeric(v) || is.null(names(v)) || anyDuplicated(names(v)) || !any(is.finite(v)))
      stop("Mixscape returned no valid pvec with unique barcodes.")
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
      if(is.null(issue)) {z$cache_source<-path;message("Reusing successful cache: ",key);return(z)}
      message("Invalid cached result; recomputing ",key,'：',issue)
    }
  }
  set.seed(job_seed(key));message("Computing: ",key);warns<-character();error_call<-''
  z<-tryCatch(withCallingHandlers({
    v<-fun();validate_pi_job(key,v);list(ok=TRUE,value=v,error='')
  },warning=function(w) {warns<<-c(warns,conditionMessage(w));invokeRestart('muffleWarning')}),
  error=function(e) {
    error_call<<-paste(deparse(conditionCall(e)),collapse=' ')
    list(ok=FALSE,value=NULL,error=conditionMessage(e))
  })
  z$key<-key;z$warnings<-unique(warns);z$cache_source<-'';z$error_call<-error_call
  saveRDS(z,fresh)
  if(!z$ok)message("Failed [",key,']：',z$error)
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
  for(j in seq_len(ncol(raw)))message("Valid component ",colnames(raw)[j],'：',sum(is.finite(raw[,j])),'/',nrow(raw))
  invisible(tab)
}


main <- function() {
  if(ANALYSIS_MODE!='annotated')stop("User cell-type DEGs support only the annotated four-component PI. Without DEGs at four resolutions, the paper's 13-component PI cannot be constructed.")
  if(!dir.exists(OUT_DIR))dir.create(OUT_DIR,recursive=TRUE,showWarnings=FALSE)
  if(!file.exists(MERGED_RDS))stop("Original expression RDS missing: ",MERGED_RDS)
  if(!requireNamespace('ggplot2',quietly=TRUE))stop("ggplot2 is missing.")
  status<-'Running';stage<-'DEG_input'
  pi_csv(data.frame(Status=status,Stage=stage,Version=USER_PI_VERSION,Time=as.character(Sys.time())),'run_status')
  on.exit({
    if(status=='Running')status<-'Error_or_interrupted'
    pi_csv(data.frame(Status=status,Stage=stage,Version=USER_PI_VERSION,Time=as.character(Sys.time())),'run_status')
    capture.output(sessionInfo(),file=pi_path('sessionInfo.txt'))
  },add=TRUE)
  message("R version: ", R.version.string)
  message("Expression input with CK2 restored: ",MERGED_RDS)
  message("Configured output directory: ",OUT_DIR)
  message("Sole DEG-count source: ",DEG_FILE)
  message("Cell-type renaming: lateral-root-related labels -> Procambium")
  message("Augur parallel mode: ", N_THREADS, ' threads (future multisession)')
  dg<-load_user_deg()
  stage<-'Other_components'
  manifest<-compute_manifest();COMPUTE_SIGNATURE<<-digest_object(manifest)
  initialize_pi_cache(manifest)
  PI_JOBS<<-list()
  dat<-get_saved_components()
  if(is.null(dat))dat<-compute_other_components()
  cm<-dat$cells;unique_rows(cm,'cell',"Current cell metadata")
  unique_rows(unique(cm[c('pi_sample','pi_condition')]),'pi_sample',"Sample-condition mapping")
  counts<-count_user_degs(dg$data,sort(unique(as.character(cm$pi_celltype))))
  pi_csv(counts,'DEG_counts_by_celltype_condition')
  if(!any(is.finite(counts$DEG_count)))stop("DEG cell types cannot be matched to expression data; see DEG_unmatched_celltypes.")
  key<-function(ct,tr)paste(ct,tr,sep='\034')
  j<-match(key(cm$pi_celltype,cm$pi_condition),key(counts$Celltype,counts$Condition))
  nr<-dat$nondeg;unique_rows(nr,'cell',"Non-DEG components")
  ix<-match(cm$cell,nr$cell);if(anyNA(ix))stop("Non-DEG components lack cells.")
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
  writeLines(c("User-DEG four-component PI: explicit adaptation of the New Phytologist Figure 4 layout.",
               paste("Sole DEG input: ",DEG_FILE),paste("Original output directory: ",OUT_DIR),
               "Rename Lateral root and its variants to Procambium in both expression metadata and DEG tables.",
               "Count unique Gene x Celltype x Treatment entries in the user's filtered DEG table; apply no new FDR/logFC filters and do not rerun DESeq2.",
               "This restriction applies to the PI DEG-count component; Mixscape still selects features internally from expression data and reads no alternative DEG files.",
               "v3 reads Mixscape tool-slot names explicitly and validates pvec treatment/barcode alignment; other successful caches remain reusable after manifest validation.",
               "v3.1 explicitly uses future::plan(multisession) to avoid RStudio multicore fork failures producing Execution halted.",
               "Use consistent thresholds across comparisons; user filtering may differ from the original PI FDR < 0.05 threshold.",
               'PI = (minmax(UserDEG_count) + minmax(Augur AUC) + minmax(scDist) + minmax(Mixscape))/4。',
               "Scale each component jointly across available Heat/Cold/PEG cells, not separately by treatment.",
               "Constant components are set to zero under the original adaptation convention and flagged in the scaling audit.",
               "A missing component prevents a complete PI. Missing DEG records do not imply a tested zero-DEG comparison; Endodermis state contrasts are excluded.",
               "Four components are not the paper's 13 components across four resolutions; cell-type DEGs are not duplicated across resolutions.",
               "The first three components are shared within each cell type/treatment; within-group cell variation is primarily driven by Mixscape.",
               "PI represents relative expression perturbation within this dataset, not damage probability, tolerance, significance or causal evidence.",
               "Violins are descriptive distributions and white points are biological-sample means; individual cells are not independent biological replicates.",
               "Grey denotes reference background or missing components; CK is the reference and is not artificially assigned PI = 0.",
               "Retain original 05_cell_PI.csv and ensemble_PI_results.rds; new outputs use PI_UserDEG_.",
               "Legacy plotting scripts reading old PI files will not switch automatically; explicitly configure the new results.rds/cell_PI.csv.",
               "Source citation: Robertson et al., New Phytologist (2026), DOI 10.1111/nph.71378, Methods 13 / Fig.4."
  ),pi_path('README.txt'),useBytes=TRUE)
  if(!any(is.finite(cm$PI))) {
    status<-'No_complete_PI'
    failed<-bind(PI_JOBS);failed<-failed[!failed$OK,,drop=FALSE]
    if(nrow(failed))message(paste(paste(failed$Job,failed$Error,sep=': '),collapse='\n'))
    absent<-colnames(raw)[colSums(is.finite(raw))==0L]
    stop("No complete four-component PI. Entirely missing components: ",if(length(absent))paste(absent,collapse=', ') else "none (valid cells do not overlap across components)",
         ". See PI_UserDEG_jobs.csv, component_coverage.csv, missing_components.csv and Mixscape_*_tool_structure.txt.")
  }
  stage<-'Publication_plots';plot_user_pi(cm,sc$scaled,dat$umap,ss)
  status<-if(all(is.finite(cm$PI)))'Completed' else 'Completed_with_missing_PI'
  message("Completed: ",sum(is.finite(cm$PI)),'/',nrow(cm)," cells have complete PI. Figure: ",pi_path('Figure4.pdf'))
  invisible(list(cells=cm,samples=ss,output_dir=OUT_DIR))
}
if(!isTRUE(getOption('soy.userpi.skip_main',FALSE)))main()

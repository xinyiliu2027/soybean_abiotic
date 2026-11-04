# Source: Figure 4 document, extracted text lines 1-1823
# Standalone script requiring a CK2-complete RDS; a legacy module table supports optional comparisons.
# See README.md and docs/KNOWN_LIMITATIONS.md before running.

# v5: three fixed CK expression clusters named P1/P2/P3, all-cell consensus WGCNA, grey-gradient UMAPs and compact GO panels.
# All-cell consensus WGCNA with Endodermis stages P1/P2/P3 ordered along the CK reference trajectory.
# Standalone script; no earlier script is sourced and the input RDS is not overwritten.
# Define P1/P2/P3 by k-means (k = 3) on CK RNA PCA and name them using an early reference, not module scores.
MERGED_RDS <- file.path(Sys.getenv("SOYBEAN_PROJECT_ROOT", unset = "data"), "CK2_added_v3_Results/seurat_obj_with_CK2_v3.rds")
OUT_BASE <- file.path(dirname(MERGED_RDS),'Ensemble_PI_Heat_Cold_PEG_with_CK2')
# Use the explicitly configured output directory rather than a stale OUT_DIR from another session.
OUT_DIR <- OUT_BASE
# Use a distinct prefix in the same directory; reruns update this version's outputs/caches, preserving legacy Nature_ files.
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

# Recluster Endodermis and assess direction using a CK-only trajectory; project all three treatments to the same reference.
ENDO_CELLTYPE <- 'Endodermis'
EARLY_ANCHOR_TYPES <- c('Meristem','QC')
ROOT_CLUSTER <- NULL # One of the new cluster IDs 1/2/3; by default name P1 using similarity to CK Meristem/QC, not legacy cluster 5.
ROOT_CELLS_FILE <- NULL # Optional CSV with a cell column containing validated early-Endodermis barcodes; takes priority over automatic anchoring.
MIN_ROOT_CELLS <- 10L
REFERENCE_MAX_CELLS_PER_CK <- 1200L # Sample CK1/CK2 equally to fit PCA, clustering and the trajectory.
TRAJECTORY_N_HVG <- 2500L
TRAJECTORY_N_PCS <- 20L
PROJECTION_DISTANCE_QUANTILE <- 0.99
PROJECTION_DISTANCE_MULTIPLIER <- 1.5
MIN_CELLS <- 20L
MIN_METACELL_GROUP <- 40L
MAX_TRAIN_CELLS_PER_SAMPLE_CT <- 500L
N_HVG <- 4000L
USE_INPUT_HVG <- FALSE # Reselect HVGs after subdividing Endodermis.
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
SOFT_POWER_FIXED <- 6L # Retain the previous setting but recheck all eight samples; do not silently relax it on failure.
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
GENE_ID_MAP_FILE <- NULL # Columns Gene, DB_ID, Package and Keytype provide explicit cross-version mappings; do not guess IDs.
MIN_GS_SIZE <- 10L
MAX_GS_SIZE <- 500L
MIN_QUERY_GENES <- 5L
OLD_MODULE_FILE <- file.path(OUT_DIR,'EndoDev_v4_network_modules.csv') # Legacy table of 11 modules, including DevProgram02.


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
    stop("Missing R packages: ", paste(missing, collapse = ', '),
         "\nInstall prerequisites first; use BiocManager for limma/GO.db/AnnotationDbi/UCell. ",
         "Install hdWGCNA using its official instructions; prefer the matching local org.Gmax.eg.db. Missing annotation packages affect GO only.")
  }
}

pick_ct <- function(md) {
  if (!is.null(CELLTYPE_COL)) {
    if (!CELLTYPE_COL %in% names(md)) stop("CELLTYPE_COL does not exist.")
    return(CELLTYPE_COL)
  }
  hits <- names(md)[canon(names(md)) %in%
                      canon(c('cell_type','celltype','celltype_final','annotation','celltype_annotated',
                              'cell_type_annotation','seurat_annotations'))]
  if (!length(hits)) stop("Cannot identify a cell-type column; set CELLTYPE_COL.")
  if (length(hits) > 1 &&
      !all(vapply(hits[-1], function(k)
        identical(as.character(md[[k]]), as.character(md[[hits[1]]])), logical(1)))) {
    stop("Annotation columns disagree; set CELLTYPE_COL explicitly: ", paste(hits, collapse = ', '))
  }
  hits[1]
}

condition_key <- function(x) {
  aliases <- list(
    CK   = c('CK','Control','CTRL','REF','Reference',"\u5bf9\u7167"),
    Heat = c('Heat','High_Temp','HighTemp','HighTemperature','HT',"\u9ad8\u6e29"),
    Cold = c('Cold','Low_Temp','LowTemp','LowTemperature','LT',"\u4f4e\u6e29"),
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
      note("Reusing this version's cache: ", name)
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
  if (!all(required %in% names(f))) stop("Soft-threshold diagnostics lack required columns.")
  if (length(SOFT_POWER_FIXED)!=1L || !is.finite(SOFT_POWER_FIXED) || SOFT_POWER_FIXED<=0)
    stop("SOFT_POWER_FIXED must be a single positive number.")
  f$Chosen <- !is.na(f$Power) & f$Power==SOFT_POWER_FIXED
  if (sum(f$Chosen)!=1L) stop("No unique row for the fixed power in the diagnostic table: ",SOFT_POWER_FIXED)
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

# Sample aliases must identify specific replicates; REF/Control indicate treatment and cannot alone identify CK1.
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
  if (!is.null(SAMPLE_COL) && !SAMPLE_COL %in% names(md)) stop("Specified SAMPLE_COL does not exist.")
  cols <- intersect(unique(c(SAMPLE_COL, 'sample_id_complete', 'sample_id', 'sample', 'orig.ident')), names(md))
  if (!length(cols)) stop("No identifiable sample column.")
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
  if (any(conflict)) stop("Metadata columns disagree on replicate identity; see input_sample_identity_audit.csv.")
  
  # Infer condition from independent-sample identity, then cross-check existing treatment columns.
  design <- data.frame(
    Sample = c('CK1','CK2','highTemp1','highTemp2','lowTemp1','lowTemp2','PEG1','PEG2'),
    Condition = rep(c('CK','Heat','Cold','PEG'), each = 2),
    Replicate = rep(1:2, 4), stringsAsFactors = FALSE)
  cond <- design$Condition[match(chosen, design$Sample)]
  if (!is.null(CONDITION_COL) && !CONDITION_COL %in% names(md)) stop("Specified CONDITION_COL does not exist.")
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
  if (any(condition_conflict)) stop("Sample identity conflicts with treatment columns; see input_sample_identity_audit.csv.")
  if (any(audit$Unresolved_target_sample)) {
    stop("Some target-condition cells lack distinguishable replicates. REF/Control cannot replace CK1/CK2; inspect the sample-identity audit.")
  }
  absent <- setdiff(design$Sample, unique(chosen))
  if (length(absent)) stop("Missing target biological samples: ", paste(absent, collapse = ', '),
                           ". Confirm the input contains CK2; do not copy or split CK1.")
  list(Sample = chosen, Condition = cond, design = design, audit = audit)
}

# JoinLayers only on an assay copy; avoid subsetting or rewriting the loaded legacy Seurat object and its reduction ordering.
extract_counts_readonly <- function(obj, cells, assay = ASSAY) {
  if (!assay %in% names(obj@assays)) stop("Missing assay: ", assay)
  a <- obj@assays[[assay]]
  if (inherits(a, 'Assay5')) {
    layer_names <- SeuratObject::Layers(a, search = '^counts')
    if (!length(layer_names)) stop("Raw RNA counts are missing.")
    a <- SeuratObject::JoinLayers(a, layers = 'counts', new = 'counts')
    counts <- SeuratObject::LayerData(a, layer = 'counts')
  } else counts <- SeuratObject::GetAssayData(a, slot = 'counts')
  if (is.null(rownames(counts)) || is.null(colnames(counts)) ||
      anyDuplicated(rownames(counts)) || anyDuplicated(colnames(counts))) stop("Counts gene/cell names are missing or duplicated.")
  if (!all(cells %in% colnames(counts))) stop("Counts lack selected cells.")
  counts <- methods::as(counts[, cells, drop = FALSE], 'dgCMatrix')
  if (any(!is.finite(counts@x)) || any(counts@x < 0) ||
      any(abs(counts@x - round(counts@x)) > 1e-6)) stop("RNA counts must be raw nonnegative integers.")
  counts
}

extract_embedding_readonly <- function(obj, reduction, cells) {
  if (!reduction %in% names(obj@reductions)) stop("Missing reduction: ", reduction)
  emb <- obj@reductions[[reduction]]@cell.embeddings
  if (is.null(rownames(emb)) || anyDuplicated(rownames(emb)) ||
      !all(cells %in% rownames(emb))) stop(reduction, " does not cover selected cells or has duplicated cell names.")
  emb <- emb[cells, , drop = FALSE]
  if (any(!is.finite(emb))) stop(reduction, " contains non-finite values.")
  emb
}


read_data <- function() {
  obj <- readRDS(MERGED_RDS)
  if (!inherits(obj, 'Seurat')) stop("Input is not a Seurat object.")
  md <- obj@meta.data
  if (is.null(rownames(md)) || anyDuplicated(rownames(md))) stop("Metadata cell names are missing or duplicated.")
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
  if (!nrow(m) || anyNA(m$Sample)) stop("No usable target cells or sample IDs are missing.")
  pairs <- unique(m[c('Sample','Condition')])
  if (anyDuplicated(pairs$Sample)) stop("A sample spans multiple conditions.")
  if (!setequal(pairs$Sample, ids$design$Sample)) stop("A target sample was lost after filtering missing annotations.")
  ns <- table(factor(pairs$Condition, CONDITIONS))
  print(ns)
  if (any(ns != EXPECTED_REPLICATES)) stop("Each condition requires two independent biological replicates.")
  if (any(grepl('^[0-9]+$', m$Celltype))) stop("The cell-type column contains numeric clusters; use cell_type.")
  pairs <- pairs[match(ids$design$Sample, pairs$Sample), , drop = FALSE]
  pairs$Replicate <- ids$design$Replicate
  pairs$Cells <- as.integer(table(m$Sample)[pairs$Sample])
  csv(pairs, 'biological_samples.csv')
  print(pairs, row.names = FALSE)
  csv(data.frame(Input = MERGED_RDS, Total_input_cells = nrow(md),
                 Identified_input_samples = length(unique(na.omit(ids$Sample))),
                 Analysis_cells = nrow(m), Analysis_samples = nrow(pairs)),
      'input_object_summary.csv')
  note("Input cells: ", nrow(md), "; retaining only CK/Heat/Cold/PEG: ", nrow(m), " cells across eight independent samples.")
  
  counts <- extract_counts_readonly(obj, m$cell)
  m$nUMI <- as.numeric(Matrix::colSums(counts))
  m$nFeature <- as.numeric(Matrix::colSums(counts > 0))
  if (any(m$nUMI <= 0)) stop("Zero-library cells detected.")
  xy <- extract_embedding_readonly(obj, UMAP_REDUCTION, m$cell)
  if (ncol(xy) < 2L) stop("Specified UMAP has fewer than two dimensions.")
  m$UMAP1 <- xy[, 1]; m$UMAP2 <- xy[, 2]
  rownames(m) <- m$cell
  if (!is.null(HARMONY_COL)) {
    if (!HARMONY_COL %in% names(md)) stop("HARMONY_COL does not exist.")
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
  # If a legacy reduction does not cover added CK2 cells, automatic mode uses RNA PCA; UMAP is not used for network construction.
  complete <- vapply(hit, function(nm) all(m$cell %in% rownames(obj@reductions[[nm]]@cell.embeddings)), logical(1))
  hit <- hit[complete]
  if (!is.null(GRAPH_REDUCTION) && !length(hit)) stop("Specified GRAPH_REDUCTION is missing or does not cover all target cells.")
  graph <- if (length(hit)) extract_embedding_readonly(obj, hit[1], m$cell) else NULL
  list(counts = counts, md = m, qc = qc,
       input_hvg = SeuratObject::VariableFeatures(obj@assays[[ASSAY]]),
       graph = graph, graph_name = if (length(hit)) hit[1] else 'RNA_PCA_new')
}


# -----------------------------------------------------------------------------
# 6. Build the consensus network.
# -----------------------------------------------------------------------------
# Validate metacell expression matrices; averaged metacells may contain fractional values and should not undergo raw-count integer checks.
validate_metacell_matrix <- function(x, cells, tag) {
  if (length(dim(x)) != 2L || !nrow(x) || !ncol(x)) stop(tag, ": empty metacell matrix.")
  if (is.null(rownames(x)) || is.null(colnames(x)) ||
      anyNA(rownames(x)) || anyNA(colnames(x)) ||
      anyDuplicated(rownames(x)) || anyDuplicated(colnames(x))) stop(tag, ": matrix row/column names are missing or duplicated.")
  if (!setequal(colnames(x), cells)) stop(tag, ": expression and metadata metacell sets disagree.")
  x <- methods::as(x[, cells, drop = FALSE], 'dgCMatrix')
  if (any(!is.finite(x@x)) || any(x@x < 0)) stop(tag, ": expression contains non-finite or negative values.")
  if (any(Matrix::colSums(x) <= 0)) stop(tag, ": zero-expression metacells detected.")
  x
}

# Extract hdWGCNA's internal metacells, join RNA counts layers, and convert to a classic Assay.
# Preserve metacell metadata, misc and grouping; do not regenerate metacells or average/round the aggregated matrix again.
prepare_metacell_rna <- function(train, nm, expected_samples) {
  mc <- hdWGCNA::GetMetacellObject(train, wgcna_name = nm)
  if (!inherits(mc, 'Seurat')) stop("hdWGCNA did not return a valid metacell Seurat object.")
  if (!'RNA' %in% names(mc@assays)) stop("Metacell object has no RNA assay.")
  md_before <- mc[[]]
  cells <- rownames(md_before)
  needed <- c('ers_sample','ers_celltype','ers_pool')
  if (!all(needed %in% names(md_before))) stop("Metacell object lost sample or cell-type grouping.")
  if (!setequal(as.character(md_before$ers_sample), expected_samples)) stop("Metacell object does not contain all eight samples.")
  a <- mc@assays[['RNA']]
  class_before <- class(a)[1]
  layers_before <- if (inherits(a, 'Assay5')) SeuratObject::Layers(a) else c('counts','data','scale.data')
  csv(data.frame(Stage = 'Before_conversion', RNA_class = class_before, Layer = layers_before),
      'metacell_layers_before.csv')
  if (inherits(a, 'Assay5')) {
    count_layers <- SeuratObject::Layers(a, search = '^counts')
    if (!length(count_layers)) stop("Metacell object lacks a counts layer.")
    # Modify only the isolated assay; discard legacy data.* layers and log-normalize consistently in the next step.
    a <- SeuratObject::JoinLayers(a, layers = 'counts', new = 'counts')
    counts <- SeuratObject::LayerData(a, layer = 'counts')
  } else counts <- SeuratObject::GetAssayData(a, slot = 'counts')
  counts <- validate_metacell_matrix(counts, cells, 'metacell counts')
  classic <- SeuratObject::CreateAssayObject(counts = counts, min.cells = 0, min.features = 0)
  mc[['RNA']] <- classic
  SeuratObject::DefaultAssay(mc) <- 'RNA'
  if (inherits(mc[['RNA']], 'Assay5')) stop("Metacell RNA was not converted to a classic Assay.")
  actual <- SeuratObject::GetAssayData(mc[['RNA']], slot = 'counts')
  if (!identical(rownames(actual), rownames(counts)) ||
      !setequal(colnames(actual), colnames(counts)) ||
      Matrix::nnzero(actual[, cells, drop = FALSE] - counts) != 0) {
    stop("Expression values or gene names changed during metacell conversion.")
  }
  keep_cols <- intersect(c(needed, 'cells_merged'), names(md_before))
  md_after <- mc[[]]
  for (co in keep_cols) {
    if (!identical(as.character(md_before[[co]]), as.character(md_after[cells, co]))) {
      stop("Grouping/membership changed during metacell conversion: ", co)
    }
  }
  methods::validObject(mc)
  csv(data.frame(RNA_class_before = class_before, RNA_class_after = class(mc[['RNA']])[1],
                 Metacells = ncol(counts), Genes = nrow(counts),
                 Values_unchanged = TRUE,
                 Fractional_values_preserved = any(abs(counts@x - round(counts@x)) > 1e-6)),
      'metacell_RNA_conversion.csv')
  train <- hdWGCNA::SetMetacellObject(train, metacell_obj = mc, wgcna_name = nm)
  note("Internal metacell RNA converted to a classic Assay; preserving all ", ncol(counts), " metacells and their expression values.")
  train
}

# Before SetMultiExpr, test its per-sample GetAssayData access to detect missing layers or samples early.
check_metacell_expression <- function(train, nm, expected_samples) {
  mc <- hdWGCNA::GetMetacellObject(train, wgcna_name = nm)
  if (inherits(mc[['RNA']], 'Assay5')) stop("Normalized metacell RNA unexpectedly became Assay5 again.")
  md <- mc[[]]
  checks <- lapply(expected_samples, function(s) {
    cells <- rownames(md)[as.character(md$ers_sample) == s]
    if (!length(cells)) stop("Metacells lack samples: ", s)
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
  note("Building a cross-sample consensus network; the first run may take time.")
  m <- dat$md
  
  # 6.1 Build a balanced training set by sampling within sample-by-cell-type groups.
  groups <- split(seq_len(nrow(m)), interaction(m$Sample, m$Celltype, drop = TRUE))
  set.seed(SEED)
  ix <- unlist(lapply(groups, function(i) {
    if (length(i) < MIN_METACELL_GROUP || any(!m$Network_eligible[i])) integer()
    else sample(i, min(length(i), MAX_TRAIN_CELLS_PER_SAMPLE_CT))
  }), use.names = FALSE)
  if (!length(ix)) stop("Too few cells to create metacells.")
  
  mt <- m[ix, ]
  if (!setequal(mt$Sample, m$Sample)) stop("Some samples have no cell types eligible for network construction.")
  
  train <- SeuratObject::CreateSeuratObject(
    counts = SeuratObject::CreateAssayObject(
      counts = dat$counts[, mt$cell, drop = FALSE]),
    meta.data = mt)
  train$ers_celltype <- mt$Celltype
  train$ers_sample   <- mt$Sample
  train$ers_pool     <- 'All'
  train <- Seurat::NormalizeData(train, verbose = FALSE)
  train <- Seurat::FindVariableFeatures(train, nfeatures = N_HVG, verbose = FALSE)
  
  # 6.2 Retain genes meeting the detection threshold in at least one cell type (the original OR rule).
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
  if (length(genes) < 200) stop("Fewer than 200 network genes; check expression and filtering.")
  csv(data.frame(Gene = genes), 'network_input_genes.csv')
  
  # 6.3 Use PCA or an external reduction.
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
  # Prepare the internal metacell RNA assay; preparing train[['RNA']] alone does not fix the layer-access issue.
  train <- prepare_metacell_rna(train, nm, sort(unique(mt$Sample)))
  train <- hdWGCNA::NormalizeMetacells(train, wgcna_name = nm, assay = 'RNA', verbose = FALSE)
  check_metacell_expression(train, nm, sort(unique(mt$Sample)))
  
  # 6.5 Construct multiExpr with one dataset per sample.
  train <- hdWGCNA::SetMultiExpr(
    train,
    group_name = 'All', group.by = 'ers_pool',
    multi.group.by = 'ers_sample',
    multi_groups = sort(unique(mt$Sample)),
    assay = 'RNA', slot = 'data', wgcna_name = nm)
  
  multi <- hdWGCNA::GetMultiExpr(train, wgcna_name = nm)
  if (length(multi) != length(unique(mt$Sample))) stop("Consensus expression matrices lack samples.")
  if (is.null(names(multi))) stop("multiExpr lacks sample names; order cannot be audited.")
  if (!setequal(names(multi), unique(mt$Sample)))
    stop("multiExpr sample names disagree; check the hdWGCNA version.")
  
  bg <- colnames(multi[[1]]$data)
  if (any(vapply(multi, function(x) !identical(colnames(x$data), bg), logical(1))))
    stop("Network gene order differs across samples.")
  
  qc <- data.frame(
    Sample    = names(multi),
    Metacells = vapply(multi, function(x) nrow(x$data), integer(1)),
    Genes     = length(bg))
  csv(qc, 'metacell_QC.csv')
  if (any(qc$Metacells < MIN_METACELLS_PER_SAMPLE))
    stop("A sample has too few metacells; see ERS_metacell_QC.csv.")
  
  # 6.6 Run per-sample diagnostics with shared fixed power = 6.
  fits <- list(); powers <- numeric(length(multi))
  for (i in seq_along(multi)) {
    x <- multi[[i]]$data
    if (any(!is.finite(as.matrix(x))) || any(apply(x, 2, stats::sd) == 0))
      stop("Consensus matrices contain non-finite values or zero-variance genes; check SetMultiExpr filtering.")
    
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
    stop("Fixed power=",SOFT_POWER_FIXED," failed R2/connectivity checks: ",
         paste(selected$Sample[!selected$QC_pass],collapse=', '),
         "; inspect Nature_soft_power_selected.csv. Thresholds were not changed automatically.")
  note("Shared power across samples=",SOFT_POWER_FIXED,"; R2 and mean-connectivity checks passed.")
  
  # 6.7 Construct the consensus network and modules.
  train <- hdWGCNA::ConstructNetwork(
    train, consensus = TRUE, soft_power = unname(powers),
    tom_outdir = out('TOM'), tom_name = 'ERS_consensus',
    overwrite_tom = TRUE, wgcna_name = nm,
    networkType = 'signed', TOMType = 'signed',
    minModuleSize = MIN_MODULE_SIZE, deepSplit = 2,
    maxBlockSize = max(5000L, length(bg) + 1L), randomSeed = SEED)
  
  mods <- hdWGCNA::GetModules(train, wgcna_name = nm)
  if (!all(c('gene_name','module') %in% names(mods)))
    stop("Module table lacks gene_name/module.")
  
  ids <- sort(unique(as.character(mods$module[
    !tolower(as.character(mods$module)) %in% c('grey','gray','0') &
      !is.na(mods$module)])))
  if (!length(ids)) stop("No non-grey modules; shared programs cannot be forced.")
  
  map <- setNames(sprintf('StageProgram%02d', seq_along(ids)), ids)
  mods$original_module <- as.character(mods$module)
  mods$Program <- unname(map[as.character(mods$module)])
  mods$module  <- factor(ifelse(is.na(mods$Program), 'grey', mods$Program),
                         levels = c(unname(map), 'grey'))
  train <- hdWGCNA::SetModules(train, mods, wgcna_name = nm)
  if (anyDuplicated(mods$gene_name)) stop("Network module table contains duplicated genes.")
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
    stop("UCell cell names do not match; actual dimensions=", paste(dim(x), collapse='x'),
         "; first row names=", paste(head(rownames(x)), collapse=', '),
         "; expected cells=", paste(head(cells), collapse=', '))
  }
  expected <- paste0(programs, '_UCell')
  if (valid_axis(colnames(x), expected)) {
    x <- x[cells, expected, drop=FALSE]
  } else if (valid_axis(colnames(x), programs)) {
    x <- x[cells, programs, drop=FALSE]
  } else {
    stop("UCell module names do not match; actual=", paste(colnames(x), collapse=', '),
         "; expected=", paste(expected, collapse=', '))
  }
  if (!is.numeric(x) || any(!is.finite(x))) stop("UCell contains nonnumeric or non-finite values.")
  colnames(x) <- programs
  x
}

score_modules <- function(dat, net) {
  z <- net$modules[!is.na(net$modules$Program), , drop=FALSE]
  sets <- lapply(split(as.character(z$gene_name), as.character(z$Program)), unique)
  if (!length(sets) || any(lengths(sets) == 0L)) stop("No modules available for scoring.")
  counts <- dat$counts
  if (is.null(colnames(counts)) || anyDuplicated(colnames(counts)))
    stop("RNA matrix lacks unique cell names.")
  if (any(!unlist(sets, use.names=FALSE) %in% rownames(counts)))
    stop("Module genes are absent from the RNA matrix.")
  chunks <- split(seq_len(ncol(counts)), ceiling(seq_len(ncol(counts))/SCORE_CHUNK))
  ans <- matrix(NA_real_, nrow=ncol(counts), ncol=length(sets),
                dimnames=list(colnames(counts), names(sets)))
  for (i in seq_along(chunks)) {
    ix <- chunks[[i]]
    note("All-cell UCell: ", i, '/', length(chunks))
    sc <- UCell::ScoreSignatures_UCell(
      matrix=counts[, ix, drop=FALSE], features=sets,
      maxRank=min(max(UCell_MAX_RANK, max(lengths(sets))+1L), nrow(counts)),
      name='_UCell', ncores=1)
    if (i == 1L) {
      note("UCell output class: ", paste(class(sc), collapse='/'),
           "; dimensions: ", paste(dim(sc), collapse=' x '))
    }
    ans[ix, ] <- align_ucell_scores(sc, colnames(counts)[ix], names(sets))
  }
  if (any(!is.finite(ans))) stop("UCell contains non-finite values.")
  ans
}


aggregate_scores <- function(scores, md) {
  if (anyDuplicated(md$cell) || anyDuplicated(rownames(scores)) ||
      !all(md$cell %in% rownames(scores))) stop("Score aggregation: cell IDs disagree.")
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
# 8. Sample-level tests and shared-program summaries.
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
    stop(label, " has mismatched cell/module names.")
  x <- x[cells, programs, drop = FALSE]
  if (any(!is.finite(x))) stop(label, " contains non-finite values.")
  x
}

all_gene_seurat_scores <- function(x) {
  mods  <- hdWGCNA::GetModules(x, wgcna_name = NETWORK_NAME)
  progs <- levels(mods$module)
  progs <- progs[progs != 'grey']
  if (!length(progs)) stop("No non-grey modules.")
  set.seed(SEED)
  
  if (length(progs) > 1) {
    x <- hdWGCNA::ModuleExprScore(
      x, n_genes = 'all', method = 'Seurat', wgcna_name = NETWORK_NAME,
      assay = 'RNA', nbin = SEURAT_NBIN, ctrl = SEURAT_CTRL,
      seed = SEED, name = 'NatureScore')
  } else {
    # Handle single-module output explicitly because drop = TRUE in some hdWGCNA versions removes column names.
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
  
  # Average log-normalized RNA explicitly to avoid AvgModuleExpr version differences.
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
  # Determine each module's color limits using all cells in the score matrix; share them across treatments and Endodermis zooms.
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
  if(!all(c('Program','gene_name') %in% names(old))) stop("Legacy module table lacks Program or gene_name.")
  genes <- unique(as.character(old$gene_name[old$Program=='DevProgram02' & !is.na(old$Program)]))
  genes <- intersect(genes,rownames(x))
  if(length(genes)<5) {
    csv(data.frame(Status='Insufficient_previous_DevProgram02_genes',Genes=length(genes)),
        'previous_Program02_status.csv'); return(invisible(NULL))
  }
  # Keep the legacy Program02 gene set fixed and evaluate it independently on new stages; never redefine P1 using scores.
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
    stop("A genuine cross-treatment batch column HARMONY_COL is required.")
  s <- unique(md[, c('Sample','Condition','Batch')])
  if (anyNA(s$Batch) || anyDuplicated(s$Sample) || length(unique(s$Batch)) < 2)
    stop("Batch is missing, not sample-level, or has only one level.")
  mm <- model.matrix(~ Condition + Batch, s)
  if (qr(mm)$rank < ncol(mm)) stop("Batch is confounded with treatment; Harmony will not run.")
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
        if (!co %in% names(modules)) stop("Missing kME columns: ", co)
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
    note("ME/Hub step incomplete: ", conditionMessage(e))
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
  # Remove transcript suffixes only from standard Glyma loci, preserving genome-version identifiers such as v2/v4.
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
  if (length(objects)!=1L) stop("Cannot uniquely identify an OrgDb object in package: ",pkg)
  objects[[1]]
}
bp_annotation <- function(db,ids,keytype) {
  cols <- AnnotationDbi::columns(db)
  if (all(c('GOALL','ONTOLOGYALL') %in% cols)) {
    a <- suppressMessages(AnnotationDbi::select(db,keys=ids,keytype=keytype,columns=c('GOALL','ONTOLOGYALL')))
    a <- a[a$ONTOLOGYALL %in% 'BP' & !is.na(a$GOALL),]
    return(unique(data.frame(DB_ID=a[[keytype]],GO=a$GOALL)))
  }
  if (!'GO' %in% cols) stop("OrgDb lacks GO and GOALL columns.")
  a <- suppressMessages(AnnotationDbi::select(db,keys=ids,keytype=keytype,columns='GO'))
  a <- unique(data.frame(DB_ID=a[[keytype]],GO=a$GO))
  a <- a[!is.na(a$GO),]
  bpkeys <- AnnotationDbi::mappedkeys(GO.db::GOBPANCESTOR)
  a <- a[a$GO %in% bpkeys,]
  if (!nrow(a)) return(a)
  # If a custom OrgDb lacks GOALL, add BP ancestors explicitly for a comparable enrichment definition.
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
      stop("Mapping table requires nonmissing Gene/DB_ID/Package/Keytype.")
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
  if (!length(candidates)) stop("Neither OrgDb provides usable matches; inspect OrgDb_mapping_diagnostics.csv. Mapping failure will not be reported as absent enrichment.")
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


# Endodermis reference developmental trajectory.
balance_cells <- function(md,types,max_n) {
  ix <- which(md$Condition=='CK' & md$Celltype %in% types)
  g <- split(ix,md$Sample[ix])
  if (!setequal(names(g),c('CK1','CK2'))) stop("Specified reference cells are incomplete across CK1/CK2.")
  n <- min(max_n,lengths(g)); if(n<MIN_CELLS) stop("Too few CK reference cells.")
  set.seed(SEED); unlist(lapply(g,sample,size=n),use.names=FALSE)
}
cluster_means <- function(pca,clusters) {
  z <- t(vapply(split(seq_len(nrow(pca)),clusters),function(i)
    colMeans(pca[i,,drop=FALSE]),numeric(ncol(pca))))
  colnames(z) <- colnames(pca); z
}
fit_reference_pca <- function(x,npcs,svd_engine=NULL) {
  # x must contain CK reference cells by genes; record variable order explicitly instead of relying on SVD dimnames.
  x <- as.matrix(x)
  if(!is.numeric(x) || any(!is.finite(x))) stop("CK PCA training matrix is nonnumeric or contains non-finite values.")
  genes <- colnames(x)
  if(is.null(genes) || anyNA(genes) || any(!nzchar(genes)) || anyDuplicated(genes))
    stop("CK PCA training gene names are missing or duplicated.")
  if(nrow(x)<5L) stop("Insufficient CK PCA reference cells.")
  mu <- colMeans(x); sigma <- apply(x,2L,stats::sd)
  keep <- is.finite(sigma) & sigma>1e-8
  audit <- data.frame(Gene=genes,Center=as.numeric(mu),Scale=as.numeric(sigma),
                      Retained=keep,Reason=ifelse(keep,'Retained','Zero_or_near_zero_variance'))
  x <- x[,keep,drop=FALSE]; genes <- colnames(x)
  mu <- mu[keep]; sigma <- sigma[keep]
  np <- min(as.integer(npcs),nrow(x)-2L,ncol(x)-2L)
  if(!is.finite(np) || np<3L) stop("Insufficient CK reference PCA dimensions after removing zero-variance genes.")
  z <- sweep(sweep(x,2L,mu,'-'),2L,sigma,'/')
  if(any(!is.finite(z))) stop("Standardizing the CK training matrix produced non-finite values.")
  # Apply truncated SVD to an explicitly standardized matrix, not prcomp_irlba center/scale fields.
  if(is.null(svd_engine)) svd_engine <- function(z,n)
    irlba::irlba(z,nv=n,nu=0L,maxit=1000L,tol=1e-6)
  set.seed(SEED)
  s <- svd_engine(z,np)
  rotation <- as.matrix(s$v)
  if(nrow(rotation)!=length(genes) || ncol(rotation)!=np)
    stop("Unexpected SVD loading dimensions: expected ",length(genes),'×',np,"; actual ",paste(dim(rotation),collapse='×'))
  if(any(!is.finite(rotation)) || length(s$d)!=np || any(!is.finite(s$d)))
    stop("SVD returned invalid loadings or singular values.")
  # SVD V rows follow input-variable order; attach gene names explicitly after checking dimensions.
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
    stop("PCA training genes disagree with loading row names.")
  if(!is.numeric(fit$center) || !is.numeric(fit$scale) ||
     length(fit$center)!=length(genes) || length(fit$scale)!=length(genes) ||
     !identical(names(fit$center),genes) || !identical(names(fit$scale),genes) ||
     any(!is.finite(fit$center)) || any(!is.finite(fit$scale)) || any(fit$scale<=0))
    stop("PCA centers/standard deviations must be finite vectors ordered by training genes.")
  if(is.null(rownames(expr)) || anyDuplicated(rownames(expr)) ||
     !all(genes %in% rownames(expr))) stop("Projection RNA matrix lacks training genes or has duplicated gene names.")
  if(is.null(colnames(expr)) || anyDuplicated(colnames(expr)) ||
     anyNA(allcells) || anyDuplicated(allcells) || !all(allcells %in% colnames(expr)))
    stop("Projection RNA cell names are missing, duplicated or mismatched.")
  scores <- matrix(NA_real_,length(allcells),ncol(rotation),
                   dimnames=list(allcells,colnames(rotation)))
  chunks <- split(seq_along(allcells),ceiling(seq_along(allcells)/chunk_size))
  for (ix in chunks) {
    y <- t(as.matrix(expr[genes,allcells[ix],drop=FALSE]))
    if(!identical(colnames(y),genes) || ncol(y)!=nrow(rotation))
      stop("Incorrect PCA projection dimensions/gene order: query=",paste(dim(y),collapse='×'),
           '；loadings=',paste(dim(rotation),collapse='×'))
    if(any(!is.finite(y))) stop("Expression matrix to project contains non-finite values.")
    y <- sweep(sweep(y,2L,fit$center,'-'),2L,fit$scale,'/')
    scores[ix,] <- y %*% rotation
  }
  if(any(!is.finite(scores))) stop("PCA projection produced non-finite results.")
  scores
}
reference_pca <- function(expr,genes,refcells,allcells) {
  if(is.null(rownames(expr)) || is.null(colnames(expr)) ||
     anyDuplicated(rownames(expr)) || anyDuplicated(colnames(expr)))
    stop("RNA expression matrix requires unique gene and cell names.")
  if(anyNA(refcells) || anyDuplicated(refcells) || !all(refcells %in% allcells))
    stop("CK reference cells are duplicated, missing or outside the projected cell set.")
  genes <- unique(as.character(genes)); genes <- genes[!is.na(genes) & nzchar(genes)]
  absent <- setdiff(genes,rownames(expr))
  genes <- genes[genes %in% rownames(expr)]
  if(length(genes)<5L || !all(refcells %in% colnames(expr))) stop("CK PCA genes or cells do not match.")
  x <- t(as.matrix(expr[genes,refcells,drop=FALSE]))
  trained <- fit_reference_pca(x,TRAJECTORY_N_PCS)
  fit <- trained$fit
  scores <- project_reference_pca(expr,allcells,fit)
  # Blockwise projection of the CK matrix must reproduce its training scores.
  error <- max(abs(scores[refcells,,drop=FALSE]-fit$x))
  if(!is.finite(error) || error>1e-8*max(1,max(abs(fit$x))))
    stop("CK training scores differ from reprojected scores.")
  csv(trained$audit,'trajectory_PCA_gene_scaling.csv')
  csv(data.frame(Gene=absent),'trajectory_PCA_absent_HVG.csv')
  csv(data.frame(Reference_cells=nrow(x),Requested_present_genes=ncol(x),
                 Retained_genes=nrow(fit$rotation),Returned_PCs=ncol(fit$rotation),
                 Projection_cells=nrow(scores),Reference_projection_max_error=error,
                 Method='Explicit CK centering/scaling + truncated SVD; named feature alignment'),
      'trajectory_PCA_projection_QC.csv')
  note("CK PCA validation passed: ",nrow(fit$rotation)," genes, ",ncol(fit$rotation),
       " PCs; projecting ",nrow(scores)," cells.")
  list(fit=fit,scores=scores)
}
# Fix three CK expression clusters; P1/P2/P3 are not pseudotime tertiles and are not defined by module scores.
nearest_centroid <- function(x,centers) {
  x <- as.matrix(x); centers <- as.matrix(centers)
  if(ncol(x)!=ncol(centers) || any(!is.finite(x)) || any(!is.finite(centers)))
    stop("Nearest-centroid classification: invalid dimensions or values.")
  d2 <- vapply(seq_len(nrow(centers)),function(j)
    rowSums(sweep(x,2,centers[j,],'-')^2),numeric(nrow(x)))
  if(is.null(dim(d2))) d2 <- matrix(d2,nrow=nrow(x))
  best <- max.col(-d2,ties.method='first')
  ord <- t(apply(d2,1,sort))
  list(cluster=rownames(centers)[best],distance=sqrt(d2[cbind(seq_len(nrow(x)),best)]),
       margin=(ord[,2]-ord[,1])/pmax(ord[,2],1e-12))
}
project_fixed_path <- function(x,path) {
  # Project onto fixed arc length along the reference curve; do not reset the origin for each query sample or block.
  x <- as.matrix(x); path <- as.matrix(path)
  if(ncol(x)!=ncol(path) || any(!is.finite(x)) || any(!is.finite(path))) stop("Invalid trajectory-projection input.")
  seg <- path[-1,,drop=FALSE]-path[-nrow(path),,drop=FALSE]
  lens <- sqrt(rowSums(seg^2)); good <- which(lens>1e-10)
  if(!length(good)) stop("Reference curve has zero length.")
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
  if(length(unique(fit$cluster))!=3L) stop("Fixed three-cluster analysis did not produce three nonempty groups.")
  distance <- vapply(seq_len(3),function(j)
    mean(sqrt(rowSums(sweep(anchor_centers,2,centers[j,],'-')^2))),numeric(1))
  root <- ROOT_CLUSTER; basis <- 'CK_Meristem_QC_expression_proximity'
  if(is.null(root) && length(root_cells)) {
    hit <- intersect(root_cells,rownames(rpc))
    if(length(hit)<MIN_ROOT_CELLS) stop("Too few specified early-cell barcodes match the CK reference.")
    tab <- sort(table(fit$cluster[match(hit,rownames(rpc))]),decreasing=TRUE)
    if(as.numeric(tab[1])/sum(tab)<=0.5) stop("Early-cell barcodes do not uniquely support one cluster.")
    root <- names(tab)[1]; basis <- 'Validated_early_barcodes'
  }
  if(is.null(root)) root <- rownames(centers)[which.min(distance)]
  else if(is.null(root_cells)) basis <- 'Explicit_root_cluster'
  root <- as.character(root)
  if(length(root)!=1L || !root %in% rownames(centers)) stop("ROOT_CLUSTER must be 1, 2 or 3 in the new clustering, not legacy cluster 5.")
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
  if(!length(endo) || !all(c('CK1','CK2') %in% m$Sample[anchors])) stop("Endodermis or early CK reference cells are missing.")
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
    if(!'cell' %in% names(root_table)) stop("ROOT_CELLS_FILE must contain a cell column.")
    roots <- as.character(root_table$cell)
  }
  cl <- fit_three_clusters(rpc,ac,roots)
  csv(cl$audit,'trajectory_root_and_three_cluster_order.csv')
  csv(as.data.frame(table(Sample=m$Sample[refix],Period=cl$map[as.character(cl$fit$cluster)])),
      'CK_three_cluster_sample_support.csv')
  note("Three fixed CK expression clusters: ",paste(paste(names(cl$map),cl$map,sep='='),collapse=', '),
       ". Stages are not defined by module scores or equal pseudotime intervals.")
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
  # The continuous trajectory is independent of module scores; record unreliable fits and use a three-centroid polyline approximation.
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
    note("Approximating the continuous trajectory with a three-centroid polyline; reason: ",reason,". Expression clusters remain unchanged; approximate coordinates are not validated developmental time.")
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
  # Cluster markers use descriptive mean differences; cell counts are not biological replicate counts.
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
  if(anyNA(j) || anyDuplicated(z$cell)) stop("Invalid cell names in the three-cluster result.")
  for(k in c('Period','Pseudotime','Pseudotime01','Endo_branch','Stage_status','Cluster_margin')) {
    m[[k]] <- z[[k]][rep(NA_integer_,nrow(m))]; m[[k]][ix] <- z[[k]][j]
  }
  m$Celltype[ix] <- paste0('Endodermis_',m$Period[ix])
  m$Network_eligible <- TRUE; m$Network_eligible[ix] <- z$Projection_pass[j]
  m$Analysis_eligible <- m$Network_eligible
  if(!setequal(unique(m$Period[ix]),c('P1','P2','P3'))) stop("Not all three Endodermis clusters were obtained.")
  csv(m,'all_cells_with_Endodermis_stages.csv')
  csv(as.data.frame(table(Sample=m$Sample,Celltype=m$Celltype,Eligible=m$Analysis_eligible)),
      'sample_stage_network_coverage.csv')
  dat$md <- m; dat$input_hvg <- character(); dat['graph'] <- list(NULL)
  dat$graph_name <- 'RNA_PCA_after_three_Endodermis_clusters'; dat
}

aggregate_lineage_scores <- function(scores,long) {
  if(anyDuplicated(paste(long$cell,long$Branch,sep='::'))) stop("Duplicated cell-by-branch entries in statistical input.")
  z <- long[long$Projection_pass & !is.na(long$Period),,drop=FALSE]
  empty <- data.frame(Branch=character(),Lineage=character(),Sample=character(),
                      Condition=character(),Period=character(),n_cells=integer(),Program=character(),Score=numeric(),
                      mean_logUMI=numeric(),mean_nFeature=numeric())
  if(!nrow(z)) return(empty)
  if(!all(z$cell %in% rownames(scores))) stop("Branch cells lack module scores.")
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
  # Adjust jointly across tested module-by-branch-by-treatment combinations; shared early cells make branches unsuitable as independent validation.
  if(nrow(tests)) tests$FDR_global <- p.adjust(tests$P_value,'BH')
  else tests <- data.frame(Branch=character(),Program=character(),Treatment=character(),
                           Delta_P2_vs_P1=numeric(),Delta_P3_vs_P1=numeric(),F_statistic=numeric(),P_value=numeric(),FDR_global=numeric())
  if(nrow(bases)) bases$FDR_global <- p.adjust(bases$P_value,'BH')
  else bases <- data.frame(Branch=character(),Program=character(),Delta_P2_vs_P1=numeric(),
                           Delta_P3_vs_P1=numeric(),F_statistic=numeric(),P_value=numeric(),FDR_global=numeric())
  csv(tests,'stage_condition_interaction_tests.csv'); csv(bases,'CK_development_module_F_tests.csv')
  csv(bind(statuses),'stage_interaction_status.csv'); csv(bind(missing),'stage_interaction_missing_groups.csv')
  # Within each branch/stage, compare two treatment replicates against two controls; test_scores records missing groups and all comparisons share BH adjustment.
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
  # Display continuous trends separately by branch rather than averaging distinct branch pseudotimes into an artificial common time.
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
  if(!dir.exists(OUT_DIR) && !dir.create(OUT_DIR,recursive=TRUE)) stop("Cannot create the output directory.")
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
  note("Building a network across target cell types with three fixed Endodermis expression clusters P1/P2/P3; output directory=",OUT_DIR)
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
  # Prefer org.Gmax.eg.db with automatic ID diagnostics; failure does not change thresholds or invalidate completed network outputs.
  go <- if(RUN_GO) tryCatch(go_analysis(net$modules,net$background,csv),error=function(e) {
    csv(data.frame(Status='Error',Message=conditionMessage(e)),'GO_error.csv'); note("GO incomplete: ",conditionMessage(e)); NULL
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
  note("Completed. Inspect the three-cluster plots, previous_Program02_P1_evidence, new_modules_P1_evidence, compact GO plots and pipeline_status.")
  csv(data.frame(Status='Core_pipeline_complete',Time=as.character(Sys.time()),Version=PIPELINE_VERSION,
                 GO=if(!RUN_GO)'Disabled' else if(is.null(go))'Error_see_GO_error' else 'Completed',
                 Note='Inspect stage_interaction_status and ME_status for optional/coverage-dependent results'),
      'pipeline_status.csv')
  finished <- TRUE
  invisible(list(trajectory=trajectory,modules=net$modules,scores=scores,GO=go))
}
if(!isTRUE(getOption('soybean.dev.skip_main',FALSE))) main()

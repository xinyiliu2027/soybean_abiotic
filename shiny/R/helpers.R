# Data validation and analysis helpers. Loading this file requires only base R.

normalise_label <- function(x) tolower(trimws(as.character(x)))

resolve_project_groups <- function(md, group_col = NULL, sample_map = setNames(character(), character())) {
  key <- function(x) gsub("[^a-z0-9]", "", normalise_label(x))
  aliases <- list(CK = c("CK", "CK1", "CK2", "CK_1", "CK_2", "control", "control1", "control2"), Heat = c("Heat", "Heat1", "Heat2", "highTemp", "highTemp1", "highTemp2", "37C", "8.37temp"), Cold = c("Cold", "Cold1", "Cold2", "lowTemp", "lowTemp1", "lowTemp2", "6C", "9.6temp"), PEG = c("PEG", "PEG1", "PEG2"))
  lookup <- setNames(rep(names(aliases), lengths(aliases)), key(unlist(aliases, use.names = FALSE)))
  if (length(sample_map)) {
    if (is.null(names(sample_map)) || any(!sample_map %in% names(aliases))) stop("Invalid exact sample-group map.")
    lookup[key(names(sample_map))] <- unname(sample_map)
  }
  candidates <- if (!is.null(group_col)) group_col else intersect(c("PublicationGroup", "group_id", "Group1", "Group", "group", "condition", "treatment", "Treatment", "orig.ident", "sample", "sample_id", "sample_name"), names(md))
  if (!length(candidates) || !all(candidates %in% names(md))) stop("Treatment metadata column not found.")
  groups <- rep(NA_character_, nrow(md)); used <- character()
  for (nm in candidates) {
    z <- unname(lookup[key(md[[nm]])])
    if (any(!is.na(groups) & !is.na(z) & groups != z)) stop("Conflicting treatment metadata; set group_col explicitly.")
    take <- is.na(groups) & !is.na(z)
    if (any(take)) {groups[take] <- z[take]; used <- c(used, nm)}
  }
  if (!any(groups == "CK", na.rm = TRUE)) stop("No CK labels identified. Configure group_col or sample_group_map; REF is not assumed to be CK.")
  list(groups = groups, keep = !is.na(groups), source = paste(used, collapse = ", "))
}

select_ck_cells <- function(md, group_col = NULL, ck_values, ck_only_object = FALSE) {
  match_ck <- function(x) !is.na(x) & normalise_label(x) %in% normalise_label(ck_values)
  candidates <- intersect(c("PublicationGroup", "group_id", "Group1", "Group", "group", "condition", "treatment", "Treatment", "orig.ident", "sample", "sample_id"), colnames(md))
  if (!is.null(group_col)) {
    if (!group_col %in% colnames(md)) stop("Configured group column does not exist.")
    keep <- match_ck(md[[group_col]])
    if (!any(keep)) stop("No configured CK values were found in the group column.")
    return(list(keep = keep, source = group_col))
  }
  hits <- candidates[vapply(candidates, function(nm) any(match_ck(md[[nm]])), logical(1))]
  if (length(hits)) {
    masks <- lapply(hits, function(nm) match_ck(md[[nm]]))
    if (!all(vapply(masks, identical, logical(1), masks[[1]]))) stop("CK selection differs between metadata columns. Set group_col explicitly.")
    return(list(keep = masks[[1]], source = hits[1]))
  }
  if (!isTRUE(ck_only_object)) stop("CK metadata not identified. Configure group_col / ck_values.")
  labels <- unlist(lapply(candidates, function(nm) normalise_label(md[[nm]])), use.names = FALSE)
  if (any(grepl("^(heat|cold|peg|hightemp|lowtemp|aba|iaa|sa|meja|ebl)([_.-]?[0-9]+)?$", labels), na.rm = TRUE)) stop("Non-CK conditions found; disable ck_only_object and select CK explicitly.")
  list(keep = rep(TRUE, nrow(md)), source = "explicit dedicated-CK-reference assumption")
}

match_features <- function(genes, available, feature_map = setNames(character(), character()), allow_missing = FALSE) {
  key <- function(x) tolower(sub("\\.[0-9]+$", "", x))
  setNames(vapply(genes, function(g) {
    if (g %in% names(feature_map)) {
      mapped <- unname(feature_map[g])
      if (!mapped %in% available) stop("Mapped feature missing: ", mapped)
      return(mapped)
    }
    if (g %in% available) return(g)
    hit <- available[key(available) == key(g)]
    if (!length(hit) && isTRUE(allow_missing)) return(NA_character_)
    if (length(hit) != 1L) stop(if (length(hit)) "Ambiguous feature: " else "Missing feature: ", g)
    hit
  }, character(1)), genes)
}

validate_de_table <- function(x) {
  required <- c("gene", "celltype", "treatment", "logFC", "padj")
  if (!all(required %in% names(x))) stop("Required CSV columns: ", paste(required, collapse = ", "))
  x <- x[, required, drop = FALSE]
  if (!nrow(x)) stop("The difference table has no rows.")
  for (nm in c("gene", "celltype", "treatment")) {
    x[[nm]] <- trimws(as.character(x[[nm]]))
    if (anyNA(x[[nm]]) || any(!nzchar(x[[nm]]))) stop("Missing identity values in ", nm)
  }
  if (any(!x$treatment %in% c("Heat", "Cold", "PEG"))) stop("treatment must be Heat, Cold, or PEG, each relative to CK.")
  for (nm in c("logFC", "padj")) {
    x[[nm]] <- suppressWarnings(as.numeric(as.character(x[[nm]])))
    if (any(!is.finite(x[[nm]]))) stop("Non-numeric or missing values in ", nm)
  }
  if (any(x$padj < 0 | x$padj > 1)) stop("padj must be between 0 and 1.")
  if (anyDuplicated(x[, c("gene", "celltype", "treatment")])) stop("Duplicate gene/celltype/treatment: separate time points or models before importing.")
  x
}

shared_response <- function(x, fdr = .05, lfc = .5, celltype = "all", same_direction = TRUE) {
  if (length(fdr) != 1L || !is.finite(fdr) || fdr < 0 || fdr > 1 || length(lfc) != 1L || !is.finite(lfc) || lfc <= 0) stop("FDR must be 0-1 and the absolute log2FC threshold must be positive.")
  x <- validate_de_table(x)
  if (celltype != "all") x <- x[x$celltype == celltype, , drop = FALSE]
  tabs <- lapply(c("Heat", "Cold", "PEG"), function(treatment) {
    d <- x[x$treatment == treatment & x$padj <= fdr & abs(x$logFC) >= lfc, c("gene", "celltype", "logFC", "padj"), drop = FALSE]
    names(d)[3:4] <- paste0(treatment, c("_logFC", "_padj"))
    d
  })
  out <- Reduce(function(a, b) merge(a, b, by = c("gene", "celltype")), tabs)
  common <- with(out, sign(Heat_logFC) == sign(Cold_logFC) & sign(Heat_logFC) == sign(PEG_logFC))
  if (same_direction) out <- out[common, , drop = FALSE]
  out$max_padj <- with(out, pmax(Heat_padj, Cold_padj, PEG_padj))
  common <- with(out, sign(Heat_logFC) == sign(Cold_logFC) & sign(Heat_logFC) == sign(PEG_logFC))
  out$direction <- ifelse(common, ifelse(out$Heat_logFC > 0, "Common up", "Common down"), "Mixed direction")
  out[order(out$max_padj, out$gene), , drop = FALSE]
}

load_ck_atlas <- function(cfg) {
  if (!file.exists(cfg$rds)) stop("Configured CK RDS is not accessible to this R process: ", cfg$rds)
  obj <- readRDS(cfg$rds)
  if (!inherits(obj, "Seurat")) stop("RDS must contain a Seurat object.")
  if (!cfg$assay %in% names(obj@assays)) stop("RNA assay missing.")
  md <- obj[[]]
  selected <- if (isTRUE(cfg$all_treatments)) resolve_project_groups(md, cfg$group_col, cfg$sample_group_map) else select_ck_cells(md, cfg$group_col, cfg$ck_values, cfg$ck_only_object)
  if (isTRUE(cfg$all_treatments)) selected$keep <- selected$keep & selected$groups %in% cfg$treatments
  selected_groups <- setNames(if (isTRUE(cfg$all_treatments)) selected$groups else rep("CK", nrow(md)), rownames(md))
  if (!is.null(cfg$df_column)) {
    if (!cfg$df_column %in% names(md)) stop("DoubletFinder classification column missing.")
    status <- normalise_label(md[[cfg$df_column]])
    if (any(!is.na(status) & !status %in% c("singlet", "doublet", ""))) stop("Unexpected DoubletFinder classification.")
    selected$keep <- selected$keep & !is.na(status) & status == "singlet"
  }
  if (!any(selected$keep)) stop("No CK cells remain.")
  obj <- subset(obj, cells = rownames(md)[selected$keep])
  SeuratObject::DefaultAssay(obj) <- cfg$assay
  ct <- cfg$celltype_col
  if (is.null(ct)) {
    candidates <- intersect(c("PublicationCellType", "celltype_final", "cell_type_final", "cell_type", "celltype", "CellType", "predicted.cell_type", "predicted.celltype", "annotation", "predicted.id"), names(obj[[]]))
    candidates <- candidates[vapply(candidates, function(nm) any(!is.na(obj[[]][[nm]]) & nzchar(trimws(as.character(obj[[]][[nm]])))), logical(1))]
    if (!length(candidates)) stop("Cell-type annotations missing; set celltype_col.")
    ct <- candidates[1]
  }
  if (!ct %in% names(obj[[]])) stop("Configured celltype_col is missing.")
  labels <- trimws(as.character(obj[[]][[ct]]))
  labels[is.na(labels) | !nzchar(labels)] <- "Unannotated"
  reduction <- cfg$reduction
  if (is.null(reduction)) {
    reductions <- intersect(c("umap", "ref.umap", "refumap"), names(obj@reductions))
    if (!length(reductions)) stop("Existing UMAP not found; set reduction explicitly.")
    reduction <- reductions[1]
  }
  if (!reduction %in% names(obj@reductions)) stop("Configured reduction not found.")
  emb <- SeuratObject::Embeddings(obj, reduction = reduction)
  if (ncol(emb) < 2L || !all(colnames(obj) %in% rownames(emb))) stop("UMAP and cell identities do not match.")
  emb <- emb[colnames(obj), 1:2, drop = FALSE]
  if (any(!is.finite(emb))) stop("Invalid UMAP coordinates.")
  if (inherits(obj[[cfg$assay]], "Assay5")) {
    obj <- SeuratObject::JoinLayers(obj, assay = cfg$assay)
    normalised_here <- !"data" %in% SeuratObject::Layers(obj[[cfg$assay]])
    if (normalised_here) obj <- Seurat::NormalizeData(obj, assay = cfg$assay, verbose = FALSE)
    mat <- SeuratObject::GetAssayData(obj, assay = cfg$assay, layer = "data")
  } else {
    mat <- SeuratObject::GetAssayData(obj, assay = cfg$assay, slot = "data")
    normalised_here <- !nrow(mat) || !ncol(mat)
    if (normalised_here) {
      obj <- Seurat::NormalizeData(obj, assay = cfg$assay, verbose = FALSE)
      mat <- SeuratObject::GetAssayData(obj, assay = cfg$assay, slot = "data")
    }
  }
  mat <- mat[, colnames(obj), drop = FALSE]
  features <- match_features(cfg$default_genes, rownames(mat), cfg$feature_map, allow_missing = TRUE)
  message("Loaded ", ncol(mat), " selected cells; group source: ", selected$source, "; annotations: ", ct, "; reduction: ", reduction)
  list(cells = data.frame(id = colnames(obj), x = emb[, 1], y = emb[, 2], celltype = labels, group = unname(selected_groups[colnames(obj)]), stringsAsFactors = FALSE), matrix = mat, features = features,
       group_source = selected$source, celltype_column = ct, reduction = reduction,
       expression_label = if (normalised_here) "RNA LogNormalize: log1p(CP10K)" else "Existing RNA normalized data; original normalization retained")
}

gene_frame <- function(atlas, feature, celltype = "all", group = "CK") {
  if (!feature %in% rownames(atlas$matrix)) stop("Gene not found in the RNA assay.")
  idx <- if (celltype == "all") seq_len(nrow(atlas$cells)) else which(atlas$cells$celltype == celltype)
  if (group != "all" && "group" %in% names(atlas$cells)) idx <- idx[atlas$cells$group[idx] == group]
  if (!length(idx)) stop("No cells match the selected cell type.")
  d <- atlas$cells[idx, , drop = FALSE]
  d$expression <- as.numeric(atlas$matrix[feature, idx, drop = TRUE])
  if (any(!is.finite(d$expression)) || any(d$expression < 0)) stop("Expression must be finite nonnegative normalized RNA data.")
  d[order(d$expression), , drop = FALSE]
}

expression_cap <- function(values, q = .99) {
  positive <- values[is.finite(values) & values > 0]
  if (length(positive)) as.numeric(stats::quantile(positive, probs = q, names = FALSE)) else 1
}

celltype_palette <- function(types) {
  colors <- c("Cortex" = "#2C80C8", "Dividing Cell" = "#F58A17", "Endodermis" = "#BD70A2", "Epidermis" = "#91746A", "Procambium" = "#098577", "Meristem" = "#79AB39", "Pericycle" = "#2B7E37", "Phloem" = "#DE3035", "QC" = "#BF2327", "Root Cap" = "#8491B2", "Xylem" = "#7F3298")
  extra <- setdiff(types, names(colors))
  if (length(extra)) colors <- c(colors, setNames(grDevices::hcl.colors(length(extra), "Dark 3"), extra))
  colors[types]
}

publication_plot <- function(atlas, feature, label = feature, celltype = "all", point_size = .3, q = .99, dpi = 600L, mode = "expression", group = "CK") {
  d <- gene_frame(atlas, feature, celltype, group)
  p <- ggplot2::ggplot(d, ggplot2::aes(x = x, y = y))
  if (mode == "celltype") {
    p <- p + ggrastr::rasterise(ggplot2::geom_point(ggplot2::aes(color = celltype), size = point_size), dpi = dpi) + ggplot2::scale_color_manual(values = celltype_palette(sort(unique(atlas$cells$celltype))))
  } else {
    cap <- expression_cap(d$expression, q)
    p <- p + ggrastr::rasterise(ggplot2::geom_point(color = "#D9DDDA", size = point_size), dpi = dpi) +
      ggrastr::rasterise(ggplot2::geom_point(data = d[d$expression > 0, , drop = FALSE], ggplot2::aes(color = expression), size = point_size), dpi = dpi) +
      ggplot2::scale_color_gradientn(colours = c("#ECE9A6", "#E7AF66", "#B75A6B", "#45316D"), limits = c(0, cap), oob = function(x, range) pmin(pmax(x, range[1]), range[2]), name = "RNA expression")
  }
  p + ggplot2::coord_fixed(xlim = range(atlas$cells$x), ylim = range(atlas$cells$y)) +
    ggplot2::labs(title = label, subtitle = group, x = "UMAP_1", y = "UMAP_2") +
    ggplot2::theme_classic(base_size = 13, base_family = "Helvetica") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = .5, face = "bold"), plot.subtitle = ggplot2::element_text(hjust = .5), legend.title = ggplot2::element_text(size = 9), legend.text = ggplot2::element_text(size = 8))
}

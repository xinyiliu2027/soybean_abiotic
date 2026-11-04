# Source: 图4.docx, extracted text lines 1989-2173
# 读取01生成的模块评分RDS及GO结果
# See README.md and docs/KNOWN_LIMITATIONS.md before running.

###
suppressPackageStartupMessages({
  library(Seurat)
  library(hdWGCNA)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(Matrix)
})

# ==============================================================================
# 1. 基础设置与加载数据 (不改输出目录)
# ==============================================================================
OUT_DIR <- file.path(Sys.getenv("SOYBEAN_PROJECT_ROOT", unset = "data"), "CK2_added_v3_Results/Ensemble_PI_Heat_Cold_PEG_with_CK2")
FILE_PREFIX <- "Endo3_v5_"
NETWORK_NAME <- "AllCells_Endo3_consensus_v5"

setwd(OUT_DIR)
seurat_obj <- readRDS(file.path(OUT_DIR, paste0(FILE_PREFIX, "all_cells_modules_scored.rds")))

# ==============================================================================
# 2. 自定义高级科研调色板 (按模块顺序 M1-M11 分配)
# ==============================================================================
# 根据你提供的 9 个颜色，并补充 M5、M7 的高级过渡色
palette_research <- c(
  '#A4373A', # M1 砖红
  '#1B6A6D', # M2 深青
  '#4B2D5D', # M3 暗紫
  '#99345C', # M4 紫红
  '#8C6D31', # M5 (原图缺，补充暗金)
  '#2B6A3E', # M6 深绿
  '#7A5B7E', # M7 (原图缺，补充灰紫)
  '#B28A3A', # M8 暗黄/金
  '#326477', # M9 青灰
  '#C35A3C', # M10 橙红
  '#788A41'  # M11 橄榄绿
)

module_color <- function(program) {
  n <- suppressWarnings(as.integer(gsub('[^0-9]', '', program)))
  if (is.na(n)) return('#5B6674') # 兜底灰
  palette_research[(n - 1L) %% length(palette_research) + 1L]
}

module_limits <- function(value) {
  value <- value[is.finite(value)]
  if(!length(value)) return(c(0,1))
  lim <- as.numeric(stats::quantile(value, c(0.01, 0.99), names = FALSE))
  if(diff(lim) < 1e-10) lim <- range(value)
  if(diff(lim) < 1e-10) lim <- lim + c(-.5, .5)
  lim
}

# ==============================================================================
# 3. 强行重写模块颜色 (解决网络图全是红色的问题)
# ==============================================================================
mods <- hdWGCNA::GetModules(seurat_obj, wgcna_name = NETWORK_NAME)
# 用新调色板覆盖模块表中的颜色列
mods$color <- sapply(as.character(mods$module), module_color)
seurat_obj <- hdWGCNA::SetModules(seurat_obj, mods, wgcna_name = NETWORK_NAME)

# ==============================================================================
# 4. 画 Top 50 模块网络图 (修复 TOM 路径 bug)
# ==============================================================================
relative_tom_path <- paste0(FILE_PREFIX, "TOM/ERS_consensus_TOM.rda")
seurat_obj@misc[[NETWORK_NAME]]$TOMFiles <- relative_tom_path
if (!is.null(seurat_obj@misc[[NETWORK_NAME]]$wgcna_net)) {
  seurat_obj@misc[[NETWORK_NAME]]$wgcna_net$TOMFiles <- relative_tom_path
}

plot_dir <- file.path(OUT_DIR, "ModuleNetworks_ResearchStyle")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

cat("开始绘制 StageProgram09 网络图 (Top 50)...\n")
p_net <- hdWGCNA::ModuleNetworkPlot(
  seurat_obj = seurat_obj,
  modules = "StageProgram09",  
  outdir = plot_dir,
  n_inner = 20,             
  n_outer = 30,             
  n_conns = Inf,            
  plot_size = c(10, 10),    
  vertex.label.cex = 0.8    
)

# ==============================================================================
# 5. 画 UMAP 图 (浅灰到模块专属色)
# ==============================================================================
sc <- hdWGCNA::GetModuleScores(seurat_obj, wgcna_name = NETWORK_NAME)
md <- seurat_obj@meta.data
programs <- setdiff(colnames(sc), c('grey', 'MEgrey'))
programs <- sort(programs) # 保证 M1 到 M11 顺序

plot_atlas_style <- function(mat, md, prefix, label) {
  mat <- as.matrix(mat)
  set.seed(20261026)
  # 抽样绘制，防止卡死
  ix <- unlist(lapply(split(seq_len(nrow(md)), md$Condition), function(i)
    i[sample.int(length(i), min(length(i), 12000))]), use.names = FALSE)
  q <- md[ix, , drop = FALSE]
  q$Condition <- factor(q$Condition, c('CK', 'Heat', 'Cold', 'PEG'))
  
  limits <- setNames(lapply(programs, function(p) module_limits(mat[, p])), programs)
  
  panel <- function(p) {
    d <- q
    d$Value <- pmin(limits[[p]][2], pmax(limits[[p]][1], mat[d$cell, p]))
    d <- d[order(d$Value), , drop = FALSE]
    
    # 背景极浅灰，高值使用模块专属色
    ggplot2::ggplot(d, ggplot2::aes(UMAP1, UMAP2, color = Value)) +
      ggplot2::geom_point(size = .2) + 
      ggplot2::coord_equal() +
      ggplot2::scale_color_gradient(low = '#F0F0F0', high = module_color(p), limits = limits[[p]]) +
      ggplot2::theme_void(base_size = 12) +
      ggplot2::labs(title = p, color = label) +
      ggplot2::facet_wrap(~Condition, nrow = 1) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 12, face = 'bold'),
        legend.title = ggplot2::element_text(size = 10),
        strip.text = ggplot2::element_text(size = 12)
      )
  }
  
  # 分页保存 UMAP (每页 4 个模块)
  pages <- split(programs, ceiling(seq_along(programs)/4))
  for(i in seq_along(pages)) {
    plots <- lapply(pages[[i]], panel)
    p_umap <- patchwork::wrap_plots(plots, ncol = 1)
    ggsave(file.path(OUT_DIR, paste0(prefix, '_page', i, '.pdf')), p_umap, width = 14, height = 2.8 * length(plots))
    ggsave(file.path(OUT_DIR, paste0(prefix, '_page', i, '.png')), p_umap, width = 14, height = 2.8 * length(plots), dpi = 300)
  }
}

cat("开始绘制 UMAP...\n")
plot_atlas_style(sc, md, 'Figure_UMAP_ResearchStyle', 'Module Score')

# ==============================================================================
# 6. 画 GO 气泡图 (模块专属色渐变高亮)
# ==============================================================================
go_file <- file.path(OUT_DIR, paste0(FILE_PREFIX, "GO_BP_significant.csv"))
if (file.exists(go_file)) {
  cat("开始绘制 GO 气泡图...\n")
  go_sig <- read.csv(go_file, stringsAsFactors = FALSE)
  
  if (nrow(go_sig) > 0) {
    # 取每个模块前 5 的条目
    q_go <- go_sig %>%
      group_by(Program) %>%
      arrange(FDR_global, desc(Count)) %>%
      slice_head(n = 5) %>%
      ungroup()
    
    panels <- lapply(split(q_go, q_go$Program), function(a) {
      a <- a[order(a$GeneRatio), , drop = FALSE]
      a$Term <- factor(a$Description, levels = unique(a$Description))
      p_name <- as.character(a$Program[1])
      
      # 气泡渐变色由浅灰到该模块的专属色
      ggplot2::ggplot(a, ggplot2::aes(GeneRatio, Term, size = Count, color = -log10(pmax(FDR_global, 1e-300)))) +
        ggplot2::geom_point() + 
        ggplot2::scale_color_gradient(low = '#E0E0E0', high = module_color(p_name)) +
        ggplot2::scale_size_continuous(range = c(2, 6)) +
        ggplot2::theme_bw(base_size = 12) +
        ggplot2::theme(
          panel.grid.minor = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_text(size = 10, color = 'black'),
          plot.title = ggplot2::element_text(size = 12, face = 'bold', color = module_color(p_name))
        ) +
        ggplot2::labs(title = paste(p_name, 'GO BP'), x = 'Gene ratio', y = NULL, color = '-log10(FDR)', size = 'Genes')
    })
    
    # 每页两个图合并
    pages_go <- split(names(panels), ceiling(seq_along(panels)/2))
    for(i in seq_along(pages_go)) {
      p_go <- patchwork::wrap_plots(panels[pages_go[[i]]], ncol = 1)
      ggsave(file.path(OUT_DIR, paste0('Figure_GO_ResearchStyle_page', i, '.pdf')), p_go, width = 8, height = 8)
      ggsave(file.path(OUT_DIR, paste0('Figure_GO_ResearchStyle_page', i, '.png')), p_go, width = 8, height = 8, dpi = 300)
    }
    cat("GO 气泡图已保存！\n")
  }
}

cat("\n✅ 所有绘图完成！请去目录下查看。\n")

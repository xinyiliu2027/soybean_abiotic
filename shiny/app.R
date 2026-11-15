required <- c("shiny", "Seurat", "SeuratObject", "Matrix", "ggplot2", "plotly", "DT", "patchwork", "ggrastr", "jsonlite", "htmltools")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Run install_packages.R first. Missing packages: ", paste(missing, collapse = ", "))
library(shiny)
source("config.R", local = TRUE)
source(file.path("R", "helpers.R"), local = TRUE)
options(shiny.maxRequestSize = 50 * 1024^2)
set.seed(config$seed)
# Load once per R process, not once per visitor. Keep the RNA matrix sparse.
shared_atlas <- tryCatch(load_ck_atlas(config), error = function(e) {
  message("CK atlas load failed: ", conditionMessage(e)); list(load_error = conditionMessage(e))
})

shared_de <- if (file.exists("data/stress_de.csv.gz")) validate_de_table(utils::read.csv(gzfile("data/stress_de.csv.gz"), stringsAsFactors=FALSE, check.names=FALSE)) else if (file.exists("data/stress_de.csv")) validate_de_table(utils::read.csv("data/stress_de.csv", stringsAsFactors=FALSE, check.names=FALSE)) else NULL
ui <- fluidPage(
  tags$head(tags$title("single_cell_abiotic_soybean"), tags$style(HTML("
    body{background:#f5f8f3;color:#223b2e;font-family:'Segoe UI','Microsoft YaHei',sans-serif}
    .container-fluid{max-width:1500px;padding:20px 35px}.hero{padding:20px 0 25px}.hero h1{font-size:31px;font-weight:600;margin:8px 0 15px}.eyebrow{letter-spacing:2px;color:#3d6e4c;font-size:11px}.hero p{color:#718072;max-width:850px}.nav-tabs{border-bottom:1px solid #dce5d9}.nav-tabs>li>a{color:#326347}.nav-tabs>li.active>a{background:#fff;border-top:3px solid #376949}.well{background:#fff;border:1px solid #dce5d9;box-shadow:none;border-radius:6px;margin-top:20px}.panel-area{background:white;border:1px solid #dce5d9;border-radius:6px;padding:20px;margin-top:20px}.btn{border-radius:4px}.btn-primary{background:#326347;border-color:#326347}.preset{display:block;width:100%;text-align:left;font-size:11px;margin-bottom:7px}.status-box{background:#eaf1e4;padding:13px 18px;border-radius:5px;margin:0 0 23px}.muted{color:#718072;font-size:12px}.downloads{display:flex;flex-wrap:wrap;gap:10px;margin:12px 0 20px}.metric{display:inline-block;margin-right:25px;font-size:12px}.section-title{font-size:20px;font-weight:600}.methods{max-width:950px;padding:25px 0}.methods li{margin-bottom:10px}.shiny-output-error-validation{color:#718072}.footer{border-top:1px solid #dce5d9;font-size:11px;color:#718072;margin-top:35px;padding:20px 0}
  "))),
  div(class = "institution-bar", tags$a(href="https://ics.caas.cn/en/",target="_blank",tags$img(src="caas.png",height="50",alt="CAAS emblem")," Institute of Crop Sciences, CAAS"),tags$a(href="https://bri.caas.cn/en/",target="_blank",tags$img(src="caas.png",height="50",alt="CAAS emblem")," Biotechnology Research Institute, CAAS")), div(class = "hero", div(class = "eyebrow", "GLYCINE MAX · SINGLE-CELL STRESS ATLAS"), h1("single_cell_abiotic_soybean"), p("Soybean stress single-cell atlas: locate CK marker expression and explore shared responses to heat, cold, and PEG.")),
  uiOutput("data_status"),
  tabsetPanel(id = "main_tab",
    tabPanel("Gene UMAP", sidebarLayout(
      sidebarPanel(width = 3,
        h4("Genes & cells"), selectInput("treatment", "Treatment", choices = c("CK", "Heat", "Cold", "PEG"), selected = "CK"), selectizeInput("gene", "Search RNA genes", choices = NULL),
        tags$label("Selected in situ gene"), lapply(seq_along(config$default_genes), function(i) actionButton(paste0("preset_", i), config$default_genes[i], class = "preset")),
        selectInput("celltype", "Cell type", choices = c("All cell types" = "all")),
        selectInput("color_mode", "Color by", choices = c("Gene expression" = "expression", "Cell type" = "celltype")),
        sliderInput("point_size", "Interactive point size", min = 1, max = 6, value = 2, step = .5),
        selectInput("quantile", "Positive-expression color limit", choices = c("99th percentile" = .99, "95th percentile" = .95, "Maximum" = 1)),
        p(class = "muted", "Gray indicates zero expression. Each gene has its own color scale; colors cannot be compared directly across genes."),
        h4("Publication plot settings"), sliderInput("pub_point", "Point size", min = .1, max = .6, value = .3, step = .05),
        p(class = "muted", "Downloaded PDFs rasterize scatter layers at 600 dpi while retaining vector text and axes.")
      ),
      mainPanel(width = 9, div(class = "panel-area", uiOutput("plot_heading"), plotly::plotlyOutput("umap", height = "610px"),
        div(class = "downloads", downloadButton("plot_pdf", "Current plot PDF · 600 dpi"), downloadButton("plot_png", "Current plot PNG · 600 dpi"), downloadButton("cell_csv", "Selected-cell expression CSV")), DT::DTOutput("expression_summary")))
    )),
    tabPanel("Five-gene overview", div(class = "panel-area", h3("CK cell types and five selected in situ genes"), p(class = "muted", "Uses all CK cells with independent gene color scales. The single-gene cell-type filter does not affect this overview."), plotOutput("overview", height = "780px"),
      div(class = "downloads", downloadButton("overview_pdf", "Six-panel PDF · 600 dpi"), downloadButton("overview_png", "Six-panel PNG · 600 dpi"), downloadButton("atlas_json", "Export atlas.json for the static website")))),
    tabPanel("Shared response", div(class = "panel-area", h3("Heat / Cold / PEG shared-response candidates"), p(class="muted", "Bundled source: All_DEG_filtered_paper_thresholds.csv. Original control label: REF; FDR: p_adj.loc. This is a prefiltered table; absent rows are not evidence of no response."), p("Import differential expression relative to CK to find the intersection of significant responses within the same cell type."),
      fluidRow(column(5, fileInput("de_file", "DE results CSV (up to 50 MB)", accept = ".csv")), column(3, numericInput("fdr", "FDR ≤", .05, min = 0, max = 1, step = .01)), column(3, numericInput("lfc", "|log2FC| ≥", .5, min = .01, step = .1))),
      fluidRow(column(5, selectInput("de_celltype", "Cell type", choices = c("All (computed separately)" = "all"))), column(5, checkboxInput("same_direction", "Require the same response direction in all three stresses", TRUE))),
      p(class = "muted", "Required fields: gene, celltype, treatment, logFC, padj. Treatment must be Heat, Cold, or PEG. Use replicate-aware pseudobulk / muscat results. Import each time point or model separately."),
      uiOutput("de_message"), DT::DTOutput("de_table"), div(class = "downloads", downloadButton("de_csv", "Export all candidates as CSV")), p(class = "muted", "Shared responses do not establish stress tolerance or breeding effects. Genetic and phenotypic validation is required."))),
    tabPanel("Data & methods", div(class = "methods", h3("Data and analysis scope"), tags$ul(
      tags$li("The configured merged Seurat object is filtered to CK, Heat, Cold, and PEG. Existing shared UMAP coordinates are retained; independently fitted sample UMAPs are not concatenated."),
      tags$li("Existing UMAP coordinates and cell-type annotations are retained. Expression uses normalized RNA data, not the integrated assay or scale.data."),
      tags$li("LogNormalize is run only if RNA data are missing. The data provider must confirm any existing normalization."),
      tags$li("DoubletFinder is not rerun. Confirmed singlets are retained only when an existing classification column is explicitly configured."),
      tags$li("Shared responses require all three treatments to pass FDR and log2FC thresholds within the same cell type. Missing comparisons are not filled with zeros.")),
      uiOutput("provenance"),
      p("The exploration concept was informed by the Zhai lab Arabidopsis hormone single-nucleus study. This site does not claim affiliation with that team or institution."),
      tags$a(href = "https://www.nature.com/articles/s41467-026-76284-y", target = "_blank", rel = "noopener", "Nature Communications · Reference")))
  ), div(class = "footer", "single_cell_abiotic_soybean · CK / Heat / Cold / PEG")
)

server <- function(input, output, session) {
  atlas <- reactive({validate(need(is.null(shared_atlas$load_error), "No dataset loaded. Configure the merged Seurat RDS on the server and restart the app.")); shared_atlas})
  output$data_status <- renderUI({
    if (!is.null(shared_atlas$load_error)) return(div(class = "status-box", "Waiting for CK data · No expression results displayed"))
    a <- atlas(); div(class = "status-box", span(class = "metric", paste(format(nrow(a$cells), big.mark = ","), "cells across CK / Heat / Cold / PEG")), span(class = "metric", paste(length(unique(a$cells$celltype)), "cell types")), span(class = "metric", "5 selected in situ genes"))
  })
  if (is.null(shared_atlas$load_error)) {
    all_features <- rownames(shared_atlas$matrix)
    updateSelectizeInput(session, "gene", choices = all_features, selected = unname(shared_atlas$features[1]), server = TRUE)
    updateSelectInput(session, "treatment", choices = intersect(c("CK", "Heat", "Cold", "PEG"), unique(shared_atlas$cells$group)), selected="CK")
    updateSelectInput(session, "celltype", choices = c("All cell types" = "all", setNames(sort(unique(shared_atlas$cells$celltype)), sort(unique(shared_atlas$cells$celltype)))))
  }
  for (i in seq_along(config$default_genes)) local({
    index <- i
    observeEvent(input[[paste0("preset_", index)]], {
      a <- atlas(); updateSelectizeInput(session, "gene", selected = unname(a$features[index]))
    })
  })
  selected_data <- reactive({req(input$gene, input$celltype); gene_frame(atlas(), input$gene, input$celltype, input$treatment)})
  output$plot_heading <- renderUI({d <- selected_data(); tagList(h3(input$gene), p(class = "muted", paste(input$treatment, "|", nrow(d), "cells |", atlas()$expression_label)))})
  output$umap <- plotly::renderPlotly({
    d <- selected_data(); a <- atlas()
    # Escape annotation/ID text before using Plotly HTML hover labels.
    d$hover <- as.character(htmltools::htmlEscape(paste0("Cell: ", d$id, " | Type: ", d$celltype, " | Gene: ", input$gene, " | Expression: ", signif(d$expression, 5))))
    if (input$color_mode == "celltype") {
      colors <- celltype_palette(sort(unique(d$celltype)))
      d$celltype <- factor(d$celltype, levels = names(colors))
      p <- plotly::plot_ly(d, x = ~x, y = ~y, type = "scattergl", mode = "markers", color = ~celltype, colors = unname(colors), text = ~hover, hoverinfo = "text", marker = list(size = input$point_size, opacity = 1))
    } else {
      cap <- expression_cap(d$expression, as.numeric(input$quantile))
      zero <- d[d$expression <= 0, , drop = FALSE]; pos <- d[d$expression > 0, , drop = FALSE]
      p <- plotly::plot_ly()
      if (nrow(zero)) p <- plotly::add_trace(p, data = zero, x = ~x, y = ~y, type = "scattergl", mode = "markers", text = ~hover, hoverinfo = "text", marker = list(color = "#D9DDDA", size = input$point_size, opacity = 1), showlegend = FALSE)
      if (nrow(pos)) p <- plotly::add_trace(p, data = pos, x = ~x, y = ~y, type = "scattergl", mode = "markers", text = ~hover, hoverinfo = "text", marker = list(color = pmin(pos$expression, cap), colorscale = list(list(0, "#ECE9A6"), list(.33, "#E7AF66"), list(.67, "#B75A6B"), list(1, "#45316D")), cmin = 0, cmax = cap, showscale = TRUE, colorbar = list(title = "RNA expression"), size = input$point_size, opacity = 1), showlegend = FALSE)
    }
    p <- plotly::layout(p, xaxis = list(title = "UMAP_1", range = range(a$cells$x), zeroline = FALSE), yaxis = list(title = "UMAP_2", range = range(a$cells$y), scaleanchor = "x", scaleratio = 1, zeroline = FALSE), paper_bgcolor = "white", plot_bgcolor = "white", margin = list(l = 55, r = 90, b = 50, t = 20), dragmode = "pan")
    plotly::config(p, scrollZoom = TRUE, displaylogo = FALSE)
  })
  output$expression_summary <- DT::renderDT({
    d <- selected_data(); rows <- lapply(split(d$expression, d$celltype), function(v) data.frame(cells = length(v), mean_expression = mean(v), pct_expressing = 100 * mean(v > 0)))
    summary <- do.call(rbind, rows); summary$celltype <- rownames(summary); rownames(summary) <- NULL
    summary <- summary[, c("celltype", "cells", "mean_expression", "pct_expressing")]
    DT::datatable(summary, rownames = FALSE, escape = TRUE, options = list(pageLength = 12, dom = "tip")) |> DT::formatRound(c("mean_expression", "pct_expressing"), 3)
  })
  current_plot <- reactive({req(input$gene); publication_plot(atlas(), input$gene, celltype = input$celltype, point_size = input$pub_point, q = as.numeric(input$quantile), dpi = config$raster_dpi, mode = input$color_mode, group = input$treatment)})
  filename <- function(suffix) paste0(input$treatment, "_", gsub("[^A-Za-z0-9_.-]", "_", input$gene), suffix)
  output$plot_pdf <- downloadHandler(filename = function() filename(".pdf"), content = function(file) ggplot2::ggsave(file, current_plot(), device = grDevices::pdf, width = 6, height = 5, useDingbats = FALSE, family = "Helvetica", bg = "white"))
  output$plot_png <- downloadHandler(filename = function() filename(".png"), content = function(file) ggplot2::ggsave(file, current_plot(), device = "png", width = 6, height = 5, dpi = config$raster_dpi, bg = "white"))
  output$cell_csv <- downloadHandler(filename = function() filename("_cells.csv"), content = function(file) utils::write.csv(selected_data(), file, row.names = FALSE, fileEncoding = "UTF-8"))
  overview <- reactive({
    a <- atlas(); ref <- publication_plot(a, unname(a$features[1]), label = "CK cell types", point_size = input$pub_point, dpi = config$raster_dpi, mode = "celltype") + ggplot2::theme(legend.position = "none")
    ps <- lapply(names(a$features), function(g) publication_plot(a, unname(a$features[g]), label = g, point_size = input$pub_point, q = as.numeric(input$quantile), dpi = config$raster_dpi))
    patchwork::wrap_plots(c(list(ref), ps), ncol = 3)
  })
  output$overview <- renderPlot({overview()}, res = 120)
  output$overview_pdf <- downloadHandler(filename = "CK_five_genes_2x3.pdf", content = function(file) ggplot2::ggsave(file, overview(), device = grDevices::pdf, width = 15, height = 9, useDingbats = FALSE, family = "Helvetica", bg = "white"))
  output$overview_png <- downloadHandler(filename = "CK_five_genes_2x3.png", content = function(file) ggplot2::ggsave(file, overview(), device = "png", width = 15, height = 9, dpi = config$raster_dpi, bg = "white"))
  output$atlas_json <- downloadHandler(filename = "atlas.json", contentType = "application/json", content = function(file) {
    a <- atlas(); cells <- a$cells
    expression <- setNames(lapply(names(a$features), function(g) I(as.numeric(a$matrix[a$features[g], ]))), names(a$features))
    jsonlite::write_json(list(schema_version = 1L, source_label = basename(config$rds), created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), expression_label = a$expression_label, cells = cells, expression = expression), file, dataframe = "rows", auto_unbox = TRUE, digits = 7)
  })
  de <- reactiveVal(shared_de); de_error <- reactiveVal(NULL)
  if (!is.null(shared_de)) updateSelectInput(session, "de_celltype", choices=c("All (computed separately)"="all", setNames(sort(unique(shared_de$celltype)),sort(unique(shared_de$celltype)))))
  observeEvent(input$de_file, {
    tryCatch({
      d <- validate_de_table(utils::read.csv(input$de_file$datapath, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8-BOM"))
      de(d); de_error(NULL); types <- sort(unique(d$celltype)); updateSelectInput(session, "de_celltype", choices = c("All (computed separately)" = "all", setNames(types, types)))
    }, error = function(e) {de(NULL); de_error(conditionMessage(e))})
  })
  candidates <- reactive({
    validate(need(is.null(de_error()), de_error()), need(!is.null(de()), "Import a differential-expression CSV first."))
    validate(need(length(input$fdr) == 1 && is.finite(input$fdr) && input$fdr >= 0 && input$fdr <= 1, "FDR must be between 0 and 1."), need(length(input$lfc) == 1 && is.finite(input$lfc) && input$lfc > 0, "The log2FC threshold must be positive."))
    shared_response(de(), input$fdr, input$lfc, input$de_celltype, input$same_direction)
  })
  output$de_message <- renderUI({d <- candidates(); p(class = "muted", paste(nrow(d), "gene / cell-type combinations meet the current thresholds."))})
  output$de_table <- DT::renderDT({DT::datatable(candidates(), rownames = FALSE, escape = TRUE, options = list(pageLength = 20, scrollX = TRUE))})
  output$de_csv <- downloadHandler(filename = "shared_stress_candidates.csv", content = function(file) utils::write.csv(candidates(), file, row.names = FALSE, fileEncoding = "UTF-8"))
  output$provenance <- renderUI({
    a <- atlas()
    tags$ul(
      tags$li(paste("Group source: ", a$group_source)),
      tags$li(paste("Annotation column: ", a$celltype_column)),
      tags$li(paste("Reduction: ", a$reduction)),
      tags$li(paste("Expression: ", a$expression_label))
    )
  })
}
shinyApp(ui, server)

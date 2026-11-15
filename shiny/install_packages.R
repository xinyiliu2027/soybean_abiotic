required <- c("shiny", "Seurat", "SeuratObject", "Matrix", "ggplot2", "plotly", "DT", "patchwork", "ggrastr", "jsonlite", "htmltools")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
remaining <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(remaining)) stop("Packages still missing: ", paste(remaining, collapse = ", "))
message("All required packages are available.")

# Connect the merged Seurat object to the English website

Live website: https://xinyiliu2027.github.io/soybean_abiotic/

The initial live dataset now contains real CK (27,042), Heat (34,404), Cold (23,923) and PEG (21,393) cells: 106,762 cells across 11 cell types. The supplied atlas contains five single-cell expression vectors: Glyma.05G103600, Glyma.15G040200, Glyma.12G234200, Glyma.08G163500 and Glyma.10G022500. It does not contain the requested Glyma.20G176100 or Glyma.11G147700. Missing vectors are displayed as unavailable rather than replaced with another gene.

The website also includes 97,265 actual filtered DE rows for Heat/Cold/PEG across 17,797 distinct genes. These are gene-level contrasts and detection frequencies, not an all-gene single-cell expression matrix. The original REF label and local FDR are retained.

## Recommended first step: export from server RStudio

The website is already interactive. You do not need to rebuild it in Shiny to inspect the five selected genes. It needs UMAP coordinates, annotations, treatment labels and selected normalized RNA expression, rather than the entire Seurat RDS.

1. Open an RStudio session on the Linux server that can read `/data01`. A Windows RStudio session cannot read this path unless the data have been transferred or mounted.
2. Save `Export_Seurat_Atlas.R` from the website to that server, or transfer the downloaded script to the server using the RStudio Files pane / your normal file transfer method.
3. Run the following in that RStudio session, replacing only the script location:

```r
source("/path/to/Export_Seurat_Atlas.R", encoding = "UTF-8")
```

The script already uses:

```r
MERGED_RDS <- "/data01/liuxy/Project/CK2_added_v3_Results/seurat_obj_with_CK2_v3.rds"
OUT_DIR <- "/data01/liuxy/Project/soybean_abiotic_exports"
```

4. Download `/data01/liuxy/Project/soybean_abiotic_exports/atlas.json` from the server to your computer.
5. On the website, select **Import atlas JSON** and choose that file. CK is selected by default. The treatment menu contains the groups actually exported; choose Heat, Cold or PEG to browse their expression.

The default export contains the five specified in situ genes and all recognized CK/Heat/Cold/PEG cells. The script uses existing UMAP coordinates and normalized RNA data; it does not fit UMAP or infer DE results from a scatter plot. It retains all selected cells without random subsampling.

The updated exporter also writes `expression_index.json` and a `gene_summaries/` directory for all RNA genes. Each file contains a small batch of genes with mean expression and detection percentage for each treatment/cell-type stratum. Download the entire export directory from RStudio, extract it locally, then choose **Import expression folder** and select the folder containing `expression_index.json`, `atlas.json`, and `gene_summaries/`.

Sparse single-cell vectors for arbitrary-gene UMAP are included when their estimated size is below the configured static-data budget. Larger transcriptomes still export all-gene summaries and the five selected single-cell vectors; use Shiny for arbitrary-gene UMAP in that case. A single huge dense JSON is not required.

If group detection stops, open `group_metadata_audit.csv` and `metadata_columns.csv` in the export directory. Set `GROUP_COL` to the real treatment column, and add exact sample mappings to `SAMPLE_GROUP_MAP` when sample IDs are coded. Conflicting or missing CK labels stop the export; the full merged object is never silently treated as CK. Check `gene_match_report.csv` for missing or ambiguous genes.

## Local import versus a published dataset

**Import atlas JSON** reads the file in the current browser session. It does not upload your data to GitHub. Reloading requires another import unless a dataset has been published.

To make a new export open automatically for all visitors, publish its `atlas.json`, `expression_index.json`, and `gene_summaries/` under `docs/`, preserving relative file paths. The current published initial atlas uses compressed small parts and a manifest; a new `docs/atlas.json` takes precedence. The viewer accepts local files up to 100 MB and 500,000 cells. Keep the same schema and cell order; sparse vectors use a matching cell-order signature.

Only publish an atlas file when you intend to make those data accessible through the public website. There is no need to publish the full `.rds`.

## When to use Shiny

Use the provided Shiny app if you need direct access to the large Seurat RDS and genome-wide RNA gene search. It filters the merged object to CK/Heat/Cold/PEG and defaults to CK; the treatment selector also shows stress expression. The six-panel marker overview remains CK-only. The shared-stress module includes the supplied filtered Heat/Cold/PEG contrasts, retaining their original REF control label.

Extract `soybean_abiotic_shiny.zip` on the data server. In server RStudio:

```r
setwd("/path/to/shiny")
source("install_packages.R") # Only needed when dependencies are missing.
Sys.setenv(SOYBEAN_ATLAS_RDS = "/data01/liuxy/Project/CK2_added_v3_Results/seurat_obj_with_CK2_v3.rds")
shiny::runApp(".", host = "0.0.0.0", port = 3838, launch.browser = FALSE)
```

Open `http://SERVER_IP:3838/` from a computer that can reach that server and port. `0.0.0.0` is a listening address, not the address to enter in the browser. The R process must remain running. For a persistent shared service, use Shiny Server / Posit Connect or your server administrator's supported service setup.

GitHub Pages hosts the static HTML/JavaScript viewer; it cannot run an R Shiny process. RStudio is where you develop and launch Shiny, while the R process/server provides the running website. Cloud hosting needs its own accessible dataset; it cannot automatically read the laboratory server's `/data01` path.

## Differential expression and interpretation

Expression UMAPs do not establish shared stress responses. Import a CSV with `gene,celltype,treatment,logFC,padj`, where treatment is Heat, Cold or PEG, logFC is log2FC, and each contrast is relative to CK. Prefer muscat / pseudobulk analyses that include biological replicates. Separate time points/models and keep FDR correction scopes explicit.

## Validation limits

The English interface, file formats, JavaScript logic, R syntax and pure-R scientific filtering helpers are checked locally. The server RDS remains inaccessible from this computer, so real metadata mappings and a full Seurat run still require confirmation on the data server.

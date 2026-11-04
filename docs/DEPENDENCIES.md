# Dependencies

The source documents do not provide complete environment versions. This inventory is not a tested lockfile.

| Workflow | Main dependencies |
| --- | --- |
| Figure 1 | Seurat with v5 JoinLayers support, SeuratObject, dplyr, ggplot2, patchwork, gridExtra, qs, glmGamPoi; the configured UMAP backend also uses Python umap-learn |
| DoubletFinder example | Seurat, dplyr, DoubletFinder |
| Figure 2 DE | SingleCellExperiment, muscat, edgeR, limma, Seurat, dplyr, tidyr, ggplot2, cowplot |
| Figure 2 PI | Seurat, SeuratObject, Matrix, Augur, scDist, mixtools, future, and dependencies called within the script |
| Figure 3 reference states | Seurat, SeuratObject, Matrix, dplyr, ggplot2, irlba, FNN, cluster, slingshot, ComplexHeatmap, circlize; optional ggrepel; edgeR/limma for inferential DE |
| Figure 3 candidate table | openxlsx |
| Figure 4 network | Seurat, SeuratObject, Matrix, hdWGCNA, WGCNA, ggplot2, patchwork, digest, princurve, irlba, FNN, limma, UCell |
| Figure 4 TF | GENIE3, UCell, BiocParallel, SeuratObject, Matrix, pheatmap, digest; doRNG/doParallel/foreach/dplyr for parallel execution; limma for exploratory tests |
| GO | clusterProfiler, AnnotationDbi, GO.db and an appropriate soybean OrgDb |
| Python velocity | numpy, pandas, scanpy, scvelo, matplotlib, listed in requirements.txt |

Use Bioconductor packages compatible with the selected R release. Install Augur, scDist, hdWGCNA and DoubletFinder according to the instructions for the selected package versions. The scripts refer to both `org.G.max.eg.db` and `org.Gmax.eg.db`; use the package and identifier type that actually match the data's annotation version. Changing a package name alone does not establish compatibility.

No bioinformatics analysis environment was installed or validated during packaging.

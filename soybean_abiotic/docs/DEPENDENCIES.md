# Dependencies

原文没有提供完整环境版本。以下按脚本中的包加载和命名空间调用整理，属于依赖清单而非经过复现的锁文件。

| 分析 | 主要依赖 |
| --- | --- |
| 图1 | Seurat（使用 v5 JoinLayers 接口）、SeuratObject、dplyr、ggplot2、patchwork、gridExtra、qs、glmGamPoi；UMAP 配置涉及 Python umap-learn |
| DoubletFinder 示例 | Seurat、dplyr、DoubletFinder |
| 图2 DE | SingleCellExperiment、muscat、edgeR、limma、Seurat、dplyr、tidyr、ggplot2、cowplot |
| 图2 PI | Seurat、SeuratObject、Matrix、Augur、scDist、mixtools、future，以及脚本中调用的依赖 |
| 图3 CK 分期 | Seurat、SeuratObject、Matrix、dplyr、ggplot2、irlba、FNN、cluster、slingshot、ComplexHeatmap、circlize；可选 ggrepel；DE 使用 edgeR/limma |
| 图3候选表 | openxlsx |
| 图4网络 | Seurat、SeuratObject、Matrix、hdWGCNA、WGCNA、ggplot2、patchwork、digest、princurve、irlba、FNN、limma、UCell |
| 图4 TF | GENIE3、UCell、BiocParallel、SeuratObject、Matrix、pheatmap、digest；并行依赖 doRNG、doParallel、foreach、dplyr；探索统计使用 limma |
| GO | clusterProfiler、AnnotationDbi、GO.db，以及与输入基因版本一致的大豆 OrgDb |
| Python velocity | requirements.txt 中的 numpy、pandas、scanpy、scvelo、matplotlib |

Bioconductor 包使用与当前 R 版本匹配的 Bioconductor 环境；Augur、scDist、hdWGCNA 和 DoubletFinder 应按各包对应版本的安装说明准备。原文出现 `org.G.max.eg.db` 和 `org.Gmax.eg.db` 两种名称，实际需要匹配你的本地注释库和 ID 类型，不能仅替换包名来假定兼容。

本次未批量安装生物信息学依赖，也未修改用户已有 R/Python 科研环境。

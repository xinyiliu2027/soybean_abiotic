# Soybean abiotic stress single-cell analysis

大豆非生物胁迫单细胞分析代码，整理自图1–图4及 RNA velocity 五份 Word 文档。代码按分析阶段拆分，保留原始算法、阈值、基因列表和主要绘图设置。

本仓库是研究脚本整理版。未包含表达矩阵、Seurat 对象、注释数据库或分析结果，尚未使用原始数据完成端到端复现。各图的输入并不完全相同，不能直接顺序执行整个目录。

## 目录

| 目录 | 内容 |
| --- | --- |
| `scripts/figure1/` | marker 点图、标签转移、样本合并与 UMAP；依赖连续 R 会话 |
| `scripts/figure2/` | muscat pseudobulk 差异分析、四分量扰动指数 PI |
| `scripts/figure3/` | CK 参考轨迹/分期、P2/P3 候选差异响应基因 |
| `scripts/figure4/` | 共识 WGCNA、模块汇总、UMAP/GO、候选 TF regulon |
| `scripts/rna_velocity/` | PEG RNA velocity 可视化 |
| `examples/` | 原文中的草莓花 DoubletFinder 参考示例，不能当作大豆输入直接运行 |
| `legacy/` | 原文图3早期 Cold/PEG 探索片段 |
| `source_text/` | 五份文档的原始正文提取文本，用于追溯；不要直接执行 |
| `docs/` | 脚本索引、来源行号、修改记录、依赖和已知限制 |

完整对应关系见 [脚本索引](docs/SCRIPT_INDEX.md) 和 [来源与修改清单](docs/source_manifest.json)。原始 Word 文件未重复上传。

## 配置与运行

先阅读 [依赖说明](docs/DEPENDENCIES.md) 和 [已知限制](docs/KNOWN_LIMITATIONS.md)。所有命令从仓库根目录执行。

R 脚本使用 `SOYBEAN_PROJECT_ROOT` 替代原服务器的 `/data01/liuxy/Project`，默认 `data`。保持原项目下的相对目录结构，或修改各脚本开头的配置。输出仍遵循原项目的相对布局；重新运行可能更新同名结果/缓存。

```r
# 在 R 中配置；路径应指向你实际存放研究数据的目录。
Sys.setenv(SOYBEAN_PROJECT_ROOT = "D:/soybean_project")
source("scripts/figure2/01_pseudobulk_de.R", encoding = "UTF-8")
```

图1不是独立命令行流程：先载入 `Gm_CK_integrated` 和各个 `Gm_*` Seurat 样本；设置 `SOYBEAN_REFERENCE_QS` 指向已注释参考对象，然后按需要在同一 R 会话运行 01–06。01 的点图不是后续标签转移的输入文件；02–06 共享 `reference`、`hormone.list` 和 `merged_samples` 等对象。不同绘图版本仍分别保留。

图2的两个脚本应在各自的新 R 会话运行。02 需要补全 CK2 的表达对象及外部汇总 DEG 表；01 仅生成逐比较文件，不自动生成 02 所需的汇总表。

图3的 01 与图4的 01 是不同的分期/分析方案，不能把两者的 P1/P2/P3 直接视为相同标签。图3的 02 需要真实阶段 DEG 标准表路径：

```r
Sys.setenv(SOYBEAN_PROJECT_ROOT = "D:/soybean_project")
Sys.setenv(SOYBEAN_STAGE_DEG_FILE = "D:/soybean_project/path/to/DE_standardized_full.csv")
source("scripts/figure3/02_stage_specific_candidates.R", encoding = "UTF-8")
```

图4先运行 01，再按需独立运行 02、03、04；这些脚本通过输出文件衔接。请核对脚本顶部的输入文件、样本列、注释列和生物学重复设置。

Python velocity：

```bash
python -m pip install -r requirements.txt
python scripts/rna_velocity/plot_velocity.py --h5ad data/PEG_velocity_2222_2223.h5ad --metadata data/PEG_metadata_with_umap.csv --output-dir results/rna_velocity
```

CSV 必须含 `cell_id`、`UMAP_1`、`UMAP_2`；细胞 ID 必须符合原文 `CELL数字_N数字` 匹配规则且匹配键唯一。细胞类型优先从 `merged_cell_type`、`cell_type`、`new_cell_type`、`Group1`、`Group` 中选择非空值最多的一列；请确认选中列的实际含义。计算 velocity 需要 h5ad 中的 spliced/unspliced 数据或已计算的 velocity 及相关结构。

## 验证范围

整理时执行静态语法检查，不运行耗时分析，不生成或声称验证生物学结果。具体检查结果见 [验证记录](docs/VALIDATION.md)。依赖版本未由原文提供，因此未编造锁定版本；首次成功运行后应保存 `sessionInfo()` 和 `pip freeze`。

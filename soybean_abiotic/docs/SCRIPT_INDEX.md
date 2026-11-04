# Script index

| Script | Purpose / prerequisites |
| --- | --- |
| [examples/doubletfinder_flower_reference.R](../examples/doubletfinder_flower_reference.R) | 原文草莓花双细胞示例；非大豆分析入口 |
| [scripts/figure1/01_marker_dotplot.R](../scripts/figure1/01_marker_dotplot.R) | 需要已有 Gm_CK_integrated 对象 |
| [scripts/figure1/02_label_transfer.R](../scripts/figure1/02_label_transfer.R) | 需要预载 Gm_* 样本和 reference_renamed.qs |
| [scripts/figure1/03_merge_and_umap.R](../scripts/figure1/03_merge_and_umap.R) | 续接标签转移会话；含原始诊断输出 |
| [scripts/figure1/04_sample_colors_and_merge.R](../scripts/figure1/04_sample_colors_and_merge.R) | 续接上一步，建立 merged_samples |
| [scripts/figure1/05_group_umap.R](../scripts/figure1/05_group_umap.R) | 续接上一步的分组绘图版本 |
| [scripts/figure1/06_combined_umap.R](../scripts/figure1/06_combined_umap.R) | 续接上一步的组合图与坐标修复版本 |
| [scripts/figure2/01_pseudobulk_de.R](../scripts/figure2/01_pseudobulk_de.R) | 独立运行；合并表达RDS和CK RData |
| [scripts/figure2/02_perturbation_index.R](../scripts/figure2/02_perturbation_index.R) | 独立运行；CK2补全RDS和外部汇总DEG表 |
| [legacy/figure3_cold_peg_exploratory.R](../legacy/figure3_cold_peg_exploratory.R) | 原文早期探索片段，止于Step 7标题；不作为主流程 |
| [scripts/figure3/01_ck_reference_states.R](../scripts/figure3/01_ck_reference_states.R) | 独立运行；含CK和三个处理的原始RNA counts |
| [scripts/figure3/02_stage_specific_candidates.R](../scripts/figure3/02_stage_specific_candidates.R) | 独立运行；阶段DEG表和全细胞类型DEG表 |
| [scripts/figure4/01_consensus_wgcna.R](../scripts/figure4/01_consensus_wgcna.R) | 独立运行；CK2补全RDS；旧模块表用于可选比较 |
| [scripts/figure4/02_module_gene_summary.R](../scripts/figure4/02_module_gene_summary.R) | 读取01生成的consensus_network_object.rds |
| [scripts/figure4/03_module_umap_go.R](../scripts/figure4/03_module_umap_go.R) | 读取01生成的模块评分RDS及GO结果 |
| [scripts/figure4/04_tf_regulon_heatmap.R](../scripts/figure4/04_tf_regulon_heatmap.R) | 读取01输出；GENIE3推断候选regulon |
| [scripts/rna_velocity/plot_velocity.py](../scripts/rna_velocity/plot_velocity.py) | 独立运行；h5ad和含UMAP的细胞元数据CSV |

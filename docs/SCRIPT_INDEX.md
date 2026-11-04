# Script index

| Script | Purpose and prerequisites |
| --- | --- |
| [examples/doubletfinder_flower_reference.R](../examples/doubletfinder_flower_reference.R) | Strawberry flower doublet example from the source; not a soybean analysis entry point. |
| [scripts/figure1/01_marker_dotplot.R](../scripts/figure1/01_marker_dotplot.R) | Requires an existing Gm_CK_integrated object. |
| [scripts/figure1/02_label_transfer.R](../scripts/figure1/02_label_transfer.R) | Requires preloaded Gm_* samples and reference_renamed.qs. |
| [scripts/figure1/03_merge_and_umap.R](../scripts/figure1/03_merge_and_umap.R) | Continues the label-transfer session; retains the original diagnostic output. |
| [scripts/figure1/04_sample_colors_and_merge.R](../scripts/figure1/04_sample_colors_and_merge.R) | Continues the preceding session and creates merged_samples. |
| [scripts/figure1/05_group_umap.R](../scripts/figure1/05_group_umap.R) | Grouped plotting variant that continues the preceding session. |
| [scripts/figure1/06_combined_umap.R](../scripts/figure1/06_combined_umap.R) | Continues the preceding session with combined plots and coordinate repair. |
| [scripts/figure2/01_pseudobulk_de.R](../scripts/figure2/01_pseudobulk_de.R) | Standalone script requiring a merged expression RDS and CK RData. |
| [scripts/figure2/02_perturbation_index.R](../scripts/figure2/02_perturbation_index.R) | Standalone script requiring a CK2-complete RDS and an external combined DEG table. |
| [legacy/figure3_cold_peg_exploratory.R](../legacy/figure3_cold_peg_exploratory.R) | Early exploratory fragment ending at the Step 7 heading; not the main workflow. |
| [scripts/figure3/01_ck_reference_states.R](../scripts/figure3/01_ck_reference_states.R) | Standalone script requiring raw RNA counts for CK and three treatments. |
| [scripts/figure3/02_stage_specific_candidates.R](../scripts/figure3/02_stage_specific_candidates.R) | Standalone script requiring stage-level and whole-cell-type DEG tables. |
| [scripts/figure4/01_consensus_wgcna.R](../scripts/figure4/01_consensus_wgcna.R) | Standalone script requiring a CK2-complete RDS; a legacy module table supports optional comparisons. |
| [scripts/figure4/02_module_gene_summary.R](../scripts/figure4/02_module_gene_summary.R) | Reads consensus_network_object.rds produced by script 01. |
| [scripts/figure4/03_module_umap_go.R](../scripts/figure4/03_module_umap_go.R) | Reads the module-scored RDS and GO outputs from script 01. |
| [scripts/figure4/04_tf_regulon_heatmap.R](../scripts/figure4/04_tf_regulon_heatmap.R) | Reads script 01 outputs and infers candidate regulons using GENIE3. |
| [scripts/rna_velocity/plot_velocity.py](../scripts/rna_velocity/plot_velocity.py) | Standalone script requiring h5ad and a cell-metadata CSV with UMAP coordinates. |

# Soybean abiotic stress single-cell analysis

R and Python research scripts organized from four figure-code documents and one RNA velocity document. Scripts retain the supplied analysis parameters, gene lists and statistical methods. Documentation, code comments, plot titles and diagnostic messages are in English.

This repository contains code, not expression data or validated analysis results. It has not been run end to end with the original datasets. The figure workflows require different inputs and are not a single sequential pipeline.

## Repository layout

| Directory | Contents |
| --- | --- |
| `scripts/figure1/` | Marker dot plots, label transfer, sample merging and UMAP; requires a shared R session |
| `scripts/figure2/` | muscat pseudobulk differential expression and a four-component perturbation index |
| `scripts/figure3/` | CK reference trajectory/stages and stage-specific DEG candidates |
| `scripts/figure4/` | Consensus WGCNA, module summaries, UMAP/GO panels and candidate TF regulons |
| `scripts/rna_velocity/` | Python scVelo visualization |
| `examples/` | Original strawberry-flower DoubletFinder example; not a soybean entry point |
| `legacy/` | Earlier exploratory Cold/PEG trajectory fragment |
| `tools/` | Syntax-checking utility |
| `docs/` | Dependencies, script index, provenance, limitations and validation |

All R code uses `.R`; Python code uses `.py`. The package does not include Word documents or `.txt` copies of code. See the [script index](docs/SCRIPT_INDEX.md) and [source manifest](docs/source_manifest.json).

## Configuration

Read [dependencies](docs/DEPENDENCIES.md) and [known limitations](docs/KNOWN_LIMITATIONS.md) before running. Run commands from the repository root.

`SOYBEAN_PROJECT_ROOT` replaces the original server project root and defaults to `data`. Preserve the expected relative directory layout or edit the configuration at the top of each script. Outputs retain the original relative layout; reruns may replace files with the same names.

```r
Sys.setenv(SOYBEAN_PROJECT_ROOT = "D:/soybean_project")
source("scripts/figure2/01_pseudobulk_de.R", encoding = "UTF-8")
```

Figure 1 requires preloaded `Gm_CK_integrated` and `Gm_*` sample objects. Set `SOYBEAN_REFERENCE_QS` to the annotated reference QS file. Script 01 creates marker plots; scripts 02–06 share session objects such as `reference`, `hormone.list` and `merged_samples`. Inspect intermediate objects and choose plotting variants as needed.

Run the two Figure 2 scripts in separate R sessions. The PI script needs a CK2-complete expression RDS and an external combined DEG table. The DE script exports individual comparison tables, not the combined input required by PI.

Figure 3 script 01 and Figure 4 script 01 use different staging methods; their P1/P2/P3 labels are not automatically interchangeable. Figure 3 script 02 requires an actual standardized stage-DEG table:

```r
Sys.setenv(SOYBEAN_PROJECT_ROOT = "D:/soybean_project")
Sys.setenv(SOYBEAN_STAGE_DEG_FILE = "D:/soybean_project/path/to/DE_standardized_full.csv")
source("scripts/figure3/02_stage_specific_candidates.R", encoding = "UTF-8")
```

For Figure 4, run script 01 first, then run scripts 02, 03 and 04 independently as needed. They read files produced by script 01. Verify input paths, sample identifiers, annotation columns and biological-replicate settings.

## RNA velocity

```bash
python -m pip install -r requirements.txt
python scripts/rna_velocity/plot_velocity.py --h5ad data/PEG_velocity_2222_2223.h5ad --metadata data/PEG_metadata_with_umap.csv --output-dir results/rna_velocity
```

The CSV requires `cell_id`, `UMAP_1` and `UMAP_2`. The original matching rule extracts `CELL<digits>_N<digits>`; these keys must be unique. The script selects the most populated column among `merged_cell_type`, `cell_type`, `new_cell_type`, `Group1` and `Group`; verify that the selected column actually represents cell types. Velocity computation requires appropriate spliced/unspliced data or an existing velocity layer and compatible supporting structures.

## Validation

Static syntax checks pass; see the [validation record](docs/VALIDATION.md). This does not establish package compatibility or successful scientific reproduction. Package versions were not supplied, so no verified lockfile is claimed. Save `sessionInfo()` and `pip freeze` after a successful run in the analysis environment.

# Changes and known limitations

## Packaging and localization

Scripts were split by analysis stage, saved as UTF-8 with LF line endings, and given `.R` or `.py` extensions. Project-root paths are configurable. Source ranges, source-document hashes and changes are recorded in `source_manifest.json`. Original extraction files are retained locally outside this upload package.

Comments, documentation, plot titles and diagnostics were translated into English. Original non-English condition aliases and the Endodermis matching pattern retain their exact runtime values through Unicode escapes. The self-contained marker-label dictionary was translated consistently at definitions and lookups. Scientific parameters, gene IDs and numerical expressions were not changed by localization.

Statements such as "complete" or "fixed" in inherited comments describe source-script intent, not proof of end-to-end validation.

## Concrete corrections during initial organization

- Closed the truncated treatment loop at the end of the Figure 2 DE section.
- Corrected the Figure 3 completion message to use `OUT` rather than another fragment's `output_dir`.
- Removed redundant `rm(SA1)` after sample-list creation, guarded temporary-object cleanup, and made multiple DoubletFinder classification columns an explicit error.
- Replaced a historical stage-DEG path with the required `SOYBEAN_STAGE_DEG_FILE` setting.
- Added Python CLI paths, duplicate-key/column-conflict checks, finite-UMAP validation and velocity-graph computation when absent.

## Before analysis

1. Figure 1 lacks complete input-loading code and includes interactive repair/plotting variants. Inspect each stage rather than running the directory unattended. Confirm sample annotations, doublet columns, UMAP prefixes and inclusion of SA samples.
2. Figure 2 DE retains the original fallback to the first two samples when CK1/CK2 names are not identified, plus relaxed pbDS retries. Verify sample identity and filtering choices before use.
3. PI and candidate-gene scripts require external combined/standardized DEG tables whose full production chain was not supplied. Verify columns, comparison direction, thresholds and untested contrasts; missing records are not zero DEGs.
4. CK2-complete and CK-reference RDS files, the annotated QS reference and legacy module tables are not included or generated automatically.
5. Figure 3 retains `BIOLOGICAL_REPLICATES_CONFIRMED = FALSE` and `BIOLOGICAL_DIRECTION_CONFIRMED = FALSE`. Change these only with actual supporting evidence. P1/P2/P3 labels alone do not establish discrete biological states.
6. The TF heatmap represents RNA-inferred candidate regulons, not ATAC/motif activity. The script may download a PlantTFDB list; verify its provenance and gene version.
7. Velocity retains the original identifier extraction, automatic cell-type selection and cell filtering. Duplicate core IDs require a sample-aware key rather than silent deduplication.
8. Scripts may overwrite their existing same-name outputs or caches. Use a separate analysis project/output location for initial runs.

Literature citations, statistical suitability and biological conclusions have not undergone an independent scientific audit in this packaging task.

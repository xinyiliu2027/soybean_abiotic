# Set these environment variables before launching the app when necessary.
config <- list(
  rds = Sys.getenv("SOYBEAN_ATLAS_RDS", "/data01/liuxy/Project/CK2_added_v3_Results/seurat_obj_with_CK2_v3.rds"),
  assay = "RNA",
  all_treatments = TRUE,
  treatments = c("CK", "Heat", "Cold", "PEG"),
  sample_group_map = setNames(character(), character()),
  reduction = NULL,
  group_col = NULL,
  ck_values = c("CK", "CK1", "CK2", "CK_1", "CK_2", "Control", "Control1", "Control2"),
  # The input is a multi-treatment merged object. Never assume every cell is CK.
  # Set group_col / ck_values explicitly if CK cannot be identified in metadata.
  ck_only_object = tolower(Sys.getenv("SOYBEAN_CK_ONLY", "false")) == "true",
  celltype_col = NULL,
  df_column = NULL,
  feature_map = setNames(character(), character()),
  default_genes = c("Glyma.05G103600", "Glyma.15G040200", "Glyma.12G234200", "Glyma.08G163500", "Glyma.10G022500"),
  raster_dpi = 600L,
  seed = 42L
)

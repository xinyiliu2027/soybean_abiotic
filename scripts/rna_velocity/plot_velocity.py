# Source: rna-velocity.docx, extracted text lines 1-251
# Standalone script requiring h5ad and a cell-metadata CSV with UMAP coordinates.
# See README.md and docs/KNOWN_LIMITATIONS.md before running.

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import argparse
from pathlib import Path
import re

import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import matplotlib
matplotlib.use("Agg")  # Render plots without a graphical display.
import matplotlib.pyplot as plt
from matplotlib import cm

# ============================================================
# 0. Configure paths.
# ============================================================
parser = argparse.ArgumentParser(description="Plot soybean PEG RNA velocity with external UMAP metadata")
parser.add_argument("--h5ad", type=Path, required=True)
parser.add_argument("--metadata", type=Path, required=True)
parser.add_argument("--output-dir", type=Path, default=Path("results/rna_velocity"))
args = parser.parse_args()
h5ad_path = args.h5ad
meta_csv = args.metadata
output_dir = args.output_dir
output_dir.mkdir(parents=True, exist_ok=True)

if not h5ad_path.is_file():
    raise FileNotFoundError(f"File not found: {h5ad_path}. Check the h5ad path.")

if not meta_csv.is_file():
    raise FileNotFoundError(f"Metadata file not found: {meta_csv}. Check the path.")

# ============================================================
# 1. Read h5ad and metadata.
# ============================================================
print(f"Reading AnnData from {h5ad_path} ...")
adata = sc.read(h5ad_path)
print("adata shape:", adata.shape)
print("Example adata.obs indices:", adata.obs_names[:5].tolist())
print()

print(f"Reading metadata with UMAP from {meta_csv} ...")
meta = pd.read_csv(meta_csv)
print("Metadata shape:", meta.shape)
print("First 10 metadata columns:", meta.columns.tolist()[:10])

print()

required_cols = ["cell_id", "UMAP_1", "UMAP_2"]
for col in required_cols:
    if col not in meta.columns:
        raise ValueError(f"Required metadata column missing: {col}. Check {meta_csv}")

# ============================================================
# 2. Extract core_id (CELLxxx_Ny).
# ============================================================
def extract_core_id(s: str):
    "Extract a substring of the form CELL<digits>_N<digits>."
    if not isinstance(s, str):
        return None
    m = re.search(r"(CELL\d+_N\d+)", s)
    return m.group(1) if m else None

print("Extracting and aligning core_id in AnnData and metadata...")

# AnnData identifiers.
adata.obs["core_id"] = [extract_core_id(str(x)) for x in adata.obs_names]
print("Example adata.obs core_id values:", adata.obs["core_id"].head().tolist())
print("Nonmissing adata.obs core_id count:", adata.obs["core_id"].notna().sum())

adata = adata[adata.obs["core_id"].notna()].copy()
print("adata.n_obs after removing missing core_id:", adata.n_obs)
print()

# Metadata identifiers.
meta["cell_id"] = meta["cell_id"].astype(str)
meta["core_id"] = meta["cell_id"].apply(extract_core_id)
print("Example metadata core_id values:", meta["core_id"].head().tolist())
print("Nonmissing metadata core_id count:", meta["core_id"].notna().sum())

meta = meta[meta["core_id"].notna()].copy()
if meta["core_id"].duplicated().any() or adata.obs["core_id"].duplicated().any():
    raise ValueError("core_id is not unique; use a sample-aware key before merging metadata.")
meta = meta.set_index("core_id")
print("Metadata rows after filtering and uniqueness checks:", meta.shape[0])
print()

# Find shared identifiers.
common_core = np.intersect1d(
    adata.obs["core_id"].values,
    meta.index.values,
)
print("Shared core_id count:", len(common_core))
if len(common_core) == 0:
    raise ValueError("No shared core_id values; check identifier formats.")

adata = adata[adata.obs["core_id"].isin(common_core)].copy()
print("adata.n_obs after retaining shared core_id:", adata.n_obs)

overlap = adata.obs.columns.intersection(meta.columns)
if len(overlap):
    raise ValueError(f"Metadata columns overlap with existing obs: {list(overlap)}. Resolve explicitly before merging.")
adata.obs = adata.obs.join(meta, on="core_id", how="left")
print("adata.obs shape after metadata merge:", adata.obs.shape)
print("Example merged obs columns:", adata.obs.columns[:10].tolist())
print()

# ============================================================
# 3. UMAP coordinates.
# ============================================================
if ("UMAP_1" not in adata.obs.columns) or ("UMAP_2" not in adata.obs.columns):
    raise ValueError("UMAP_1/UMAP_2 missing from adata.obs after merging.")

adata.obsm["X_umap"] = adata.obs[["UMAP_1", "UMAP_2"]].to_numpy()
if not np.isfinite(adata.obsm["X_umap"]).all():
    raise ValueError("UMAP coordinates contain missing or non-finite values.")
print("UMAP coordinates assigned to adata.obsm['X_umap'].")
print()

# ============================================================
# 4. Select a populated cell-type column.
# ============================================================
celltype_col = None
candidates = [
    "merged_cell_type",
    "cell_type",
    "new_cell_type",
    "Group1",
    "Group",
]

best_col = None
best_nnz = -1

for cand in candidates:
    if cand in adata.obs.columns:
        nnz = adata.obs[cand].notna().sum()
        print(f"Candidate cell-type column {cand}: nonmissing count = {nnz}")
        if nnz > best_nnz:
            best_nnz = nnz
            best_col = cand

if best_col is None or best_nnz == 0:
    print("WARNING: No populated candidate column; labeling all cells as unknown.")
    adata.obs["celltype"] = "unknown"
else:
    celltype_col = best_col
    print(f"Using '{celltype_col}' for celltype (nonmissing count = {best_nnz}).")
    adata.obs["celltype"] = adata.obs[celltype_col]

# Convert NaN and the string "nan" to unknown, then cast to categorical.
adata.obs["celltype"] = (
    adata.obs["celltype"]
    .astype(str)
    .replace({"nan": np.nan})
    .fillna("unknown")
    .astype("category")
)
cats = adata.obs["celltype"].cat.categories
print("Final cell-type categories:", list(cats))
print()

# ============================================================
# 5. Compute velocity if necessary.
# ============================================================
scv.settings.verbosity = 3
scv.settings.presenter_view = False
scv.settings.figdir = str(output_dir)
os.makedirs(scv.settings.figdir, exist_ok=True)
scv.set_figure_params("scvelo", dpi=150, frameon=False)

if "velocity" not in adata.layers.keys():
    print("No adata.layers['velocity']; computing velocity now (this may take time)...")
    # Basic scVelo workflow.
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    scv.tl.velocity(adata, mode="stochastic")
    scv.tl.velocity_graph(adata)
    # Save the processed h5ad after computing velocity.
    adata.write(output_dir / "PEG_velocity_2222_2223_processed.h5ad", compression="gzip")
    print("Velocity computed; processed h5ad saved in the output directory.")
else:
    print("Existing adata.layers['velocity'] detected; using it for plotting.")

# ============================================================
# 6. Assign colors to cell types.
# ============================================================
if "velocity_graph" not in adata.uns:
    scv.tl.velocity_graph(adata)

cats = adata.obs["celltype"].cat.categories
cmap = cm.get_cmap("tab20", len(cats) if len(cats) > 0 else 1)
palette = {cat: cmap(i) for i, cat in enumerate(cats)}

print("Cell types and assigned colors:")
for c in cats:
    print("  ", c, "->", palette[c])
print()

# ============================================================
# 7. Plot and save figures.
# ============================================================
# 1. UMAP colored by cell type.
print("Plotting colored UMAP by celltype ...")
ax = scv.pl.scatter(
    adata,
    basis="umap",
    color="celltype",
    palette=palette,
    legend_loc="right",
    size=8,
    alpha=0.9,
    show=False,
)
ax.figure.suptitle("PEG UMAP (cell types)", fontsize=18)
ax.figure.savefig(
    output_dir / "PEG_celltype_umap_color.png",
    dpi=300,
    bbox_inches="tight",
)
plt.close(ax.figure)

# 2. Velocity streamlines colored by cell type.
print("Plotting velocity stream colored by celltype ...")
ax = scv.pl.velocity_embedding_stream(
    adata,
    basis="umap",
    color="celltype",
    palette=palette,
    legend_loc="right",
    title="PEG RNA velocity (celltype, stream)",
    density=2,
    linewidth=1,
    arrow_size=1.5,
    show=False,
)
ax.figure.savefig(
    output_dir / "PEG_velocity_stream_celltype_color.png",
    dpi=300,
    bbox_inches="tight",
)
plt.close(ax.figure)

# 3. Velocity arrows colored by cell type.
print("Plotting velocity arrows colored by celltype ...")
ax = scv.pl.velocity_embedding(
    adata,
    basis="umap",
    color="celltype",
    palette=palette,
    arrow_size=1.5,
    arrow_length=3,
    dpi=150,
    legend_loc="right",
    title="PEG RNA velocity (celltype, arrows)",
    show=False,
)
ax.figure.savefig(
    output_dir / "PEG_velocity_embedding_celltype_color.png",
    dpi=300,
    bbox_inches="tight",
)
plt.close(ax.figure)

print(f"All done. Plots saved in {output_dir}")

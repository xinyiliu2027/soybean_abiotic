# Source: rna-velocity.docx, extracted text lines 1-251
# 独立运行；h5ad和含UMAP的细胞元数据CSV
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
matplotlib.use("Agg")  # 后台画图
import matplotlib.pyplot as plt
from matplotlib import cm

# ============================================================
# 0. 路径设置
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
    raise FileNotFoundError(f"找不到 {h5ad_path}，请确认 h5ad 文件位置。")

if not meta_csv.is_file():
    raise FileNotFoundError(f"找不到 meta 文件 {meta_csv}，请确认路径。")

# ============================================================
# 1. 读 h5ad 和 meta
# ============================================================
print(f"Reading AnnData from {h5ad_path} ...")
adata = sc.read(h5ad_path)
print("adata shape:", adata.shape)
print("adata.obs index 示例:", adata.obs_names[:5].tolist())
print()

print(f"Reading metadata with UMAP from {meta_csv} ...")
meta = pd.read_csv(meta_csv)
print("meta 形状:", meta.shape)
print("meta 列名前 10:", meta.columns.tolist()[:10])

print()

required_cols = ["cell_id", "UMAP_1", "UMAP_2"]
for col in required_cols:
    if col not in meta.columns:
        raise ValueError(f"meta 中缺少必须列：{col} ，请检查 {meta_csv}")

# ============================================================
# 2. 提取 core_id（CELLxxx_Ny）
# ============================================================
def extract_core_id(s: str):
    """提取 'CELL数字_N数字' 子串"""
    if not isinstance(s, str):
        return None
    m = re.search(r"(CELL\d+_N\d+)", s)
    return m.group(1) if m else None

print("提取 adata 和 meta 的 core_id 并对齐 ...")

# ---- adata 这边 ----
adata.obs["core_id"] = [extract_core_id(str(x)) for x in adata.obs_names]
print("adata.obs core_id 示例:", adata.obs["core_id"].head().tolist())
print("adata.obs 非 NA core_id 数量:", adata.obs["core_id"].notna().sum())

adata = adata[adata.obs["core_id"].notna()].copy()
print("过滤 NA core_id 后 adata.n_obs:", adata.n_obs)
print()

# ---- meta 这边 ----
meta["cell_id"] = meta["cell_id"].astype(str)
meta["core_id"] = meta["cell_id"].apply(extract_core_id)
print("meta core_id 示例:", meta["core_id"].head().tolist())
print("meta 非 NA core_id 数量:", meta["core_id"].notna().sum())

meta = meta[meta["core_id"].notna()].copy()
if meta["core_id"].duplicated().any() or adata.obs["core_id"].duplicated().any():
    raise ValueError("core_id is not unique; use a sample-aware key before merging metadata.")
meta = meta.set_index("core_id")
print("meta 去 NA + 去重后行数:", meta.shape[0])
print()

# ---- 找交集 ----
common_core = np.intersect1d(
    adata.obs["core_id"].values,
    meta.index.values,
)
print("两边共同的 core_id 数量:", len(common_core))
if len(common_core) == 0:
    raise ValueError("没有任何共同 core_id，ID 规则对不上，请检查。")

adata = adata[adata.obs["core_id"].isin(common_core)].copy()
print("过滤到共同 core_id 后 adata.n_obs:", adata.n_obs)

overlap = adata.obs.columns.intersection(meta.columns)
if len(overlap):
    raise ValueError(f"Metadata columns overlap with existing obs: {list(overlap)}. Resolve explicitly before merging.")
adata.obs = adata.obs.join(meta, on="core_id", how="left")
print("合并 meta 后 adata.obs 形状:", adata.obs.shape)
print("合并后 obs 示例列:", adata.obs.columns[:10].tolist())
print()

# ============================================================
# 3. UMAP 坐标
# ============================================================
if ("UMAP_1" not in adata.obs.columns) or ("UMAP_2" not in adata.obs.columns):
    raise ValueError("合并后在 adata.obs 中找不到 'UMAP_1'/'UMAP_2'。")

adata.obsm["X_umap"] = adata.obs[["UMAP_1", "UMAP_2"]].to_numpy()
if not np.isfinite(adata.obsm["X_umap"]).all():
    raise ValueError("UMAP coordinates contain missing or non-finite values.")
print("UMAP 坐标已写入 adata.obsm['X_umap']。")
print()

# ============================================================
# 4. 选择一个“真的有值”的 celltype 列
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
        print(f"候选细胞类型列 {cand}: 非 NA 数量 = {nnz}")
        if nnz > best_nnz:
            best_nnz = nnz
            best_col = cand

if best_col is None or best_nnz == 0:
    print("WARNING: 所有候选列几乎都是 NA，全部标记为 'unknown'")
    adata.obs["celltype"] = "unknown"
else:
    celltype_col = best_col
    print(f"选择 '{celltype_col}' 作为 celltype（非 NA 个数 = {best_nnz}）。")
    adata.obs["celltype"] = adata.obs[celltype_col]

# 把 NaN 和字符串 "nan" 都归为 unknown，再转为分类
adata.obs["celltype"] = (
    adata.obs["celltype"]
    .astype(str)
    .replace({"nan": np.nan})
    .fillna("unknown")
    .astype("category")
)
cats = adata.obs["celltype"].cat.categories
print("最终 celltype categories:", list(cats))
print()

# ============================================================
# 5. 如有必要，重新计算 velocity
# ============================================================
scv.settings.verbosity = 3
scv.settings.presenter_view = False
scv.settings.figdir = str(output_dir)
os.makedirs(scv.settings.figdir, exist_ok=True)
scv.set_figure_params("scvelo", dpi=150, frameon=False)

if "velocity" not in adata.layers.keys():
    print("没有发现 adata.layers['velocity']，现在重新计算 velocity（会稍微久一点）...")
    # 基本 scvelo pipeline
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    scv.tl.velocity(adata, mode="stochastic")
    scv.tl.velocity_graph(adata)
    # 计算完可以顺便保存一份 processed 的 h5ad
    adata.write(output_dir / "PEG_velocity_2222_2223_processed.h5ad", compression="gzip")
    print("velocity 计算完成，已保存到 PEG_velocity_2222_2223_processed.h5ad")
else:
    print("检测到已有 adata.layers['velocity']，直接用于作图。")

# ============================================================
# 6. 为每个 celltype 生成颜色 palette
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
# 7. 画图并保存
# ============================================================
# 1) 纯 UMAP（celltype 上色）
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

# 2) 流线 velocity（celltype 上色）
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

# 3) 箭头 velocity（celltype 上色）
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

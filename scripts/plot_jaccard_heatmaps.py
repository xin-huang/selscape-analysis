# Copyright 2026 Xin Huang and Simon Chen
#
# GNU General Public License v3.0
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, please see
#
#    https://www.gnu.org/licenses/gpl-3.0.en.html


import sys
from pathlib import Path
 
sys.path.insert(0, snakemake.scriptdir)
from jaccard_common import pooled_gene_sets, jaccard_matrix
 
import matplotlib
 
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
 
log_fh = open(snakemake.log[0], "w") if snakemake.log else sys.stdout
sys.stderr = log_fh
sys.stdout = log_fh
 
def parse_positive_path(rel_path):
    parts = rel_path.split("/")
    if len(parts) == 7 and parts[0] == "selscan" and parts[3] in ("1pop", "2pop"):
        method = parts[5].split("_", 1)[0]  # folder is "{method}_{maf}", e.g. "ihs_0.05"
        return dict(category="positive_selection", tool="selscan", method=method,
                    dataset=parts[2], unit=parts[4])
    if len(parts) == 8 and parts[0] == "scikit-allel" and parts[3] in ("1pop", "2pop"):
        return dict(category="positive_selection", tool="scikit-allel", method=parts[5],
                    dataset=parts[2], unit=parts[4])
    return None
 
 
def parse_balancing_path(rel_path):
    parts = rel_path.split("/")
    if len(parts) == 7 and parts[0] == "scikit-allel":
        return dict(category="balancing_selection", tool="scikit-allel", method=parts[3],
                    dataset=parts[2], unit=parts[4])
    if len(parts) == 6 and parts[0] == "betascan":
        return dict(category="balancing_selection", tool="betascan", method="betascan_b1",
                    dataset=parts[2], unit=parts[3])
    return None
 
 
def plot_heatmap(matrix, title, output_path):
    fig, ax = plt.subplots(figsize=(0.6 * len(matrix) + 2, 0.6 * len(matrix) + 2), dpi=300)
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("#e8e8e8")
    masked = np.ma.masked_invalid(matrix.values)
    im = ax.imshow(masked, vmin=0, vmax=1, cmap=cmap)
 
    ax.set_xticks(range(len(matrix.columns)))
    ax.set_xticklabels(matrix.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(matrix.index)))
    ax.set_yticklabels(matrix.index)
 
    for i in range(len(matrix.index)):
        for j in range(len(matrix.columns)):
            val = matrix.values[i, j]
            if np.isnan(val):
                continue
            color = "white" if val < 0.5 else "black"
            ax.text(j, i, f"{val:.2f}", ha="center", va="center", color=color, fontsize=8)
 
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Jaccard similarity")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)
 
 
# --- main ---
 
datasets = snakemake.params.datasets
assert datasets, "config['datasets_to_compare'] is empty -- nothing to compare."
 
outdir = Path(snakemake.output.outdir)
outdir.mkdir(parents=True, exist_ok=True)
 
sources = [
    ("results/positive_selection", parse_positive_path),
    ("results/balancing_selection", parse_balancing_path),
]
 
for root_dir, parse_path in sources:
    pooled = pooled_gene_sets(Path(root_dir), parse_path)
    if not pooled:
        print(f"no outlier.genes files found under {root_dir}, skipping", file=log_fh)
        continue
 
    combos = set()
    for category, tool, method, dataset in pooled:
        combos.add((category, tool, method))
    combos = sorted(combos)
 
    for category, tool, method in combos:
        matrix = jaccard_matrix(pooled, category, tool, method, datasets)
        stem = f"{category}.{tool}.{method}"
        matrix.to_csv(outdir / f"{stem}.jaccard_matrix.tsv", sep="\t")
        title = f"{category.replace('_', ' ')} - {tool} ({method})"
        plot_heatmap(matrix, title=title,
                     output_path=outdir / f"{stem}.jaccard_heatmap.svg")
        print(f"wrote {stem}.jaccard_heatmap.svg ({len(matrix)}x{len(matrix)})", file=log_fh)

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
from jaccard_common import (
    pooled_gene_sets,
    jaccard_matrix,
    population_gene_sets,
    population_jaccard_table,
)
from compare_common import parse_positive_path, parse_balancing_path
 
import matplotlib
 
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
 
log_fh = open(snakemake.log[0], "w") if snakemake.log else sys.stdout
sys.stderr = log_fh
sys.stdout = log_fh
 
 
def plot_heatmap(matrix, title, output_path):
    n_rows = len(matrix.index)
    n_cols = len(matrix.columns)
    fig, ax = plt.subplots(figsize=(0.6 * n_cols + 3, 0.35 * n_rows + 2), dpi=300)
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("#e8e8e8")
    masked = np.ma.masked_invalid(matrix.values)
    im = ax.imshow(masked, vmin=0, vmax=1, cmap=cmap, aspect="auto")
 
    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(matrix.columns, rotation=45, ha="right")
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(matrix.index)
 
    for i in range(n_rows):
        for j in range(n_cols):
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
 
datasets = snakemake.params.datasets
assert datasets, "config['datasets_to_compare'] is empty -- nothing to compare."
 
populations = snakemake.params.populations
population_datasets = snakemake.params.population_datasets
 
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
 
    if population_datasets and populations:
        gene_sets = population_gene_sets(Path(root_dir), parse_path)
        for category, tool, method in combos:
            table = population_jaccard_table(
                gene_sets, category, tool, method, populations, population_datasets
            )
            stem = f"{category}.{tool}.{method}.population_level"
            table.to_csv(outdir / f"{stem}.jaccard_matrix.tsv", sep="\t")
            title = f"{category.replace('_', ' ')} - {tool} ({method}) - population level"
            plot_heatmap(table, title=title,
                         output_path=outdir / f"{stem}.jaccard_heatmap.svg")
            print(f"wrote {stem}.jaccard_heatmap.svg ({len(table)}x{len(table.columns)})", file=log_fh)

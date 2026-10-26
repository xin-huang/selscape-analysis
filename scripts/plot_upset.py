# Copyright 2025 Xin Huang and Simon Chen
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

"""UpSet plot: candidate genes shared between population groups (XP and within-population)."""

import os
import sys
 
import matplotlib
 
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from upsetplot import UpSet, from_contents
 
log_fh = open(snakemake.log[0], "w")
sys.stderr = log_fh
sys.stdout = log_fh
 
population_groups = snakemake.params.population_groups or {}
max_intersections = int(snakemake.params.max_intersections)
plot_title = snakemake.params.title
 
BAR = "#0072B2"
FALLBACK_COLORS = ["#E69F00", "#56B4E9", "#009E73", "#0072B2", "#CC79A7", "#D55E00"]
 
 
def no_results(message):
    with open(snakemake.output.table, "w") as out:
        out.write("groups\tdegree\tn_genes\n")
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.text(0.5, 0.55, "No results", ha="center", va="center", fontsize=16, color="#666")
    ax.text(0.5, 0.40, message, ha="center", va="center", fontsize=9, color="#999")
    ax.axis("off")
    plt.savefig(snakemake.output.plot, bbox_inches="tight")
    plt.close()
    print(f"no results: {message}")
    log_fh.close()
    sys.exit(0)
 
 
# (gene file, population) pairs. Population IDs come from the rule, never from file names,
# so IDs containing "." or "_" stay intact.
genes_input = snakemake.input.genes
if isinstance(genes_input, str) and os.path.isdir(genes_input):
    # XP: directory with one gene list per (pair, focal population)
    entries = [(os.path.join(genes_input, name), focal) for name, focal, _reference in snakemake.params.focal_files]
else:
    # within-population: one gene list per population
    gene_files = [genes_input] if isinstance(genes_input, str) else list(genes_input)
    populations = list(snakemake.params.populations)
    if len(gene_files) != len(populations):
        raise ValueError("number of gene files and populations differ")
    entries = list(zip(gene_files, populations))
 
missing = [path for path, _ in entries if not os.path.exists(path)]
if missing:
    raise FileNotFoundError("expected gene lists not found:\n" + "\n".join(missing))
 
group_of = {
    population: group
    for group, cfg in population_groups.items()
    for population in (cfg or {}).get("populations", [])
}
 
population_genes = {}
for path, population in entries:
    with open(path) as handle:
        genes = {g for g in (line.strip() for line in handle) if g and g != "Gene"}
    population_genes.setdefault(population, set()).update(genes)
 
# one set per group; if the dataset spans a single group (one great ape species), one set per population
by_group = len({group_of.get(p) for p, genes in population_genes.items() if genes} - {None}) >= 2
sets, ungrouped = {}, set()
for population, genes in population_genes.items():
    key = group_of.get(population) if by_group else population
    if key is None:
        key = population
        ungrouped.add(population)
    sets.setdefault(key, set()).update(genes)
 
if ungrouped:
    print("not in population_groups, shown as own category:", ", ".join(sorted(ungrouped)))
 
order = list(population_groups) if by_group else [
    p for g in population_groups for p in (population_groups[g] or {}).get("populations", [])
]
groups = [name for name in order if sets.get(name)] + sorted(n for n in sets if sets[n] and n not in order)
if len(groups) < 2:
    no_results("fewer than two categories have candidates")
 
memberships = from_contents({group: sets[group] for group in groups})
intersections = memberships.index.value_counts()
 
with open(snakemake.output.table, "w") as out:
    out.write("groups\tdegree\tn_genes\n")
    for pattern, count in intersections.items():
        members = [group for group, present in zip(groups, pattern) if present]
        out.write(f"{'&'.join(members)}\t{len(members)}\t{count}\n")
 
upset = UpSet(
    memberships,
    subset_size="count",
    sort_by="cardinality",
    sort_categories_by=None,
    max_subset_rank=max_intersections,
    show_counts=True,
    facecolor=BAR,
    element_size=None,
)
for i, group in enumerate(groups):
    color_group = group if by_group else group_of.get(group)
    color = (population_groups.get(color_group) or {}).get("color") or FALLBACK_COLORS[i % len(FALLBACK_COLORS)]
    upset.style_categories(group, bar_facecolor=color)
 
fig = plt.figure(figsize=(max(6.0, 0.55 * min(len(intersections), max_intersections) + 3.0), 4.5), dpi=200)
axes = upset.plot(fig=fig)
axes["intersections"].grid(False)
axes["totals"].grid(False)
axes["totals"].set_xlabel("Total outlier genes")
fig.suptitle(plot_title, fontsize=11, fontweight="bold")
plt.savefig(snakemake.output.plot, bbox_inches="tight")
plt.close()
 
print(f"{len(intersections)} intersections, showing at most {max_intersections}")
log_fh.close()

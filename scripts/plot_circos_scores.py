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



import pandas as pd
import numpy as np
from pycirclize import Circos
from pycirclize.utils import load_eukaryote_example_dataset
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os

# Get parameters
population = snakemake.params.population
plot_type = snakemake.params.plot_type

# load score file
def load_score_file(score_file, score_column):
    df = pd.read_csv(score_file, sep="\t")
    df["chr"] = "chr" + df["CHR"].astype(str)
    df["pos"] = df["BP"].astype(int) if "BP" in df.columns else df["window_start"].astype(int)
    df["score"] = df[score_column].astype(float)
    df = df[~df["chr"].isin(["chrX", "chrY"])]
    df = df.dropna(subset=["score"])
    return df[["chr", "pos", "score"]]

# Track configurations
TRACK_CONFIGS = {
    "positive_selection": [
        ("iHS", "ihs_scores", "normalized_ihs", (61, 75), "#1f77b4"),
        ("nSL", "nsl_scores", "normalized_nsl", (48, 61), "#ff7f0e"),
        ("Moving Tajima's D (−)", "mtjd_scores", "tajima_d", (35, 48), "#2ca02c"),
        ("Windowed Tajima's D (−)", "wtjd_scores", "tajima_d", (22, 35), "#d62728")
    ],
    "balancing_selection": [
        ("B1", "b1_scores", "B1", (48, 75), "#1f77b4"),
        ("Moving Tajima's D (+)", "mtjd_bal_scores", "tajima_d", (35, 48), "#2ca02c"),
        ("Windowed Tajima's D (+)", "wtjd_bal_scores", "tajima_d", (22, 35), "#d62728")
    ]
}

# Load all data
track_data = []
for name, input_key, score_col, r_range, color in TRACK_CONFIGS[plot_type]:
    score_file = getattr(snakemake.input, input_key)
    df = load_score_file(score_file, score_col)
    
    vmin = df["score"].min()
    vmax = df["score"].max()
    
    track_data.append((name, df, r_range, color, vmin, vmax))

# Setup circos
chr_bed_file, cytoband_file, _ = load_eukaryote_example_dataset("hg38")
chr_bed_df = pd.read_csv(chr_bed_file, sep="\t", header=None)
autosomes = [f"chr{i}" for i in range(1, 23)]
chr_bed_df = chr_bed_df[chr_bed_df[0].isin(autosomes)]

temp_dir = os.path.dirname(snakemake.output.plot)
os.makedirs(temp_dir, exist_ok=True)
filtered_chr_bed = os.path.join(temp_dir, f"{population}_{plot_type}_chr.bed")
chr_bed_df.to_csv(filtered_chr_bed, sep="\t", header=False, index=False)

circos = Circos.initialize_from_bed(filtered_chr_bed, space=2)
circos.add_cytoband_tracks((81, 85), cytoband_file)

# Plot each chromosome
for sector in circos.sectors:
    # Chromosome axis
    axis_track = sector.add_track((81, 85), r_pad_ratio=0.0)
    axis_track.axis(lw=0.6)
    sector.text(sector.name, r=88, size=8)
    
    # Plot each statistic track
    for stat_name, df, (r_inner, r_outer), color, vmin, vmax in track_data:
        track = sector.add_track((r_inner, r_outer), r_pad_ratio=0.0)
        track.axis(lw=0.3)
        
        chr_data = df[df["chr"] == sector.name].sort_values("pos")
        
        if len(chr_data) == 0:
            continue
        
        track.line(
            chr_data["pos"].values,
            chr_data["score"].values,
            color=color,
            lw=0.3,
            vmin=vmin,
            vmax=vmax
        )

# Title
circos.text(population, size=12, r=0, weight="bold")

# Legend with score ranges
fig = circos.plotfig()
legend_elements = [
    mpatches.Patch(
        color=color,
        label=f"{name} range=[{vmin:.2f}, {vmax:.2f}]"
    )
    for name, df, _, color, vmin, vmax in track_data
]

fig.legend(
    handles=legend_elements,
    loc="center left",
    bbox_to_anchor=(1.02, 0.5),
    fontsize=9,
    frameon=True,
    title="Statistics"
)

fig.savefig(snakemake.output.plot, dpi=300, bbox_inches="tight")
plt.close()

if os.path.exists(filtered_chr_bed):
    os.remove(filtered_chr_bed)

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
from pycirclize import Circos
from pycirclize.utils import load_eukaryote_example_dataset
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os

def load_bed_without_sex_chr(bed_file):
    """Load BED file and remove sex chromosomes"""
    df = pd.read_csv(bed_file, sep="\t", header=None, dtype={0: str})
    df.columns = ["chr", "start", "end", "gene"]
    df["chr"] = df["chr"].astype(str).str.strip()
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    df = df[~df["chr"].isin(["chrX", "chrY"])]
    return df

# Track configurations for different plot types
TRACK_CONFIGS = {
    "positive_selection": [
        {
            "name": "iHS",
            "bed_input": "ihs_bed",
            "r_range": (61, 75),
            "color": "#1f77b4"  # Blue
        },
        {
            "name": "nSL",
            "bed_input": "nsl_bed",
            "r_range": (48, 61),
            "color": "#ff7f0e"  # Orange
        },
        {
            "name": "Moving Tajima's D (−)",
            "bed_input": "mtjd_bed",
            "r_range": (35, 48),
            "color": "#2ca02c"  # Green
        },
        {
            "name": "Windowed Tajima's D (−)",
            "bed_input": "wtjd_bed",
            "r_range": (22, 35),
            "color": "#d62728"  # Red
        }
    ],
    "balancing_selection": [
        {
            "name": "B1",
            "bed_input": "b1_bed",
            "r_range": (48, 75),
            "color": "#1f77b4"  # Blue
        },
        {
            "name": "Moving Tajima's D (+)",
            "bed_input": "mtjd_bal_bed",
            "r_range": (35, 48),
            "color": "#2ca02c"  # Green
        },
        {
            "name": "Windowed Tajima's D (+)",
            "bed_input": "wtjd_bal_bed",
            "r_range": (22, 35),
            "color": "#d62728"  # Red
        }
    ]
}

# Get parameters
population = snakemake.params.population
plot_type = snakemake.params.plot_type

# Get track configuration for this plot type
track_specs = TRACK_CONFIGS[plot_type]

# Load data and build track configs
track_configs = []
for spec in track_specs:
    bed_file = getattr(snakemake.input, spec["bed_input"])
    df = load_bed_without_sex_chr(bed_file)
    track_configs.append((
        spec["name"],
        df,
        spec["r_range"],
        spec["color"]
    ))

# Initialize circos sectors from chromosome BED
chr_bed_file, cytoband_file, _ = load_eukaryote_example_dataset("hg38")

# Read and filter to only autosomes (chr1-22)
chr_bed_df = pd.read_csv(chr_bed_file, sep="\t", header=None)
autosomes = [f"chr{i}" for i in range(1, 23)]
chr_bed_df = chr_bed_df[chr_bed_df[0].isin(autosomes)]

# Create temporary filtered file
temp_dir = os.path.dirname(snakemake.output.plot)
os.makedirs(temp_dir, exist_ok=True)
filtered_chr_bed = os.path.join(temp_dir, f"{population}_{plot_type}_chr_filtered.bed")
chr_bed_df.to_csv(filtered_chr_bed, sep="\t", header=False, index=False)

# Initialize circos with filtered chromosomes
circos = Circos.initialize_from_bed(filtered_chr_bed, space=2)

# Add cytoband tracks
circos.add_cytoband_tracks((81, 85), cytoband_file)

# Plot tracks
for sector in circos.sectors:
    # Outer chromosome axis track
    axis_track = sector.add_track((81, 85), r_pad_ratio=0.0)
    axis_track.axis(lw=0.6)
    sector.text(sector.name, r=88, size=8)

    # Add tracks for each statistic
    for stat_name, df, (r_inner, r_outer), color in track_configs:
        gene_track = sector.add_track((r_inner, r_outer), r_pad_ratio=0.0)
        gene_track.axis(lw=0.3)

        # Filter genes on this chromosome
        sub = df[df["chr"] == sector.name]

        # Draw gene blocks
        for s, e in zip(sub["start"].to_numpy(), sub["end"].to_numpy()):
            gene_track.rect(float(s), float(e), fc=color, ec="none", lw=1)

# Add only population name as title
circos.text(population, size=12, r=0, weight="bold")

# Create figure and add legend
fig = circos.plotfig()

# Create legend with unique gene counts
legend_elements = []
for stat_name, df, _, color in track_configs:
    legend_elements.append(
        mpatches.Patch(color=color, label=f"{stat_name}")
#        mpatches.Patch(color=color, label=f"{name} (n={len(df)})
    )

# Position legend outside the plot area
legend = fig.legend(
    handles=legend_elements,
    loc="center left",
    bbox_to_anchor=(1.02, 0.5),
    fontsize=10,
    frameon=True,
    fancybox=True,
    shadow=True,
    title="Statistics",
    title_fontsize=10
)

# Save figure with extra space for legend
fig.savefig(snakemake.output.plot, dpi=300, bbox_inches="tight")
plt.close()

# Clean up temporary file
if os.path.exists(filtered_chr_bed):
    os.remove(filtered_chr_bed)

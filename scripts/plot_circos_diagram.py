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

# Function to load and filter BED data
def load_bed_without_sex_chr(bed_file):
    """Load BED file and remove sex chromosomes"""
    df = pd.read_csv(bed_file, sep="\t", header=None, dtype={0: str})
    df.columns = ["chr", "start", "end", "gene"]
    df["chr"] = df["chr"].astype(str).str.strip()
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    # Remove sex chromosomes
    df = df[~df["chr"].isin(["chrX", "chrY"])]
    return df

# Read gene BED files from each statistic
ihs_df = load_bed_without_sex_chr(snakemake.input.ihs_bed)
nsl_df = load_bed_without_sex_chr(snakemake.input.nsl_bed)
mtjd_df = load_bed_without_sex_chr(snakemake.input.mtjd_bed)
wtjd_df = load_bed_without_sex_chr(snakemake.input.wtjd_bed)
b1_df = load_bed_without_sex_chr(snakemake.input.b1_bed)

# Initialize circos sectors from chromosome BED
chr_bed_file, cytoband_file, _ = load_eukaryote_example_dataset("hg38")

# Read and filter to only autosomes (chr1-22)
chr_bed_df = pd.read_csv(chr_bed_file, sep="\t", header=None)
autosomes = [f"chr{i}" for i in range(1, 23)]
chr_bed_df = chr_bed_df[chr_bed_df[0].isin(autosomes)]

# Create temporary filtered file
import os
temp_dir = os.path.dirname(snakemake.output.plot)
os.makedirs(temp_dir, exist_ok=True)
filtered_chr_bed = os.path.join(temp_dir, f"{snakemake.params.population}_chr_filtered.bed")
chr_bed_df.to_csv(filtered_chr_bed, sep="\t", header=False, index=False)

# Initialize circos with filtered chromosomes
circos = Circos.initialize_from_bed(filtered_chr_bed, space=2)

# Add cytoband tracks
circos.add_cytoband_tracks((81, 85), cytoband_file)

# Track positions (inner to outer) - each track represents one statistic
# 5 tracks now: positive selection (iHS, nSL, 2x Tajima's D) + balancing selection (B1)
track_configs = [
    ("iHS", ihs_df, (61, 75), "#1f77b4"),                    # Blue
    ("nSL", nsl_df, (48, 61), "#ff7f0e"),                    # Orange
    ("Moving Tajima's D", mtjd_df, (35, 48), "#2ca02c"),     # Green
    ("Windowed Tajima's D", wtjd_df, (22, 35), "#d62728"),   # Red
    ("B1", b1_df, (9, 22), "#9467bd"),                       # Purple
]

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

# Add population name as title
population = snakemake.params.population
title = f"{population}"
circos.text(title, size=10, r=0)

# Create figure and add legend
fig = circos.plotfig()

# Create legend with unique gene counts
legend_elements = [
    mpatches.Patch(color="#1f77b4", label=f"iHS (n={ihs_df['gene'].nunique()})"),
    mpatches.Patch(color="#ff7f0e", label=f"nSL (n={nsl_df['gene'].nunique()})"),
    mpatches.Patch(color="#2ca02c", label=f"Moving Tajima's D (n={mtjd_df['gene'].nunique()})"),
    mpatches.Patch(color="#d62728", label=f"Windowed Tajima's D (n={wtjd_df['gene'].nunique()})"),
    mpatches.Patch(color="#9467bd", label=f"B1 (n={b1_df['gene'].nunique()})"),
]

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

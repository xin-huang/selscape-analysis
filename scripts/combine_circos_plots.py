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


import svgutils.transform as sg

scores_plot =  snakemake.input.scores_plot
genes_plot = snakemake.input.genes_plot

# Load SVGs
fig1 = sg.fromfile(scores_plot)
fig2 = sg.fromfile(genes_plot)

# Get original sizes
fig1_root = fig1.getroot()
fig2_root = fig2.getroot()

# Create combined figure
fig_combined = sg.SVGFigure("2000px", "1000px")

# Position plots side by side
plot1 = fig1_root
plot2 = fig2_root
plot2.moveto(1000, 0) 

fig_combined.append([plot1, plot2])

fig_combined.save(snakemake.output.combined)

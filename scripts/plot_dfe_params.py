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

import matplotlib.pyplot as plt
import pandas as pd
import matplotlib
matplotlib.use("Agg")

# Read merged DFE data
data = pd.read_csv(snakemake.input.data, sep="\t")

# Population groups
population_list = {
    'AFR': ['ACB', 'ASW', 'ESN', 'GWD', 'LWK', 'MSL', 'YRI'],
    'AMR': ['CLM', 'MXL', 'PEL', 'PUR'],
    'EAS': ['CDX', 'CHB', 'CHS', 'JPT', 'KHV'],
    'EUR': ['CEU', 'FIN', 'GBR', 'IBS', 'TSI'],
    'SAS': ['BEB', 'GIH', 'ITU', 'PJL', 'STU'],
}

# Population order for plotting
sorter = [
    'ACB', 'ASW', 'ESN', 'GWD', 'LWK', 'MSL', 'YRI',  # AFR
    'CLM', 'MXL', 'PEL', 'PUR',                        # AMR
    'CDX', 'CHB', 'CHS', 'JPT', 'KHV',                 # EAS
    'CEU', 'FIN', 'GBR', 'IBS', 'TSI',                 # EUR
    'BEB', 'GIH', 'ITU', 'PJL', 'STU',                 # SAS
]

# Colors by population group
colors = [
    'black', 'black', 'black', 'black', 'black', 'black', 'black',  # AFR
    'green', 'green', 'green', 'green',                              # AMR
    'gold', 'gold', 'gold', 'gold', 'gold',                          # EAS
    'blue', 'blue', 'blue', 'blue', 'blue',                          # EUR
    'brown', 'brown', 'brown', 'brown', 'brown',                     # SAS
]

# Sort data by population order
df = data.sort_values(by="Pop", key=lambda column: column.map(lambda e: sorter.index(e)))

# Create figure with 2x2 subplots, but only use left column
fig, axs = plt.subplots(nrows=2, ncols=2, constrained_layout=True, figsize=(7.5, 4), dpi=350)
gridspec = axs[0, 0].get_subplotspec().get_gridspec()
# Remove right column subplots (will use for legend)
for a in axs[:, 1]:
    a.remove()

# Plot μ (log mean) - upper plot
axs[0, 0].scatter(sorter, df['mu'].values, color=colors, zorder=2)
axs[0, 0].scatter(sorter, df['mu_lb'].values, marker='_', color='grey')
axs[0, 0].scatter(sorter, df['mu_ub'].values, marker='_', color='grey')
# Draw confidence interval lines
for i in range(len(sorter)):
    axs[0, 0].plot([i, i], [df['mu'].values[i], df['mu_ub'].values[i]],
                   linestyle='dashed', color='grey', zorder=1)
    axs[0, 0].plot([i, i], [df['mu_lb'].values[i], df['mu'].values[i]],
                   linestyle='dashed', color='grey', zorder=1)
axs[0, 0].set_ylim([-1, 7])
axs[0, 0].set_ylabel('$\mu$')
axs[0, 0].set_xticks(list(range(len(sorter))), sorter, rotation=90)

# Plot σ (log std dev) - lower plot
axs[1, 0].scatter(sorter, df['sigma'].values, color=colors, zorder=2)
axs[1, 0].scatter(sorter, df['sigma_lb'].values, marker='_', color='grey')
axs[1, 0].scatter(sorter, df['sigma_ub'].values, marker='_', color='grey')
# Draw confidence interval lines
for i in range(len(sorter)):
    axs[1, 0].plot([i, i], [df['sigma'].values[i], df['sigma_ub'].values[i]],
                   linestyle='dashed', color='grey', zorder=1)
    axs[1, 0].plot([i, i], [df['sigma_lb'].values[i], df['sigma'].values[i]],
                   linestyle='dashed', color='grey', zorder=1)
axs[1, 0].set_xticks(list(range(len(sorter))), sorter, rotation=90)
axs[1, 0].set_ylim([0, 40])
axs[1, 0].set_ylabel('$\sigma$')

# Add legend in right column
subfig = fig.add_subfigure(gridspec[:, 1])
handles, labels = subfig.gca().get_legend_handles_labels()
afr = axs[0, 1].scatter([0], [0], label='AFR', color='black')
amr = axs[0, 1].scatter([0], [0], label='AMR', color='green')
eas = axs[0, 1].scatter([0], [0], label='EAS', color='gold')
eur = axs[0, 1].scatter([0], [0], label='EUR', color='blue')
sas = axs[0, 1].scatter([0], [0], label='SAS', color='brown')
handles.extend([afr, amr, eas, eur, sas])
subfig.legend(handles=handles, fontsize=8, handlelength=1.5, loc='upper left')

# Save figure
fig.set_constrained_layout_pads(w_pad=4 / 72, h_pad=4 / 72, hspace=0, wspace=0.1)
plt.savefig(snakemake.output.plot, bbox_inches='tight')
plt.close()

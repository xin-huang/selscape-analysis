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


import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
matplotlib.use("Agg")

df = pd.read_csv(snakemake.input.data, sep="\t")

human_sorter = [
    'ACB', 'ASW', 'ESN', 'GWD', 'LWK', 'MSL', 'YRI',
    'CLM', 'MXL', 'PEL', 'PUR',
    'CDX', 'CHB', 'CHS', 'JPT', 'KHV',
    'CEU', 'FIN', 'GBR', 'IBS', 'TSI',
    'BEB', 'GIH', 'ITU', 'PJL', 'STU',
]
ape_sorter = ['PPA', 'PTE', 'PTS', 'PTT', 'PTV', 'GBB', 'GBG', 'GGG', 'PA', 'PP']
sorter = human_sorter + ape_sorter

human_colors = (
    ['black'] * 7 +
    ['green'] * 4 +
    ['gold']  * 5 +
    ['blue']  * 5 +
    ['brown'] * 5
)
ape_colors = (
    ['#FF7F00'] * 5 +
    ['#9B59B6'] * 3 +
    ['#E91E63'] * 2
)
colors = human_colors + ape_colors

df = df.sort_values(
    by="Pop", key=lambda col: col.map(lambda p: sorter.index(p))
).reset_index(drop=True)

fig, axs = plt.subplots(nrows=2, ncols=2, constrained_layout=True, figsize=(10, 4), dpi=350)
gridspec = axs[0, 0].get_subplotspec().get_gridspec()
for a in axs[:, 1]:
    a.remove()

def plot_param(ax, param, ylim):
    lb, ub = f"{param}_lb", f"{param}_ub"
    x = list(range(len(sorter)))
    ax.scatter(x, df[param].values, color=colors, zorder=2)
    ax.scatter(x, df[lb].values, marker='_', color='grey')
    ax.scatter(x, df[ub].values, marker='_', color='grey')
    for i in range(len(sorter)):
        ax.plot([i, i], [df[param].values[i], df[ub].values[i]],
                linestyle='dashed', color='grey', zorder=1)
        ax.plot([i, i], [df[lb].values[i], df[param].values[i]],
                linestyle='dashed', color='grey', zorder=1)
    ax.set_ylim(ylim)
    ax.set_xticks(x, sorter, rotation=90)

plot_param(axs[0, 0], "mu",    [-6, 13])
plot_param(axs[1, 0], "sigma", [-10, 40])
axs[0, 0].set_ylabel('$\\mu$')
axs[1, 0].set_ylabel('$\\sigma$')

subfig = fig.add_subfigure(gridspec[:, 1])
handles = [
    mpatches.Patch(color=color, label=label)
    for label, color in [
        ('AFR', 'black'), ('AMR', 'green'), ('EAS', 'gold'),
        ('EUR', 'blue'),  ('SAS', 'brown'),
        ('Pan', '#FF7F00'), ('Gorilla', '#9B59B6'), ('Pongo', '#E91E63'),
    ]
]
subfig.legend(handles=handles, fontsize=8, handlelength=1.5, loc='upper left')

fig.set_constrained_layout_pads(w_pad=4/72, h_pad=4/72, hspace=0, wspace=0.1)
plt.savefig(snakemake.output.plot, bbox_inches='tight')
plt.close()


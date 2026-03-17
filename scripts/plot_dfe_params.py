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
import pandas as pd
matplotlib.use("Agg")

data = pd.read_csv(snakemake.input.data, sep="\t")
sorter = snakemake.params.sorter
colors = snakemake.params.colors
legend = snakemake.params.legend
mu_ylim = snakemake.params.mu_ylim
sigma_ylim = snakemake.params.sigma_ylim

df = data.sort_values(by="Pop", key=lambda col: col.map(lambda p: sorter.index(p)))

fig, axs = plt.subplots(nrows=2, ncols=2, constrained_layout=True, figsize=(7.5, 4), dpi=350)
gridspec = axs[0, 0].get_subplotspec().get_gridspec()
for a in axs[:, 1]:
    a.remove()

def plot_param(ax, param, ylim):
    lb, ub = f"{param}_lb", f"{param}_ub"
    ax.scatter(sorter, df[param].values, color=colors, zorder=2)
    ax.scatter(sorter, df[lb].values, marker='_', color='grey')
    ax.scatter(sorter, df[ub].values, marker='_', color='grey')
    for i in range(len(sorter)):
        ax.plot([i, i], [df[param].values[i], df[ub].values[i]], linestyle='dashed', color='grey', zorder=1)
        ax.plot([i, i], [df[lb].values[i], df[param].values[i]], linestyle='dashed', color='grey', zorder=1)
    ax.set_ylim(ylim)
    ax.set_xticks(list(range(len(sorter))), sorter, rotation=90)

plot_param(axs[0, 0], "mu",    mu_ylim)
plot_param(axs[1, 0], "sigma", sigma_ylim)
axs[0, 0].set_ylabel('$\\mu$')
axs[1, 0].set_ylabel('$\\sigma$')

subfig = fig.add_subfigure(gridspec[:, 1])
handles, labels = subfig.gca().get_legend_handles_labels()
for label, color in legend:
    handles.append(axs[0, 1].scatter([0], [0], label=label, color=color))
subfig.legend(handles=handles, fontsize=8, handlelength=1.5, loc='upper left')

fig.set_constrained_layout_pads(w_pad=4/72, h_pad=4/72, hspace=0, wspace=0.1)
plt.savefig(snakemake.output.plot, bbox_inches='tight')
plt.close()

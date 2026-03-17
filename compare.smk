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


rule all:
    input:
        expand(
            "results/plots/dfe/{species}.two_epoch.lognormal.dfe_params.combined.svg",
            species=config["species"],
        )

# DFE comparison
rule plot_dfe_apes:
    input:
        data="resources/ape_dfe_params.tsv",
    output:
        plot="results/plots/dfe/{species}.two_epoch.lognormal.dfe_params.apes.svg",
    params:
        sorter=['PPA', 'PTE', 'PTS', 'PTT', 'PTV', 'GBB', 'GBG', 'GGG', 'PA', 'PP'],
        colors=(
            ['gold'] * 5 +
            ['blue'] * 3 +
            ['red'] * 2
        ),
        legend=[
            ('Pan',     'gold'),
            ('Gorilla', 'blue'),
            ('Pongo',   'red'),
        ],
        mu_ylim=[-6, 13],
        sigma_ylim=[-10, 30],
    script:
        "scripts/plot_dfe_params.py"


rule combine_dfe_plots:
    input:
        human="results/plots/dfe/{species}.two_epoch.lognormal.dfe_params.svg",
        apes=rules.plot_dfe_apes.output.plot,
    output:
        combined="results/plots/dfe/{species}.two_epoch.lognormal.dfe_params.combined.svg",
    script:
        "scripts/combine_dfe_plots.py"

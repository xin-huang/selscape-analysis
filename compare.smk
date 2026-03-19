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
        ),
        expand(
            "results/comparison/{species}/outlier_gene_overlap/", 
            species=config["species"]
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

# outlier comparison
rule compare_outlier_genes:
    input:
        ihs=expand("results/positive_selection/selscan/{{species}}/1pop/{ppl}/ihs_0.05/{ppl}.normalized.ihs.maf_0.05.top_0.0005.outlier.genes", ppl=config["populations"]),
        nsl=expand("results/positive_selection/selscan/{{species}}/1pop/{ppl}/nsl_0.05/{ppl}.normalized.nsl.maf_0.05.top_0.0005.outlier.genes", ppl=config["populations"]),
        mtjd_pos=expand("results/positive_selection/scikit-allel/{{species}}/1pop/{ppl}/moving_tajima_d/100_1/{ppl}.moving_tajima_d.top_0.05.outlier.genes", ppl=config["populations"]),
        wtjd_pos=expand("results/positive_selection/scikit-allel/{{species}}/1pop/{ppl}/windowed_tajima_d/100000_1/{ppl}.windowed_tajima_d.top_0.05.outlier.genes", ppl=config["populations"]),
        betascan=expand("results/balancing_selection/betascan/{{species}}/{ppl}/m_0.15/{ppl}.hg38.m_0.15.b1.top_0.0005.outlier.genes", ppl=config["populations"]),
        mtjd_bal=expand("results/balancing_selection/scikit-allel/{{species}}/moving_tajima_d/{ppl}/100_1/{ppl}.moving_tajima_d.top_0.05.outlier.genes", ppl=config["populations"]),
        wtjd_bal=expand("results/balancing_selection/scikit-allel/{{species}}/windowed_tajima_d/{ppl}/100000_1/{ppl}.windowed_tajima_d.top_0.05.outlier.genes", ppl=config["populations"]),
        study_dir="resources/candidates",
    output:
        outdir=directory("results/comparison/{species}/outlier_gene_overlap/"),
    script:
        "scripts/compare_outlier_genes.py"

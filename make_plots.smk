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

DEMOG = "two_epoch"
DFE = "lognormal"

# plot_group -> list of (species, dataset, ref_genome, population)
DFE_PLOT_GROUPS = {
    "1kg_high_hg38": [
        ("Human", "1kg_high_hg38", "hg38", pop) for pop in [
            "ACB", "ASW", "BEB", "CDX", "CEU", "CHB", "CHS", "CLM", "ESN",
            "FIN", "GBR", "GIH", "GWD", "IBS", "ITU", "JPT", "KHV", "LWK",
            "MSL", "MXL", "PEL", "PJL", "PUR", "STU", "TSI", "YRI",
        ]
    ],
    "1kg_low_hg38": [
        ("Human", "1kg_low_hg38", "hg38", pop) for pop in [
            "ACB", "ASW", "BEB", "CDX", "CEU", "CHB", "CHS", "CLM", "ESN",
            "FIN", "GBR", "GIH", "GWD", "IBS", "ITU", "JPT", "KHV", "LWK",
            "MSL", "MXL", "PEL", "PJL", "PUR", "STU", "TSI", "YRI",
        ]
    ],
    "1kg_low_hg19": [
        ("Human", "1kg_low_hg19", "hg19", pop) for pop in [
            "ACB", "ASW", "BEB", "CDX", "CEU", "CHB", "CHS", "CLM", "ESN",
            "FIN", "GBR", "GIH", "GWD", "IBS", "ITU", "JPT", "KHV", "LWK",
            "MSL", "MXL", "PEL", "PJL", "PUR", "STU", "TSI", "YRI",
        ]
    ],
    "greatape": (
        [("Gorilla", "gorilla", "hg38", pop) for pop in ["GBB", "GBG", "GGG"]]
        + [("Pan", "pan", "hg38", pop) for pop in ["PPA", "PTE", "PTS", "PTV", "PTT"]]
        + [("Pongo", "pongo", "hg38", pop) for pop in ["PP", "PA"]]
    ),
}

POPULATION_GROUPS = {
    "AFR":  {"color": "black",   "populations": ["YRI", "ACB", "ASW", "ESN", "GWD", "LWK", "MSL"]},
    "AMR":  {"color": "green",   "populations": ["CLM", "MXL", "PEL", "PUR"]},
    "EAS":  {"color": "gold",    "populations": ["CDX", "CHB", "CHS", "JPT", "KHV"]},
    "EUR":  {"color": "blue",    "populations": ["CEU", "FIN", "GBR", "IBS", "TSI"]},
    "SAS":  {"color": "brown",   "populations": ["BEB", "GIH", "ITU", "PJL", "STU"]},
    "Gorilla": {"color": "#9B59B6", "populations": ["GBB", "GBG", "GGG"]},
    "Pan":     {"color": "#FF7F00", "populations": ["PPA", "PTE", "PTS", "PTT", "PTV"]},
    "Pongo":   {"color": "#E91E63", "populations": ["PA", "PP"]},
}

rule all_visualization:
    input:
        expand(
            "results/plots/dfe/{group}/{group}.dfe_params.svg",
            group=DFE_PLOT_GROUPS.keys(),
        ),


def get_dfe_bestfit_files(wildcards):
    return [
        f"results/dadi/{sp}/{ds}/dfe/{pop}/InferDFE/{pop}.{rg}.{DEMOG}.{DFE}.InferDFE.bestfits"
        for sp, ds, rg, pop in DFE_PLOT_GROUPS[wildcards.group]
    ]


def get_dfe_ci_files(wildcards):
    return [
        f"results/dadi/{sp}/{ds}/dfe/{pop}/StatDFE/{pop}.{rg}.{DEMOG}.{DFE}.godambe.ci"
        for sp, ds, rg, pop in DFE_PLOT_GROUPS[wildcards.group]
    ]


rule merge_dfe_confidence_intervals:
    input:
        bestfit_files=get_dfe_bestfit_files,
        ci_files=get_dfe_ci_files,
    output:
        merged="results/plots/dfe/{group}/{group}.dfe_params.tsv",
    params:
        populations=lambda wc: [pop for sp, ds, rg, pop in DFE_PLOT_GROUPS[wc.group]],
        datasets=lambda wc: [ds for sp, ds, rg, pop in DFE_PLOT_GROUPS[wc.group]],
    log:
        "logs/make_plots/merge_dfe_confidence_intervals.{group}.log",
    script:
        "scripts/merge_dfe_ci.py"


rule plot_dfe_confidence_intervals:
    input:
        data=rules.merge_dfe_confidence_intervals.output.merged,
    output:
        plot="results/plots/dfe/{group}/{group}.dfe_params.svg",
    params:
        populations=lambda wc: [pop for sp, ds, rg, pop in DFE_PLOT_GROUPS[wc.group]],
        population_groups=POPULATION_GROUPS,
        mu_ylim=None,
        sigma_ylim=None,
    log:
        "logs/make_plots/plot_dfe_confidence_intervals.{group}.log",
    script:
        "scripts/plot_dfe_params.py"

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


import os
from itertools import combinations


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

POPULATIONS = [
    "ACB", "ASW", "BEB", "CDX", "CEU", "CHB", "CHS", "CLM", "ESN", "FIN",
    "GBR", "GIH", "GWD", "IBS", "ITU", "JPT", "KHV", "LWK", "MSL", "MXL",
    "PEL", "PJL", "PUR", "STU", "TSI", "YRI",
]

DATASET_SPECIES = {
    "1kg_high_hg38": "Human",
    "pan": "Pan",
    "gorilla": "Gorilla",
    "pongo": "Pongo",
}

# same order as in the dataset configs, otherwise the pair names do not match
DATASET_POPULATIONS = {
    "1kg_high_hg38": POPULATIONS,
    "pan": ["PPA", "PTE", "PTS", "PTV", "PTT"],
    "gorilla": ["GBB", "GBG", "GGG"],
    "pongo": ["PP", "PA"],
}

XP_METHODS = {
    "xpehh": dict(
        prefix="results/positive_selection/selscan/{species}/{dataset}/2pop/{pair}/xpehh_0.05/{pair}.normalized.xpehh.maf_0.05.top_0.0005",
        score_column="normalized_xpehh",
        pop1_sign="positive",
        title="XP-EHH (MAF=0.05, Top 0.05%)",
    ),
    "xpnsl": dict(
        prefix="results/positive_selection/selscan/{species}/{dataset}/2pop/{pair}/xpnsl_0.05/{pair}.normalized.xpnsl.maf_0.05.top_0.0005",
        score_column="normalized_xpnsl",
        pop1_sign="positive",
        title="XP-nSL (MAF=0.05, Top 0.05%)",
    ),
    "delta_tajima_d": dict(
        prefix="results/positive_selection/scikit-allel/{species}/{dataset}/2pop/{pair}/delta_tajima_d/100_1/{pair}.delta_tajima_d.top_0.0005",
        score_column="delta_tajima_d",
        pop1_sign="negative",
        title="Delta Tajima's D (Window size=100 SNPs, Step size=100 SNPs, Top 0.05%)",
    ),
}

# within-population statistics: one outlier gene list per population
WP_METHODS = {
    "ihs": dict(
        path="results/positive_selection/selscan/{species}/{dataset}/1pop/{ppl}/ihs_0.05/{ppl}.normalized.ihs.maf_0.05.top_0.0005.outlier.genes",
        title="iHS (MAF=0.05, Top 0.05%)",
    ),
    "nsl": dict(
        path="results/positive_selection/selscan/{species}/{dataset}/1pop/{ppl}/nsl_0.05/{ppl}.normalized.nsl.maf_0.05.top_0.0005.outlier.genes",
        title="nSL (MAF=0.05, Top 0.05%)",
    ),
    "moving_tajima_d": dict(
        path="results/positive_selection/scikit-allel/{species}/{dataset}/1pop/{ppl}/moving_tajima_d/100_1/{ppl}.moving_tajima_d.top_0.0005.outlier.genes",
        title="Moving Tajima's D, positive selection (Window size=100 SNPs, Step size=100 SNPs, Top 0.05%)",
    ),
    "windowed_tajima_d": dict(
        path="results/positive_selection/scikit-allel/{species}/{dataset}/1pop/{ppl}/windowed_tajima_d/100000_1/{ppl}.windowed_tajima_d.top_0.0005.outlier.genes",
        title="Windowed Tajima's D, positive selection (Window size=100000 bp, Step size=100000 bp, Top 0.05%)",
    ),
    "b1": dict(
        path="results/balancing_selection/betascan/{species}/{dataset}/{ppl}/m_0.15/{ppl}.hg38.m_0.15.b1.top_0.0005.outlier.genes",
        title="B1 (Core Freq=0.15, Top 0.05%)",
    ),
    "moving_tajima_d_balancing": dict(
        path="results/balancing_selection/scikit-allel/{species}/{dataset}/moving_tajima_d/{ppl}/100_1/{ppl}.moving_tajima_d.top_0.0005.outlier.genes",
        title="Moving Tajima's D, balancing selection (Window size=100 SNPs, Step size=100 SNPs, Top 0.05%)",
    ),
    "windowed_tajima_d_balancing": dict(
        path="results/balancing_selection/scikit-allel/{species}/{dataset}/windowed_tajima_d/{ppl}/100000_1/{ppl}.windowed_tajima_d.top_0.0005.outlier.genes",
        title="Windowed Tajima's D, balancing selection (Window size=100000 bp, Step size=100000 bp, Top 0.05%)",
    ),
}


wildcard_constraints:
    dataset="|".join(DATASET_SPECIES),
    method="|".join(XP_METHODS),
    wp_method="|".join(WP_METHODS),


DFE_PLOTS = expand("results/plots/dfe/{group}/{group}.dfe_params.svg", group=DFE_PLOT_GROUPS.keys())
XP_PLOTS = expand(
    "results/comparisons/xp/{plot}/{dataset}.{method}.{plot}.svg",
    plot=["upset", "matrix"], dataset=DATASET_SPECIES, method=XP_METHODS,
)
WP_PLOTS = expand(
    "results/comparisons/wp/{plot}/{dataset}.{wp_method}.{plot}.svg",
    plot=["upset", "matrix"], dataset=DATASET_SPECIES, wp_method=WP_METHODS,
)
WP_GO_PLOTS = expand(
    "results/comparisons/wp_go/matrix/{dataset}.{wp_method}.go_matrix.svg",
    dataset=DATASET_SPECIES, wp_method=WP_METHODS,
)


rule all_visualization:
    input:
        DFE_PLOTS,
        XP_PLOTS,
        WP_PLOTS,
        WP_GO_PLOTS,


rule all_dfe:
    input:
        DFE_PLOTS,


rule all_xp:
    input:
        XP_PLOTS,


rule all_wp:
    input:
        WP_PLOTS,


rule all_wp_go:
    input:
        WP_GO_PLOTS,


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


def get_xp_files(wildcards, suffix):
    prefix = XP_METHODS[wildcards.method]["prefix"]
    return [
        prefix.format(species=DATASET_SPECIES[wildcards.dataset], dataset=wildcards.dataset, pair=f"{a}_{b}") + suffix
        for a, b in combinations(DATASET_POPULATIONS[wildcards.dataset], 2)
    ]

def get_xp_focal_files(wildcards):
    """[file name in the split_xp_focal_genes directory, focal, reference] for both directions of every pair."""
    return [
        [f"{a}_{b}.focal_{focal}.outlier.genes", focal, reference]
        for a, b in combinations(DATASET_POPULATIONS[wildcards.dataset], 2)
        for focal, reference in ((a, b), (b, a))
    ]


def wp_path(dataset, wp_method, ppl):
    return WP_METHODS[wp_method]["path"].format(
        species=DATASET_SPECIES[dataset], dataset=dataset, ppl=ppl
    )


def get_wp_populations(wildcards):
    """Populations with an existing outlier gene list: not every method ran for every dataset (e.g. iHS needs ancestral alleles)."""
    return [
        ppl
        for ppl in DATASET_POPULATIONS[wildcards.dataset]
        if os.path.exists(wp_path(wildcards.dataset, wildcards.wp_method, ppl))
    ]


def get_wp_files(wildcards):
    return [
        wp_path(wildcards.dataset, wildcards.wp_method, ppl)
        for ppl in get_wp_populations(wildcards)
    ]


def wp_enrichment_path(dataset, wp_method, ppl):
    path = WP_METHODS[wp_method]["path"].replace(".outlier.genes", ".gowinda.enrichment.tsv")
    return path.format(species=DATASET_SPECIES[dataset], dataset=dataset, ppl=ppl)


def get_wp_go_populations(wildcards):
    """Populations with a Gowinda result."""
    return [
        ppl
        for ppl in DATASET_POPULATIONS[wildcards.dataset]
        if os.path.exists(wp_enrichment_path(wildcards.dataset, wildcards.wp_method, ppl))
    ]


def get_wp_go_files(wildcards):
    """One significant GO term list per population, for populations with a Gowinda result."""
    return [
        f"results/comparisons/wp_go/terms/{wildcards.dataset}.{wildcards.wp_method}/{ppl}.significant_go.txt"
        for ppl in get_wp_go_populations(wildcards)
    ]



rule split_xp_focal_genes:
    input:
        scores=lambda wc: get_xp_files(wc, ".outliers.scores"),
        annotated=lambda wc: get_xp_files(wc, ".annotated.outliers"),
    output:
        outdir=directory("results/comparisons/xp/genes/{dataset}.{method}"),
    params:
        score_column=lambda wc: XP_METHODS[wc.method]["score_column"],
        pop1_sign=lambda wc: XP_METHODS[wc.method]["pop1_sign"],
    log:
        "logs/comparisons/split_xp_focal_genes.{dataset}.{method}.log",
    script:
        "scripts/split_xp_focal_genes.py"


rule plot_xp_upset:
    input:
        genes=rules.split_xp_focal_genes.output.outdir,
    output:
        plot="results/comparisons/xp/upset/{dataset}.{method}.upset.svg",
        table="results/comparisons/xp/upset/{dataset}.{method}.upset.tsv",
    params:
        focal_files=get_xp_focal_files,
        population_groups=POPULATION_GROUPS,
        max_intersections=15,
        title=lambda wc: f"{wc.dataset} {XP_METHODS[wc.method]['title']}",
    log:
        "logs/comparisons/plot_xp_upset.{dataset}.{method}.log",
    conda:
        "compare_tmp.yaml"
    script:
        "scripts/plot_upset.py"


rule plot_xp_matrix:
    input:
        genes=rules.split_xp_focal_genes.output.outdir,
    output:
        plot="results/comparisons/xp/matrix/{dataset}.{method}.matrix.svg",
        table="results/comparisons/xp/matrix/{dataset}.{method}.matrix.tsv",
    params:
        focal_files=get_xp_focal_files,
        population_groups=POPULATION_GROUPS,
        title=lambda wc: f"{wc.dataset} {XP_METHODS[wc.method]['title']}",
    log:
        "logs/comparisons/plot_xp_matrix.{dataset}.{method}.log",
    conda:
        "compare_tmp.yaml"
    script:
        "scripts/plot_xp_matrix.py"


rule plot_wp_upset:
    input:
        genes=get_wp_files,
    output:
        plot="results/comparisons/wp/upset/{dataset}.{wp_method}.upset.svg",
        table="results/comparisons/wp/upset/{dataset}.{wp_method}.upset.tsv",
    params:
        populations=get_wp_populations,
        population_groups=POPULATION_GROUPS,
        max_intersections=15,
        title=lambda wc: f"{wc.dataset} {WP_METHODS[wc.wp_method]['title']}",
    log:
        "logs/comparisons/plot_wp_upset.{dataset}.{wp_method}.log",
    conda:
        "compare_tmp.yaml"
    script:
        "scripts/plot_upset.py"


rule plot_wp_matrix:
    input:
        genes=get_wp_files,
    output:
        plot="results/comparisons/wp/matrix/{dataset}.{wp_method}.matrix.svg",
        table="results/comparisons/wp/matrix/{dataset}.{wp_method}.matrix.tsv",
    params:
        populations=get_wp_populations,
        population_groups=POPULATION_GROUPS,
        title=lambda wc: f"{wc.dataset} {WP_METHODS[wc.wp_method]['title']}",
    log:
        "logs/comparisons/plot_wp_matrix.{dataset}.{wp_method}.log",
    conda:
        "compare_tmp.yaml"
    script:
        "scripts/plot_wp_matrix.py"


rule extract_wp_go_terms:
    input:
        enrichment=lambda wc: wp_enrichment_path(wc.dataset, wc.wp_method, wc.ppl),
    output:
        go_terms="results/comparisons/wp_go/terms/{dataset}.{wp_method}/{ppl}.significant_go.txt",
    params:
        fdr=0.05,
    log:
        "logs/comparisons/extract_wp_go_terms.{dataset}.{wp_method}.{ppl}.log",
    shell:
        """
        awk -F'\\t' -v fdr={params.fdr} 'BEGIN {{ print "GO_ID" }} $1 != "GO_ID" && $5 != "" && $5 + 0 < fdr {{ print $1 }}' \\
            {input.enrichment} > {output.go_terms} 2> {log}
        """


rule plot_wp_go_matrix:
    input:
        genes=get_wp_go_files,
    output:
        plot="results/comparisons/wp_go/matrix/{dataset}.{wp_method}.go_matrix.svg",
        table="results/comparisons/wp_go/matrix/{dataset}.{wp_method}.go_matrix.tsv",
    params:
        populations=get_wp_go_populations,
        population_groups=POPULATION_GROUPS,
        unit="GO terms",
        title=lambda wc: f"{wc.dataset} GO terms, {WP_METHODS[wc.wp_method]['title']}",
    log:
        "logs/comparisons/plot_wp_go_matrix.{dataset}.{wp_method}.log",
    conda:
        "compare_tmp.yaml"
    script:
        "scripts/plot_wp_matrix.py"

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

configfile: "config/compare.yaml"

COMPARISONS = config.get("comparisons", [])
TOP_N = config.get("top_n", 20)

P_ADJUSTED_THRESHOLDS = config.get("p_adjusted_threshold", [None])
if not isinstance(P_ADJUSTED_THRESHOLDS, list):
    P_ADJUSTED_THRESHOLDS = [P_ADJUSTED_THRESHOLDS]


def _thr_label(t):
    return "none" if t is None else str(t)


def _thr_value(label):
    return None if label == "none" else float(label)


POPULATIONS = [
    "ACB", "ASW", "BEB", "CDX", "CEU", "CHB", "CHS", "CLM", "ESN", "FIN",
    "GBR", "GIH", "GWD", "IBS", "ITU", "JPT", "KHV", "LWK", "MSL", "MXL",
    "PEL", "PJL", "PUR", "STU", "TSI", "YRI",
]


rule all_comparisons:
    input:
        [f"results/comparisons/{a}_vs_{b}/positive_selection/thr_{_thr_label(t)}"
         for a, b in COMPARISONS for t in P_ADJUSTED_THRESHOLDS],
        [f"results/comparisons/{a}_vs_{b}/balancing_selection/thr_{_thr_label(t)}"
         for a, b in COMPARISONS for t in P_ADJUSTED_THRESHOLDS],
        "results/comparisons/1kg_high_hg38/outlier_gene_overlap/",
        "results/comparisons/jaccard/",


rule compare_positive_selection:
    output:
        directory("results/comparisons/{dataset_a}_vs_{dataset_b}/positive_selection/thr_{threshold}"),
    params:
        top_n=TOP_N,
        p_adjusted_threshold=lambda wc: _thr_value(wc.threshold),
    log:
        "logs/comparisons/compare_positive_selection.{dataset_a}_vs_{dataset_b}.thr_{threshold}.log",
    script:
        "scripts/compare_positive_selection.py"


rule compare_balancing_selection:
    output:
        directory("results/comparisons/{dataset_a}_vs_{dataset_b}/balancing_selection/thr_{threshold}"),
    params:
        top_n=TOP_N,
        p_adjusted_threshold=lambda wc: _thr_value(wc.threshold),
    log:
        "logs/comparisons/compare_balancing_selection.{dataset_a}_vs_{dataset_b}.thr_{threshold}.log",
    script:
        "scripts/compare_balancing_selection.py"



rule compare_outlier_genes:
    input:
        ihs=expand(
            "results/positive_selection/selscan/Human/1kg_high_hg38/1pop/{ppl}/ihs_0.05/{ppl}.normalized.ihs.maf_0.05.top_0.0005.outlier.genes",
            ppl=POPULATIONS,
        ),
        nsl=expand(
            "results/positive_selection/selscan/Human/1kg_high_hg38/1pop/{ppl}/nsl_0.05/{ppl}.normalized.nsl.maf_0.05.top_0.0005.outlier.genes",
            ppl=POPULATIONS,
        ),
        mtjd_pos=expand(
            "results/positive_selection/scikit-allel/Human/1kg_high_hg38/1pop/{ppl}/moving_tajima_d/100_1/{ppl}.moving_tajima_d.top_0.0005.outlier.genes",
            ppl=POPULATIONS,
        ),
        wtjd_pos=expand(
            "results/positive_selection/scikit-allel/Human/1kg_high_hg38/1pop/{ppl}/windowed_tajima_d/100000_1/{ppl}.windowed_tajima_d.top_0.0005.outlier.genes",
            ppl=POPULATIONS,
        ),
        betascan=expand(
            "results/balancing_selection/betascan/Human/1kg_high_hg38/{ppl}/m_0.15/{ppl}.hg38.m_0.15.b1.top_0.0005.outlier.genes",
            ppl=POPULATIONS,
        ),
        mtjd_bal=expand(
            "results/balancing_selection/scikit-allel/Human/1kg_high_hg38/moving_tajima_d/{ppl}/100_1/{ppl}.moving_tajima_d.top_0.0005.outlier.genes",
            ppl=POPULATIONS,
        ),
        wtjd_bal=expand(
            "results/balancing_selection/scikit-allel/Human/1kg_high_hg38/windowed_tajima_d/{ppl}/100000_1/{ppl}.windowed_tajima_d.top_0.0005.outlier.genes",
            ppl=POPULATIONS,
        ),
        study_dir="resources/candidates",
    output:
        outdir=directory("results/comparisons/1kg_high_hg38/outlier_gene_overlap/"),
    script:
        "scripts/compare_outlier_genes.py"


DATASETS_TO_COMPARE = config.get("datasets_to_compare", [])

rule plot_gene_jaccard_heatmaps:
    output:
        outdir=directory("results/comparisons/jaccard/"),
    params:
        datasets=DATASETS_TO_COMPARE,
    log:
        "logs/comparisons/plot_gene_jaccard_heatmaps.log",
    script:
        "scripts/plot_jaccard_heatmaps.py"

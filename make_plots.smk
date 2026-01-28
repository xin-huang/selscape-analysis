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


# Configuration
POPULATIONS = config["populations"]

# Statistics paths
STATISTICS = {
    "ihs": {
        "genes": "results/positive_selection/selscan/{species}/1pop/{ppl}/ihs_0.05/{ppl}.normalized.ihs.maf_0.05.top_0.0005.outlier.genes",
        "candidates": "results/positive_selection/selscan/{species}/1pop/{ppl}/ihs_0.05/{ppl}.normalized.ihs.maf_0.05.top_0.0005.annotated.outliers",
        "scores": "results/positive_selection/selscan/{species}/1pop/{ppl}/ihs_0.05/{ppl}.normalized.ihs.scores"
    },
    "nsl": {
        "genes": "results/positive_selection/selscan/{species}/1pop/{ppl}/nsl_0.05/{ppl}.normalized.nsl.maf_0.05.top_0.0005.outlier.genes",
        "candidates": "results/positive_selection/selscan/{species}/1pop/{ppl}/nsl_0.05/{ppl}.normalized.nsl.maf_0.05.top_0.0005.annotated.outliers",
        "scores": "results/positive_selection/selscan/{species}/1pop/{ppl}/nsl_0.05/{ppl}.normalized.nsl.scores"
    },
    "moving_tajima_d": {
        "genes": "results/positive_selection/scikit-allel/{species}/1pop/{ppl}/moving_tajima_d/100_1/{ppl}.moving_tajima_d.top_0.05.outlier.genes",
        "candidates": "results/positive_selection/scikit-allel/{species}/1pop/{ppl}/moving_tajima_d/100_1/{ppl}.moving_tajima_d.top_0.05.annotated.outliers",
        "scores": "results/positive_selection/scikit-allel/{species}/1pop/{ppl}/moving_tajima_d/100_1/{ppl}.moving_tajima_d.scores"
    },
    "windowed_tajima_d": {
        "genes": "results/positive_selection/scikit-allel/{species}/1pop/{ppl}/windowed_tajima_d/100000_1/{ppl}.windowed_tajima_d.top_0.05.outlier.genes",
        "candidates": "results/positive_selection/scikit-allel/{species}/1pop/{ppl}/windowed_tajima_d/100000_1/{ppl}.windowed_tajima_d.top_0.05.annotated.outliers",
        "scores": "results/positive_selection/scikit-allel/{species}/1pop/{ppl}/windowed_tajima_d/100000_1/{ppl}.windowed_tajima_d.scores"
    },
    "betascan_b1": {
        "genes": "results/balancing_selection/betascan/{species}/{ppl}/m_0.15/{ppl}.hg38.m_0.15.b1.top_0.0005.outlier.genes",
        "candidates": "results/balancing_selection/betascan/{species}/{ppl}/m_0.15/{ppl}.hg38.m_0.15.b1.top_0.0005.annotated.outliers",
        "scores": "results/balancing_selection/betascan/{species}/{ppl}/m_0.15/{ppl}.hg38.m_0.15.b1.scores"
    },
    "moving_tajima_d_bal": {
        "genes": "results/balancing_selection/scikit-allel/{species}/moving_tajima_d/{ppl}/100_1/{ppl}.moving_tajima_d.top_0.05.outlier.genes",
        "candidates": "results/balancing_selection/scikit-allel/{species}/moving_tajima_d/{ppl}/100_1/{ppl}.moving_tajima_d.top_0.05.annotated.outliers",
        "scores": "results/balancing_selection/scikit-allel/{species}/moving_tajima_d/{ppl}/100_1/{ppl}.moving_tajima_d.merged.scores"
    },
    "windowed_tajima_d_bal": {
        "genes": "results/balancing_selection/scikit-allel/{species}/windowed_tajima_d/{ppl}/100000_1/{ppl}.windowed_tajima_d.top_0.05.outlier.genes",
        "candidates": "results/balancing_selection/scikit-allel/{species}/windowed_tajima_d/{ppl}/100000_1/{ppl}.windowed_tajima_d.top_0.05.annotated.outliers",
        "scores": "results/balancing_selection/scikit-allel/{species}/windowed_tajima_d/{ppl}/100000_1/{ppl}.windowed_tajima_d.merged.scores"
    }
}

rule all_visualization:
    input:
        expand("results/plots/circos/{ppl}/{ppl}_positive_selection_combined_circos.svg", ppl=POPULATIONS),
        expand("results/plots/circos/{ppl}/{ppl}_balancing_selection_combined_circos.svg", ppl=POPULATIONS),
        "results/plots/dfe/Human.two_epoch.lognormal.dfe_params.svg"


rule extract_gene_regions:
    input:
        refgene="resources/tools/annovar/hg38_db/hg38_refGene.txt",
        genes=lambda wildcards: STATISTICS[wildcards.stat]["genes"].format(
            species=config["species"],
            ppl=wildcards.ppl
        ),
        candidates=lambda wildcards: STATISTICS[wildcards.stat]["candidates"].format(
            species=config["species"],
            ppl=wildcards.ppl
        )
    output:
        bed="results/plots/circos/{ppl}/genes/{stat}.genes.bed"
    log:
        "logs/plots/extract_gene_regions.{ppl}.{stat}.log"
    conda:
        "env.yaml"
    script:
        "scripts/extract_gene_regions.py"


rule make_positive_selection_circos:
    input:
        ihs_bed="results/plots/circos/{ppl}/genes/ihs.genes.bed",
        nsl_bed="results/plots/circos/{ppl}/genes/nsl.genes.bed",
        mtjd_bed="results/plots/circos/{ppl}/genes/moving_tajima_d.genes.bed",
        wtjd_bed="results/plots/circos/{ppl}/genes/windowed_tajima_d.genes.bed"
    output:
        plot="results/plots/circos/{ppl}/{ppl}_positive_selection_circos.svg"
    params:
        population="{ppl}",
        plot_type="positive_selection"
    log:
        "logs/plots/make_positive_selection_circos.{ppl}.log"
    conda:
        "env.yaml"
    script:
        "scripts/plot_circos_diagram.py"


rule make_balancing_selection_circos:
    input:
        b1_bed="results/plots/circos/{ppl}/genes/betascan_b1.genes.bed",
        mtjd_bal_bed="results/plots/circos/{ppl}/genes/moving_tajima_d_bal.genes.bed",
        wtjd_bal_bed="results/plots/circos/{ppl}/genes/windowed_tajima_d_bal.genes.bed"
    output:
        plot="results/plots/circos/{ppl}/{ppl}_balancing_selection_circos.svg"
    params:
        population="{ppl}",
        plot_type="balancing_selection"
    log:
        "logs/plots/make_balancing_selection_circos.{ppl}.log"
    conda:
        "env.yaml"
    script:
        "scripts/plot_circos_diagram.py"



def get_score_file(stat_name):
    def _get_file(wildcards):
        return STATISTICS[stat_name]["scores"].format(
            species=config["species"],
            ppl=wildcards.ppl
        )
    return _get_file


rule make_positive_selection_circos_scores:
    input:
        ihs_scores=get_score_file("ihs"),
        nsl_scores=get_score_file("nsl"),
        mtjd_scores=get_score_file("moving_tajima_d"),
        wtjd_scores=get_score_file("windowed_tajima_d")
    output:
        plot="results/plots/circos/{ppl}/{ppl}_positive_selection_circos_scores.svg"
    params:
        population="{ppl}",
        plot_type="positive_selection"
    resources:
        mem_gb=32,
    log:
        "logs/plots/make_positive_selection_circos_scores.{ppl}.log"
    conda:
        "env.yaml"
    script:
        "scripts/plot_circos_scores.py"


rule make_balancing_selection_circos_scores:
    input:
        b1_scores=get_score_file("betascan_b1"),
        mtjd_bal_scores=get_score_file("moving_tajima_d_bal"),
        wtjd_bal_scores=get_score_file("windowed_tajima_d_bal")
    output:
        plot="results/plots/circos/{ppl}/{ppl}_balancing_selection_circos_scores.svg"
    params:
        population="{ppl}",
        plot_type="balancing_selection"
    log:
        "logs/plots/make_balancing_selection_circos_scores.{ppl}.log"
    conda:
        "env.yaml"
    script:
        "scripts/plot_circos_scores.py"


rule make_combined_circos:
    input:
        scores_plot="results/plots/circos/{ppl}/{ppl}_{selection_type}_circos_scores.svg",
        genes_plot="results/plots/circos/{ppl}/{ppl}_{selection_type}_circos.svg"
    output:
        combined="results/plots/circos/{ppl}/{ppl}_{selection_type}_combined_circos.svg"
    params:
        population="{ppl}",
        selection_type="{selection_type}"
    resources:
        mem_gb=16,
    log:
        "logs/plots/make_combined_circos.{ppl}.{selection_type}.log"
    conda:
        "env.yaml"
    script:
        "scripts/combine_circos_plots.py"


rule merge_dfe_confidence_intervals:
    input:
        bestfit_files=expand("results/dadi/{{species}}/dfe/{ppl}/InferDFE/{ppl}.hg38.two_epoch.lognormal.InferDFE.bestfits", ppl=POPULATIONS),
        ci_files=expand("results/dadi/{{species}}/dfe/{ppl}/StatDFE/{ppl}.hg38.two_epoch.lognormal.godambe.ci", ppl=POPULATIONS),
    output:
        merged="results/plots/dfe/{species}.two_epoch.lognormal.dfe_params.tsv"
    params:
        populations=POPULATIONS
    log:
        "logs/plots/dfe/merge_dfe_confidence_intervals.{species}.log"
    conda:
        "env.yaml"
    script:
        "scripts/merge_dfe_ci.py"


rule plot_dfe_confidence_intervals:
    input:
        data=rules.merge_dfe_confidence_intervals.output.merged
    output:
        plot="results/plots/dfe/{species}.two_epoch.lognormal.dfe_params.svg"
    log:
        "logs/plots/dfe/plot_dfe_confidence_intervals.{species}.log"
    conda:
        "env.yaml"
    script:
        "scripts/plot_dfe_params.py"

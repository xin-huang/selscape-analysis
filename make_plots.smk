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
        "candidates": "results/positive_selection/selscan/{species}/1pop/{ppl}/ihs_0.05/{ppl}.normalized.ihs.maf_0.05.top_0.0005.annotated.outliers"
    },
    "nsl": {
        "genes": "results/positive_selection/selscan/{species}/1pop/{ppl}/nsl_0.05/{ppl}.normalized.nsl.maf_0.05.top_0.0005.outlier.genes",
        "candidates": "results/positive_selection/selscan/{species}/1pop/{ppl}/nsl_0.05/{ppl}.normalized.nsl.maf_0.05.top_0.0005.annotated.outliers"
    },
    "moving_tajima_d": {
        "genes": "results/positive_selection/scikit-allel/{species}/1pop/{ppl}/moving_tajima_d/100_1/{ppl}.moving_tajima_d.top_0.05.outlier.genes",
        "candidates": "results/positive_selection/scikit-allel/{species}/1pop/{ppl}/moving_tajima_d/100_1/{ppl}.moving_tajima_d.top_0.05.annotated.outliers"
    },
    "windowed_tajima_d": {
        "genes": "results/positive_selection/scikit-allel/{species}/1pop/{ppl}/windowed_tajima_d/100000_1/{ppl}.windowed_tajima_d.top_0.05.outlier.genes",
        "candidates": "results/positive_selection/scikit-allel/{species}/1pop/{ppl}/windowed_tajima_d/100000_1/{ppl}.windowed_tajima_d.top_0.05.annotated.outliers"
    },
    "betascan_b1": {
        "genes": "results/balancing_selection/betascan/{species}/{ppl}/m_0.15/{ppl}.hg38.m_0.15.b1.top_0.0005.outlier.genes",
        "candidates": "results/balancing_selection/betascan/{species}/{ppl}/m_0.15/{ppl}.hg38.m_0.15.b1.top_0.0005.annotated.outliers"
    }
}

rule all_circos_plots:
    input:
        expand("results/plots/circos/{ppl}/{ppl}_circos.svg", ppl=POPULATIONS)

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

rule make_circos_plot:
    input:
        ihs_bed="results/plots/circos/{ppl}/genes/ihs.genes.bed",
        nsl_bed="results/plots/circos/{ppl}/genes/nsl.genes.bed",
        mtjd_bed="results/plots/circos/{ppl}/genes/moving_tajima_d.genes.bed",
        wtjd_bed="results/plots/circos/{ppl}/genes/windowed_tajima_d.genes.bed",
        b1_bed="results/plots/circos/{ppl}/genes/betascan_b1.genes.bed"
    output:
        plot="results/plots/circos/{ppl}/{ppl}_circos.svg"
    params:
        population="{ppl}"
    log:
        "logs/plots/make_circos.{ppl}.log"
    conda:
        "env.yaml"
    script:
        "scripts/plot_circos_diagram.py"

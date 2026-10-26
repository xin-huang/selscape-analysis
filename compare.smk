# Copyright 2026 Xin Huang and Simon Chen
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


GOWINDA_COMPARISONS = config.get("gowinda_comparisons", [])


GOWINDA_COMPARISON_METHODS = {
    ("positive_selection", "selscan", "ihs"):
        "results/positive_selection/selscan/{species}/{dataset}/1pop/{ppl}/ihs_0.05/{ppl}.normalized.ihs.maf_0.05.top_0.0005",
    ("positive_selection", "selscan", "nsl"):
        "results/positive_selection/selscan/{species}/{dataset}/1pop/{ppl}/nsl_0.05/{ppl}.normalized.nsl.maf_0.05.top_0.0005",
    ("positive_selection", "scikit-allel", "moving_tajima_d"):
        "results/positive_selection/scikit-allel/{species}/{dataset}/1pop/{ppl}/moving_tajima_d/100_1/{ppl}.moving_tajima_d.top_0.0005",
    ("positive_selection", "scikit-allel", "windowed_tajima_d"):
        "results/positive_selection/scikit-allel/{species}/{dataset}/1pop/{ppl}/windowed_tajima_d/100000_1/{ppl}.windowed_tajima_d.top_0.0005",
    ("balancing_selection", "betascan", "betascan_b1"):
        "results/balancing_selection/betascan/{species}/{dataset}/{ppl}/m_0.15/{ppl}.hg38.m_0.15.b1.top_0.0005",
    ("balancing_selection", "scikit-allel", "moving_tajima_d"):
        "results/balancing_selection/scikit-allel/{species}/{dataset}/moving_tajima_d/{ppl}/100_1/{ppl}.moving_tajima_d.top_0.0005",
    ("balancing_selection", "scikit-allel", "windowed_tajima_d"):
        "results/balancing_selection/scikit-allel/{species}/{dataset}/windowed_tajima_d/{ppl}/100000_1/{ppl}.windowed_tajima_d.top_0.0005",
}

DATASET_SPECIES = {
    "1kg_high_hg38": "Human",
    "pan": "Pan",
    "gorilla": "Gorilla",
    "pongo": "Pongo",
}

DATASET_POPULATIONS = {
    "1kg_high_hg38": POPULATIONS,
    "pan": ["PPA", "PTE", "PTS", "PTV", "PTT"],
    "gorilla": ["GBB", "GBG", "GGG"],
    "pongo": ["PP", "PA"],
}

rule all_comparisons:
    input:
        [f"results/comparisons/{a}_vs_{b}/positive_selection/thr_{_thr_label(t)}"
         for a, b in COMPARISONS for t in P_ADJUSTED_THRESHOLDS],
        [f"results/comparisons/{a}_vs_{b}/balancing_selection/thr_{_thr_label(t)}"
         for a, b in COMPARISONS for t in P_ADJUSTED_THRESHOLDS],
        "results/comparisons/1kg_high_hg38/outlier_gene_overlap/",
        "results/comparisons/jaccard/plots/",
        "results/comparisons/jaccard/gene_lists/",	
        [f"results/comparisons/jaccard/gowinda/{c}.{t}.{m}.{a}_vs_{b}.{s}.gowinda.enrichment.png"
         for c, t, m in GOWINDA_COMPARISON_METHODS
         for a, b in GOWINDA_COMPARISONS
         for s in ["overlap", f"unique_{a}", f"unique_{b}"]],


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
        outdir=directory("results/comparisons/jaccard/plots/"),
    params:
        datasets=DATASETS_TO_COMPARE,
        populations=POPULATIONS,
        population_datasets=config.get("population_level_datasets", []),
    log:
        "logs/comparisons/plot_gene_jaccard_heatmaps.log",
    script:
        "scripts/plot_jaccard_heatmaps.py"


rule write_jaccard_gene_lists:
    output:
        outdir=directory("results/comparisons/jaccard/gene_lists/"),
    params:
        comparisons=GOWINDA_COMPARISONS,
    log:
        "logs/comparisons/write_jaccard_gene_lists.log",
    script:
        "scripts/jaccard_gene_lists.py"


wildcard_constraints:
    category="positive_selection|balancing_selection",
    tool="selscan|scikit-allel|betascan",
    method="[A-Za-z0-9_]+",
    dataset_a="[A-Za-z0-9_]+",
    dataset_b="[A-Za-z0-9_]+",
    subset="overlap|unique_[A-Za-z0-9_]+",
 
 
def get_gowinda_comparison_files(wildcards, suffix):
    if wildcards.subset == "overlap":
        datasets = [wildcards.dataset_a, wildcards.dataset_b]
    else:
        datasets = [wildcards.subset.removeprefix("unique_")]
 
    prefix = GOWINDA_COMPARISON_METHODS[(wildcards.category, wildcards.tool, wildcards.method)]
    return [
        prefix.format(species=DATASET_SPECIES[dataset], dataset=dataset, ppl=ppl) + suffix
        for dataset in datasets
        for ppl in DATASET_POPULATIONS[dataset]
    ]
 
 
rule extract_gowinda_comparison_snps:
    input:
        gene_lists=rules.write_jaccard_gene_lists.output.outdir,
        outliers=lambda wc: get_gowinda_comparison_files(wc, ".annotated.outliers"),
        total=lambda wc: get_gowinda_comparison_files(wc, ".total.snps.tsv"),
    output:
        genes=temp("results/comparisons/jaccard/gowinda/{category}.{tool}.{method}.{dataset_a}_vs_{dataset_b}.{subset}.snp.keys"),
        candidate_snps="results/comparisons/jaccard/gowinda/{category}.{tool}.{method}.{dataset_a}_vs_{dataset_b}.{subset}.candidate.snps.tsv",
        total_snps="results/comparisons/jaccard/gowinda/{category}.{tool}.{method}.{dataset_a}_vs_{dataset_b}.{subset}.total.snps.tsv",
    params:
        gene_list="{category}.{tool}.{method}.{dataset_a}_vs_{dataset_b}.{subset}.genes",
    resources:
        mem_mb=16000,
    log:
        "logs/comparisons/extract_gowinda_comparison_snps.{category}.{tool}.{method}.{dataset_a}_vs_{dataset_b}.{subset}.log",
    shell:
        r"""
        ( awk 'NR==FNR {{ if (FNR>1) genes[$1]; next }}
               FNR>1 && ($7 in genes) {{ print $1":"$2 }}' \
            {input.gene_lists}/{params.gene_list} {input.outliers} | sort -u > {output.genes} ) 2> {log}
 
        ( cat {input.total} | sort -u > {output.total_snps} ) 2>> {log}
 
        ( awk 'NR==FNR {{ keys[$1]; next }}
               {{ chrom=$1; sub(/^chr/, "", chrom); if ((chrom":"$2) in keys) print $1"\t"$2 }}' \
            {output.genes} {output.total_snps} > {output.candidate_snps} ) 2>> {log}
        """
 
 
rule enrichment_gowinda_comparisons:
    input:
        gowinda="resources/tools/gowinda/Gowinda-1.12.jar",
        go2gene="results/annotated_data/Human/1kg_high_hg38.gowinda.go2gene",
        gtf="results/annotated_data/Human/1kg_high_hg38.gowinda.gtf",
        candidate_snps=rules.extract_gowinda_comparison_snps.output.candidate_snps,
        total_snps=rules.extract_gowinda_comparison_snps.output.total_snps,
    output:
        enrichment="results/comparisons/jaccard/gowinda/{category}.{tool}.{method}.{dataset_a}_vs_{dataset_b}.{subset}.gowinda.enrichment.tsv",
    resources:
        mem_mb=32000,
        cpus=8,
    log:
        "logs/comparisons/enrichment_gowinda_comparisons.{category}.{tool}.{method}.{dataset_a}_vs_{dataset_b}.{subset}.log",
    shell:
        r"""
        touch {output.enrichment} 2> {log}
 
        if [ -s {input.candidate_snps} ]; then
            java -Xmx{resources.mem_mb}m -jar {input.gowinda} \
                --snp-file {input.total_snps} \
                --candidate-snp-file {input.candidate_snps} \
                --gene-set-file {input.go2gene} \
                --annotation-file {input.gtf} \
                --simulations 1000000 \
                --min-significance 1 \
                --gene-definition gene \
                --threads {resources.cpus} \
                --output-file {output.enrichment} \
                --mode gene \
                --min-genes 1 >> {log} 2>&1 || true
        else
            echo "no candidate SNPs for {wildcards.subset} -- skipping Gowinda" >> {log}
        fi
 
        sed -i '1iGO_ID\tavg_genes_sim\tgenes_found\tp_value\tp_adjusted\tgenes_uniq\tgenes_max\tgenes_total\tdescription\tgene_list' {output.enrichment} 2>> {log}
        """
 
 
rule plot_gowinda_enrichment_comparisons:
    input:
        enrichment=rules.enrichment_gowinda_comparisons.output.enrichment,
    output:
        plot="results/comparisons/jaccard/gowinda/{category}.{tool}.{method}.{dataset_a}_vs_{dataset_b}.{subset}.gowinda.enrichment.png",
    params:
        title="{dataset_a} vs {dataset_b}: {subset} genes ({method}, {category})",
    resources:
        mem_mb=8000,
    log:
        "logs/comparisons/plot_gowinda_enrichment_comparisons.{category}.{tool}.{method}.{dataset_a}_vs_{dataset_b}.{subset}.log",
    script:
        "workflow/scripts/visualization/plot_gowinda_enrichment.py"

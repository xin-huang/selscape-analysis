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


import pandas as pd


refgene_txt = "resources/tools/annovar/hg38_db/hg38_refGene.txt"
gene_list = "results/positive_selection/scikit-allel/Human/1pop/YRI/moving_tajima_d/100_1/YRI.moving_tajima_d.top_0.005.candidate.genes"
out_bed = "mtjd.genes.bed"

# Read gene list
genes = pd.read_csv(gene_list, sep="\t")["Gene"].astype(str)
gene_set = set(g.strip() for g in genes.tolist() if g.strip())

# Read ANNOVAR refGene (genePred)
rg = pd.read_csv(
    refgene_txt,
    sep="\t",
    header=None,
    usecols=[2, 4, 5, 12],
    names=["chr", "start", "end", "gene"],
    dtype={"chr": str, "start": int, "end": int, "gene": str},
)

rg["gene"] = rg["gene"].astype(str).str.strip()
rg["chr"] = rg["chr"].astype(str).str.strip()

# Keep only chr1..chr22 (drop chrX, chrY, chrM, and alternative contigs)
autosomes = {f"chr{i}" for i in range(1, 23)}
rg = rg[rg["chr"].isin(autosomes)].copy()

rg = rg[rg["gene"].isin(gene_set)].copy()

bed = (
    rg.groupby(["chr", "gene"], as_index=False)
        .agg(start=("start", "min"), end=("end", "max"))
        .sort_values(["chr", "start"])
)
bed = bed[["chr", "start", "end", "gene"]]

bed.to_csv(out_bed, sep="\t", header=False, index=False)

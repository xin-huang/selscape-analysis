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


import numpy as np

include: "download_resources.smk"

rule download_1kg_low_cov_hg38_data:
    input:
        # 1KG low cov (hg38)
        expand(
            "resources/data/Human/1KG_low_cov_hg38/chr{i}.vcf.gz", i=np.arange(1, 23)
        ),
        expand(
            "resources/data/Human/ancestral_alleles/homo_sapiens_ancestor_GRCh38_bare/homo_sapiens_ancestor.{i}.bed.gz",
            i=np.arange(1, 23),
        ),
        "resources/data/Human/1KG_low_cov_hg38/samples.txt",


# 1kg low cov hg38


rule download_1KG_low_cov_hg38_vcf:
    output:
        vcf="resources/data/Human/1KG_low_cov_hg38/chr{i}.vcf.gz",
        idx="resources/data/Human/1KG_low_cov_hg38/chr{i}.vcf.gz.tbi",
    shell:
        """
        wget -c https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr{wildcards.i}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz -O {output.vcf}
        wget -c https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr{wildcards.i}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz.tbi -O {output.idx}
        """


rule create_1KG_low_cov_hg38_metadata:
    input:
        panel=rules.download_1KG_info.output.panel,
    output:
        samples="resources/data/Human/1KG_low_cov_hg38/samples.txt",
    shell:
        r"""
        sed '1d' {input.panel} | awk '{{print $1"\t"$2}}' | sed '1iSample\tPopulation' > {output.samples}
        """


rule extract_anc_info_hg38_bare:
    input:
        anc_alleles=rules.download_ensembl_ancestral_alleles_hg38.output.anc_alleles,
    output:
        bed=temp(
            "resources/data/Human/ancestral_alleles/homo_sapiens_ancestor_GRCh38_bare/homo_sapiens_ancestor.{i}.bed"
        ),
    params:
        fasta="resources/data/Human/ancestral_alleles/homo_sapiens_ancestor_GRCh38/homo_sapiens_ancestor_{i}.fa",
    run:
        import pysam
        import re

        fasta = pysam.FastaFile(params.fasta)
        with open(output.bed, "wt") as out:
            for raw_chrom in fasta.references:
                match = re.search(r"GRCh\d+:(\d+|X|Y)", raw_chrom)
                if not match:
                    print(f"Skipping unrecognized chromosome name: {raw_chrom}")
                    continue
                chrom = match.group(1)
                seq = fasta.fetch(raw_chrom).upper()
                for pos, base in enumerate(seq):
                    if base in "ACGT":
                        out.write(f"{chrom}\t{pos}\t{pos+1}\t{base}\n")
        fasta.close()


rule compress_anc_info_hg38_bare:
    input:
        bed=rules.extract_anc_info_hg38_bare.output.bed,
    output:
        bed="resources/data/Human/ancestral_alleles/homo_sapiens_ancestor_GRCh38_bare/homo_sapiens_ancestor.{i}.bed.gz",
    shell:
        """
        bgzip -c {input.bed} > {output.bed}
        tabix -p bed {output.bed}
        """

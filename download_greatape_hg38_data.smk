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


import numpy as np


rule download_greatape_hg38_data:
    input:
        # Great ape
        expand(
            "resources/data/greatape/Gorilla/chr{i}.filteranno.vcf.gz",
            i=np.arange(1, 23),
        ),
        expand(
            "resources/data/greatape/Pan/chr{i}.filteranno.vcf.gz", i=np.arange(1, 23)
        ),
        expand(
            "resources/data/greatape/Pongo/chr{i}.filteranno.vcf.gz",
            i=np.arange(1, 23),
        ),
        "resources/data/greatape/metadata.txt",
        expand(
            "resources/data/greatape/rhemac10_anc_alleles/hg38.chr{i}.bed.gz",
            i=np.arange(1, 23),
        ),
        "resources/data/greatape/refgenomes/hg38.fa",
        "resources/data/greatape/refgenomes/rheMac10.fa",
        "resources/data/greatape/refgenomes/hg38ToRheMac10.over.chain",


# great ape


rule download_greatape_vcf:
    output:
        vcf="resources/data/greatape/{species}/chr{i}.filteranno.vcf.gz",
        idx="resources/data/greatape/{species}/chr{i}.filteranno.vcf.gz.tbi",
    shell:
        """
        wget -c https://phaidra.univie.ac.at/pfsa/o_2066302/merged_segregating/{wildcards.species}/{wildcards.species}_all_filtered/chr{wildcards.i}.filteranno.vcf.gz -O {output.vcf}
        bcftools index -t {output.vcf}
        """


rule download_greatape_metadata:
    output:
        raw=temp("resources/data/greatape/metadata_full.txt"),
    shell:
        """
        wget -c https://phaidra.univie.ac.at/pfsa/o_2066302/metadata_full.txt -O {output.raw}
        """


rule create_greatape_metadata:
    input:
        metadata=rules.download_greatape_metadata.output.raw,
    output:
        metadata="resources/data/greatape/metadata.txt",
    shell:
        """
        grep -v captive {input.metadata} | awk 'NR>1 {{print $4"\t"$2}}' > {output.metadata}
        sed -i '1iSample\tPopulation' {output.metadata}
        """


rule download_greatape_reference_genomes:
    output:
        hg38="resources/data/greatape/refgenomes/hg38.fa",
        rheMac10="resources/data/greatape/refgenomes/rheMac10.fa",
        chain="resources/data/greatape/refgenomes/hg38ToRheMac10.over.chain",
    params:
        dir="resources/data/greatape/refgenomes",
    shell:
        """
        wget -c https://hgdownload.soe.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz -O {params.dir}/hg38.fa.gz
        gzip -d {params.dir}/hg38.fa.gz
        wget -c https://hgdownload.soe.ucsc.edu/goldenPath/rheMac10/bigZips/rheMac10.fa.gz -O {params.dir}/rheMac10.fa.gz
        gzip -d {params.dir}/rheMac10.fa.gz
        wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToRheMac10.over.chain.gz -O {params.dir}/hg38ToRheMac10.over.chain.gz
        gzip -d {params.dir}/hg38ToRheMac10.over.chain.gz
        """


rule extract_anc_info_greatape:
    input:
        hg38=rules.download_greatape_reference_genomes.output.hg38,
        rheMac10=rules.download_greatape_reference_genomes.output.rheMac10,
        chain=rules.download_greatape_reference_genomes.output.chain,
    output:
        anc_alleles="resources/data/greatape/rhemac10_anc_alleles/hg38.chr{i}.bed.gz",
        index="resources/data/greatape/rhemac10_anc_alleles/hg38.chr{i}.bed.gz.tbi",
    params:
        out="resources/data/greatape/rhemac10_anc_alleles/hg38.chr{i}.bed",
    resources:
        mem_mb=128000,
    shell:
        """
        python scripts/get_ancestral_info.py \
            --src-fasta {input.hg38} \
            --tgt-fasta {input.rheMac10} \
            --liftover-chain {input.chain} \
            --chr-name chr{wildcards.i} \
            --output {params.out}
        bgzip -c {params.out} > {output.anc_alleles}
        tabix -p bed {output.anc_alleles}
        rm {params.out}
        """

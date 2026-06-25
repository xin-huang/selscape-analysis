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


rule download_resources:
    input:
        directory("resources/tools/selscape"),
        # 1KG high cov (hg38)
        expand(
            "resources/data/Human/1KG_high_cov_hg38/full_chr{i}.vcf.gz", i=np.arange(1, 23)
        ),
        expand(
            "resources/data/Human/ancestral_alleles/homo_sapiens_ancestor_GRCh38/homo_sapiens_ancestor.chr{i}.bed.gz",
            i=np.arange(1, 23),
        ),
        "resources/data/Human/1KG_high_cov_hg38/samples.txt",
        "resources/data/Human/annotation/Human.hg38.gtf.gz",
        "resources/data/Human/annotation/gene2go.gz",
        "resources/data/Human/repeats/hg38.rmsk.autosomes.bed",
        "resources/data/Human/repeats/hg38.seg.dups.autosomes.bed",
        "resources/data/Human/repeats/hg38.simple.repeats.autosomes.bed",
        # 1KG low cov (hg19)
        expand("resources/data/Human/1KG_low_cov_hg19/{i}.vcf.gz", i=np.arange(1, 23)),
        expand(
            "resources/data/Human/ancestral_alleles/homo_sapiens_ancestor_GRCh37/homo_sapiens_ancestor.{i}.bed.gz",
            i=np.arange(1, 23),
        ),
        "resources/data/Human/1KG_low_cov_hg19/samples.txt",
        "resources/data/Human/annotation/Human.hg19.gtf.gz",
        "resources/data/Human/repeats/hg19.rmsk.autosomes.bed",
        "resources/data/Human/repeats/hg19.seg.dups.autosomes.bed",
        "resources/data/Human/repeats/hg19.simple.repeats.autosomes.bed",
        # 1KG low cov (hg38)
        expand(
            "resources/data/Human/1KG_low_cov_hg38/chr{i}.vcf.gz", i=np.arange(1, 23)
        ),
        "resources/data/Human/1KG_low_cov_hg38/samples.txt",
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


rule download_selscape:
    output:
        dir=directory("resources/tools/selscape"),
    shell:
        """
        if [ -d {output.dir} ]; then
            rm -rf {output.dir}
        fi
        git clone https://github.com/xin-huang/selscape {output.dir}
        ln -sfn resources/tools/selscape/workflow workflow
        """


# 1kg high cov


rule download_1KG_high_cov_vcf:
    output:
        vcf="resources/data/Human/1KG_high_cov_hg38/full_chr{i}.vcf.gz",
        idx="resources/data/Human/1KG_high_cov_hg38/full_chr{i}.vcf.gz.tbi",
    shell:
        """
        wget -c https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20190425_NYGC_GATK/CCDG_13607_B01_GRM_WGS_2019-02-19_chr{wildcards.i}.recalibrated_variants.vcf.gz -O {output.vcf}
        wget -c https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20190425_NYGC_GATK/CCDG_13607_B01_GRM_WGS_2019-02-19_chr{wildcards.i}.recalibrated_variants.vcf.gz.tbi -O {output.idx}
        """


rule download_1KG_info:
    output:
        panel="resources/data/Human/integrated_call_samples_v3.20130502.ALL.panel",
    shell:
        """
        wget -c https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel -O {output.panel}
        """


rule create_1KG_high_cov_metadata:
    input:
        panel=rules.download_1KG_info.output.panel,
    output:
        samples="resources/data/Human/1KG_high_cov_hg38/samples.txt",
    shell:
        r"""
        sed '1d' {input.panel} | awk '{{print $1"\t"$2}}' | sed '1iSample\tPopulation' > {output.samples}
        """


rule download_ncbi_annotation_hg38:
    output:
        gtf="resources/data/Human/annotation/Human.hg38.gtf.gz",
        gene2go="resources/data/Human/annotation/gene2go.gz",
    shell:
        """
        wget -c https://ftp.ncbi.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/GCF_000001405.40-RS_2024_08/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz -O {output.gtf}
        wget -c https://ftp.ncbi.nih.gov/gene/DATA/gene2go.gz -O {output.gene2go}
        """


rule download_ensembl_ancestral_alleles_hg38:
    output:
        anc_alleles=temp(
            "resources/data/Human/ancestral_alleles/homo_sapiens_ancestor_GRCh38.tar.gz"
        ),
    shell:
        """
        wget -c https://ftp.ensembl.org/pub/release-115/fasta/ancestral_alleles/homo_sapiens_ancestor_GRCh38.tar.gz -O {output.anc_alleles}
        mkdir -p resources/data/Human/ancestral_alleles/homo_sapiens_ancestor_GRCh38
        tar -xvzf {output.anc_alleles} -C resources/data/Human/ancestral_alleles
        """


rule extract_anc_info_hg38:
    input:
        anc_alleles=rules.download_ensembl_ancestral_alleles_hg38.output.anc_alleles,
    output:
        bed=temp(
            "resources/data/Human/ancestral_alleles/homo_sapiens_ancestor_GRCh38/homo_sapiens_ancestor.chr{i}.bed"
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
                        out.write(f"chr{chrom}\t{pos}\t{pos+1}\t{base}\n")
        fasta.close()


rule compress_anc_info_hg38:
    input:
        bed=rules.extract_anc_info_hg38.output.bed,
    output:
        bed="resources/data/Human/ancestral_alleles/homo_sapiens_ancestor_GRCh38/homo_sapiens_ancestor.chr{i}.bed.gz",
    shell:
        """
        bgzip -c {input.bed} > {output.bed}
        tabix -p bed {output.bed}
        """


rule download_hg38_repeat_files:
    output:
        rmsk="resources/data/Human/repeats/hg38.rmsk.txt.gz",
        segdup="resources/data/Human/repeats/hg38.genomicSuperDups.txt.gz",
        simrep="resources/data/Human/repeats/hg38.simpleRepeat.txt.gz",
    shell:
        """
        wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz -O {output.rmsk}
        wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz -O {output.segdup}
        wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/simpleRepeat.txt.gz -O {output.simrep}
        """


rule convert_hg38_repeat_files:
    input:
        rmsk=rules.download_hg38_repeat_files.output.rmsk,
        segdup=rules.download_hg38_repeat_files.output.segdup,
        simrep=rules.download_hg38_repeat_files.output.simrep,
    output:
        rmsk="resources/data/Human/repeats/hg38.rmsk.autosomes.bed",
        segdup="resources/data/Human/repeats/hg38.seg.dups.autosomes.bed",
        simrep="resources/data/Human/repeats/hg38.simple.repeats.autosomes.bed",
    shell:
        r"""
        zcat {input.rmsk} | awk 'BEGIN{{OFS="\t"}}$6!~/chr(X|Y|Un|M|[0-9]_|[0-9][0-9]_)/{{print $6,$7,$8,$11,$2,$10}}' | sed 's/^chr//' | sort -k1,1n -k2,2n -k3,3n > {output.rmsk}
        zcat {input.segdup} | awk 'BEGIN{{OFS="\t"}}$2!~/chr(X|Y|Un|M|[0-9]_|[0-9][0-9]_)/{{print $2,$3,$4,$5,$6,$7}}' | sed 's/^chr//' | sort -k1,1n -k2,2n -k3,3n > {output.segdup}
        zcat {input.simrep} | awk 'BEGIN{{OFS="\t"}}$2!~/chr(X|Y|Un|M|[0-9]_|[0-9][0-9]_)/{{print $2,$3,$4,$5,$11}}' | sed 's/^chr//' | sort -k1,1n -k2,2n -k3,3n > {output.simrep}
        """


# 1kg low cov


rule download_1KG_low_cov_vcf:
    output:
        vcf="resources/data/Human/1KG_low_cov_hg19/{i}.vcf.gz",
        idx="resources/data/Human/1KG_low_cov_hg19/{i}.vcf.gz.tbi",
    shell:
        """
        wget -c https://ftp.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.chr{wildcards.i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz -O {output.vcf}
        wget -c https://ftp.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.chr{wildcards.i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi -O {output.idx}
        """


rule create_1KG_low_cov_metadata:
    input:
        panel=rules.download_1KG_info.output.panel,
    output:
        samples="resources/data/Human/1KG_low_cov_hg19/samples.txt",
    shell:
        r"""
        sed '1d' {input.panel} | awk '{{print $1"\t"$2}}' | sed '1iSample\tPopulation' > {output.samples}
        """


rule download_ncbi_annotation_hg19:
    output:
        gtf="resources/data/Human/annotation/Human.hg19.gtf.gz",
    shell:
        """
        wget -c https://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz -O {output.gtf}
        """


rule download_ensembl_ancestral_alleles_hg19:
    output:
        anc_alleles=temp(
            "resources/data/Human/ancestral_alleles/homo_sapiens_ancestor_GRCh37_e71.tar.bz2"
        ),
    shell:
        """
        wget -c https://ftp.ensembl.org/pub/release-75/fasta/ancestral_alleles/homo_sapiens_ancestor_GRCh37_e71.tar.bz2 -O {output.anc_alleles}
        mkdir -p resources/data/Human/ancestral_alleles/homo_sapiens_ancestor_GRCh37
        tar -xjf {output.anc_alleles} -C resources/data/Human/ancestral_alleles
        """


rule extract_anc_info_hg19:
    input:
        anc_alleles=rules.download_ensembl_ancestral_alleles_hg19.output.anc_alleles,
    output:
        bed=temp(
            "resources/data/Human/ancestral_alleles/homo_sapiens_ancestor_GRCh37/homo_sapiens_ancestor.{i}.bed"
        ),
    params:
        fasta="resources/data/Human/ancestral_alleles/homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor_{i}.fa",
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


rule compress_anc_info_hg19:
    input:
        bed=rules.extract_anc_info_hg19.output.bed,
    output:
        bed="resources/data/Human/ancestral_alleles/homo_sapiens_ancestor_GRCh37/homo_sapiens_ancestor.{i}.bed.gz",
    shell:
        """
        bgzip -c {input.bed} > {output.bed}
        tabix -p bed {output.bed}
        """


rule download_hg19_repeat_files:
    output:
        rmsk="resources/data/Human/repeats/hg19.rmsk.txt.gz",
        segdup="resources/data/Human/repeats/hg19.genomicSuperDups.txt.gz",
        simrep="resources/data/Human/repeats/hg19.simpleRepeat.txt.gz",
    shell:
        """
        wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz -O {output.rmsk}
        wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/genomicSuperDups.txt.gz -O {output.segdup}
        wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/simpleRepeat.txt.gz -O {output.simrep}
        """


rule convert_hg19_repeat_files:
    input:
        rmsk=rules.download_hg19_repeat_files.output.rmsk,
        segdup=rules.download_hg19_repeat_files.output.segdup,
        simrep=rules.download_hg19_repeat_files.output.simrep,
    output:
        rmsk="resources/data/Human/repeats/hg19.rmsk.autosomes.bed",
        segdup="resources/data/Human/repeats/hg19.seg.dups.autosomes.bed",
        simrep="resources/data/Human/repeats/hg19.simple.repeats.autosomes.bed",
    shell:
        r"""
        zcat {input.rmsk} | awk 'BEGIN{{OFS="\t"}}$6!~/chr(X|Y|Un|M|[0-9]_|[0-9][0-9]_)/{{print $6,$7,$8,$11,$2,$10}}' | sed 's/^chr//' | sort -k1,1n -k2,2n -k3,3n > {output.rmsk}
        zcat {input.segdup} | awk 'BEGIN{{OFS="\t"}}$2!~/chr(X|Y|Un|M|[0-9]_|[0-9][0-9]_)/{{print $2,$3,$4,$5,$6,$7}}' | sed 's/^chr//' | sort -k1,1n -k2,2n -k3,3n > {output.segdup}
        zcat {input.simrep} | awk 'BEGIN{{OFS="\t"}}$2!~/chr(X|Y|Un|M|[0-9]_|[0-9][0-9]_)/{{print $2,$3,$4,$5,$11}}' | sed 's/^chr//' | sort -k1,1n -k2,2n -k3,3n > {output.simrep}
        """


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

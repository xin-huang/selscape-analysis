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
        # selscape
        directory("resources/tools/selscape"),
        # 1KG panel
        "resources/data/Human/1KG_info/integrated_call_samples_v3.20130502.ALL.panel",
        # annotation files
        "resources/data/Human/annotation/Human.hg38.gtf.gz",
        "resources/data/Human/annotation/gene2go.gz",
        "resources/data/Human/annotation/Human.hg19.gtf.gz",
        # repeat files
        "resources/data/Human/repeats/hg38.rmsk.autosomes.bed",
        "resources/data/Human/repeats/hg38.seg.dups.autosomes.bed",
        "resources/data/Human/repeats/hg38.simple.repeats.autosomes.bed",
        "resources/data/Human/repeats/hg19.rmsk.autosomes.bed",
        "resources/data/Human/repeats/hg19.seg.dups.autosomes.bed",
        "resources/data/Human/repeats/hg19.simple.repeats.autosomes.bed",
        # circos plots files
        "resources/data/Human/genome/hg38.chrom.sizes.bed",
        "resources/data/Human/genome/hg38.cytoBand.txt.gz",
        "resources/data/Human/genome/hg19.chrom.sizes.bed",
        "resources/data/Human/genome/hg19.cytoBand.txt.gz",


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

rule download_1KG_info:
    output:
        panel="resources/data/Human/1KG_info/integrated_call_samples_v3.20130502.ALL.panel",
    shell:
        """
        wget -c https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel -O {output.panel}
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

# annotation files
rule download_ncbi_annotation_hg38:
    output:
        gtf="resources/data/Human/annotation/Human.hg38.gtf.gz",
        gene2go="resources/data/Human/annotation/gene2go.gz",
    shell:
        """
        wget -c https://ftp.ncbi.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/GCF_000001405.40-RS_2024_08/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz -O {output.gtf}
        wget -c https://ftp.ncbi.nih.gov/gene/DATA/gene2go.gz -O {output.gene2go}
        """


rule download_ncbi_annotation_hg19:
    output:
        gtf="resources/data/Human/annotation/Human.hg19.gtf.gz",
    shell:
        """
        wget -c https://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz -O {output.gtf}
        """


# repeat files
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


# circos plots


rule download_hg38_genome_files:
    output:
        chrom_sizes=temp("resources/data/Human/genome/hg38.chrom.sizes"),
        cytoband="resources/data/Human/genome/hg38.cytoBand.txt.gz",
    shell:
        """
        wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes -O {output.chrom_sizes}
        wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz -O {output.cytoband}
        """


rule convert_hg38_chrom_sizes_to_bed:
    input:
        chrom_sizes="resources/data/Human/genome/hg38.chrom.sizes",
    output:
        bed="resources/data/Human/genome/hg38.chrom.sizes.bed",
    shell:
        r"""
        awk 'BEGIN{{OFS="\t"}}{{print $1, 0, $2}}' {input.chrom_sizes} > {output.bed}
        """


rule download_hg19_genome_files:
    output:
        chrom_sizes=temp("resources/data/Human/genome/hg19.chrom.sizes"),
        cytoband="resources/data/Human/genome/hg19.cytoBand.txt.gz",
    shell:
        """
        wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes -O {output.chrom_sizes}
        wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz -O {output.cytoband}
        """


rule convert_hg19_chrom_sizes_to_bed:
    input:
        chrom_sizes="resources/data/Human/genome/hg19.chrom.sizes",
    output:
        bed="resources/data/Human/genome/hg19.chrom.sizes.bed",
    shell:
        r"""
        awk 'BEGIN{{OFS="\t"}}{{print $1, 0, $2}}' {input.chrom_sizes} > {output.bed}
        """

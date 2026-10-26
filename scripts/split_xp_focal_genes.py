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


import os
import sys

import pandas as pd

log_fh = open(snakemake.log[0], "w")
sys.stderr = log_fh
sys.stdout = log_fh

outdir = snakemake.output.outdir
score_column = snakemake.params.score_column
pop1_sign = snakemake.params.pop1_sign
os.makedirs(outdir, exist_ok=True)


def read_tsv(path):
    try:
        return pd.read_csv(path, sep="\t")
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def bare_chrom(values):
    return values.astype(str).str.removeprefix("chr")


for scores_path, annotated_path in zip(snakemake.input.scores, snakemake.input.annotated):
    pair = os.path.basename(scores_path).split(".")[0]
    assert os.path.basename(annotated_path).startswith(pair + "."), f"{annotated_path} is not {pair}"
    pop1, pop2 = pair.split("_", 1)
    scores = read_tsv(scores_path)
    annotated = read_tsv(annotated_path)

    for focal in (pop1, pop2):
        genes = set()
        if not scores.empty and not annotated.empty:
            sign = pd.to_numeric(scores[score_column], errors="coerce")
            keep_positive = (pop1_sign == "positive") == (focal == pop1)
            hits = scores[sign > 0] if keep_positive else scores[sign < 0]

            # selscan outliers are single SNPs, delta Tajima's D outliers are windows
            starts = hits["window_start"] if "window_start" in hits else hits["BP"]
            ends = hits["window_end"] if "window_end" in hits else hits["BP"]
            intervals = {}
            for chrom, start, end in zip(bare_chrom(hits["CHR"]), starts, ends):
                intervals.setdefault(chrom, []).append((int(start), int(end)))

            for chrom, pos, gene in zip(bare_chrom(annotated["Chr"]), annotated["Start"], annotated["Gene.refGene"]):
                gene = str(gene)
                if ";" in gene:
                    continue
                if any(start <= int(pos) <= end for start, end in intervals.get(chrom, [])):
                    genes.add(gene)

        with open(os.path.join(outdir, f"{pair}.focal_{focal}.outlier.genes"), "w") as out:
            out.write("Gene\n")
            for gene in sorted(genes):
                out.write(f"{gene}\n")
        print(f"{pair} focal={focal}: {len(genes)} genes")

log_fh.close()

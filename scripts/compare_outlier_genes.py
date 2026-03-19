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


from pathlib import Path
import pandas as pd

POP_ANCESTRY = {
    "ACB": "AFR", "ASW": "AFR", "ESN": "AFR", "GWD": "AFR",
    "LWK": "AFR", "MSL": "AFR", "YRI": "AFR",
    "CLM": "AMR", "MXL": "AMR", "PEL": "AMR", "PUR": "AMR",
    "CDX": "EAS", "CHB": "EAS", "CHS": "EAS", "JPT": "EAS", "KHV": "EAS",
    "CEU": "EUR", "FIN": "EUR", "GBR": "EUR", "IBS": "EUR", "TSI": "EUR",
    "BEB": "SAS", "GIH": "SAS", "ITU": "SAS", "PJL": "SAS", "STU": "SAS",
}

SUPERPOP_NAMES = {
    "African": "AFR", "European": "EUR", "Asian": "EAS",
    "American": "AMR", "SouthAsian": "SAS",
}

ALL_POPS = list(POP_ANCESTRY.keys())

def load_genes(path):
    return {g for g in Path(path).read_text().splitlines() if g and g != "Gene"}

def build_pop_genes(file_list):
    pop_genes = {p: set() for p in ALL_POPS}
    for f in file_list:
        pop = Path(f).name.split(".")[0]
        pop_genes[pop].update(load_genes(f))
    return pop_genes

def resolve_pops(pop_token):
    if pop_token is None:
        return ALL_POPS, "ALL"
    if pop_token in POP_ANCESTRY:
        return [pop_token], pop_token
    if pop_token in SUPERPOP_NAMES:
        sp = SUPERPOP_NAMES[pop_token]
        return [p for p, a in POP_ANCESTRY.items() if a == sp], pop_token
    if pop_token in set(POP_ANCESTRY.values()):
        return [p for p, a in POP_ANCESTRY.items() if a == pop_token], pop_token
    return ALL_POPS, "ALL"

pos_genes = build_pop_genes(
    list(snakemake.input.ihs) + list(snakemake.input.nsl) +
    list(snakemake.input.mtjd_pos) + list(snakemake.input.wtjd_pos)
)
bal_genes = build_pop_genes(
    list(snakemake.input.betascan) + list(snakemake.input.mtjd_bal) +
    list(snakemake.input.wtjd_bal)
)

outdir = Path(snakemake.output.outdir)
outdir.mkdir(parents=True, exist_ok=True)

study_dir = Path(snakemake.input.study_dir)

for study_file in sorted(study_dir.rglob("*.list")):
    sel_type = "balancing" if any("balancing" in p for p in study_file.parts) else "positive"
    src = bal_genes if sel_type == "balancing" else pos_genes

    study_name = study_file.stem.split(".")[0]

    stem_parts = study_file.stem.split(".")
    pop_token = stem_parts[1] if stem_parts[1] not in ("candidate", "genes") else None

    pops, pop_label = resolve_pops(pop_token)
    study_genes = sorted(load_genes(study_file))
    our_genes = set().union(*[src[p] for p in pops])

    rows = [
        {
            "study":   study_name,
            "pop":     pop_label,
            "gene":    gene,
            "overlap": "Y" if gene in our_genes else "N",
        }
        for gene in study_genes
    ]

    out_name = study_file.stem.replace(".candidate.genes", "") + f".{sel_type}.overlap.tsv"
    pd.DataFrame(rows).to_csv(outdir / out_name, sep="\t", index=False)

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


 
from itertools import combinations
from pathlib import Path
 
import pandas as pd
 
 
def jaccard(a, b):
    union = a | b
    if len(union) == 0:
        return float("nan")
    return len(a & b) / len(union)
 
 
def genes_in(path):
    try:
        df = pd.read_csv(path, sep="\t")
    except pd.errors.EmptyDataError:
        return set()
    if df.empty:
        return set()
    return set(df.iloc[:, 0].astype(str))
 
 
def pooled_gene_sets(root_dir, parse_path):
 
    pooled = {}
    for filepath in Path(root_dir).rglob("*.outlier.genes"):
        rel_path = filepath.relative_to(root_dir).as_posix()
        meta = parse_path(rel_path)
        if meta is None:
            continue
        key = (meta["category"], meta["tool"], meta["method"], meta["dataset"])
        if key not in pooled:
            pooled[key] = set()
        pooled[key].update(genes_in(filepath))
    return pooled
 
 
def jaccard_matrix(pooled, category, tool, method, datasets):
    """Jaccard-similarity matrix for one (category, tool, method), across
    all `datasets` -- every pair is computed."""
    matrix = pd.DataFrame(float("nan"), index=datasets, columns=datasets)
 
    for dataset in datasets:
        genes = pooled.get((category, tool, method, dataset), set())
        matrix.loc[dataset, dataset] = jaccard(genes, genes)
 
    for dataset_a, dataset_b in combinations(datasets, 2):
        genes_a = pooled.get((category, tool, method, dataset_a), set())
        genes_b = pooled.get((category, tool, method, dataset_b), set())
        value = jaccard(genes_a, genes_b)
        matrix.loc[dataset_a, dataset_b] = value
        matrix.loc[dataset_b, dataset_a] = value
 
    return matrix
 
 
def population_gene_sets(root_dir, parse_path):
    """(category, tool, method, dataset, population) -> set(genes).
 
    Unlike pooled_gene_sets(), this does NOT union genes across
    populations -- each population keeps its own gene set. Use this when
    you want to compare the SAME population across different datasets
    (e.g. is YRI's candidate gene list reproducible between
    1kg_high_hg38, 1kg_low_hg38 and 1kg_low_hg19?), where pooling all 26
    populations together would wash out exactly the signal you want to
    see.
    """
    gene_sets = {}
    for filepath in Path(root_dir).rglob("*.outlier.genes"):
        rel_path = filepath.relative_to(root_dir).as_posix()
        meta = parse_path(rel_path)
        if meta is None:
            continue
        key = (meta["category"], meta["tool"], meta["method"], meta["dataset"], meta["unit"])
        gene_sets[key] = genes_in(filepath)
    return gene_sets
 
 
def overlap_and_differences(a, b):
    """Split two gene sets into (overlap, unique_to_a, unique_to_b) --
    the genes they share, and the genes that belong only to each."""
    overlap = a & b
    unique_a = a - b
    unique_b = b - a
    return overlap, unique_a, unique_b
 
 
def population_jaccard_table(gene_sets, category, tool, method, populations, datasets):
    """Table with one row per population and one column per dataset pair
    (every pair among `datasets`) for one (category, tool, method). Each
    cell is the Jaccard similarity of that population's own gene set
    between the two datasets in that column -- i.e. how reproducible that
    population's candidate genes are across dataset versions.
    """
    pairs = list(combinations(datasets, 2))
    columns = [f"{dataset_a} vs {dataset_b}" for dataset_a, dataset_b in pairs]
    table = pd.DataFrame(float("nan"), index=populations, columns=columns)
 
    for population in populations:
        for (dataset_a, dataset_b), column in zip(pairs, columns):
            genes_a = gene_sets.get((category, tool, method, dataset_a, population), set())
            genes_b = gene_sets.get((category, tool, method, dataset_b, population), set())
            table.loc[population, column] = jaccard(genes_a, genes_b)
 
    return table

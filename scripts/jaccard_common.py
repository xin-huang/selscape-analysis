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

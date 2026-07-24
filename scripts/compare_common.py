from pathlib import Path
import pandas as pd
 
 
def top_terms(path, top_n, p_threshold):
    """(GO_ID, description) pairs, sorted like plot_gowinda_enrichment.py
    (p_adjusted asc, genes_found desc), optionally pre-filtered to
    p_adjusted < p_threshold."""
    try:
        df = pd.read_csv(path, sep="\t")
    except pd.errors.EmptyDataError:
        return []
    if df.empty:
        return []
    if p_threshold is not None:
        df = df[df["p_adjusted"] < p_threshold]
    df = df.sort_values(["p_adjusted", "genes_found"], ascending=[True, False]).head(top_n)
    return list(zip(df["GO_ID"], df["description"]))
 
 
def build_hits(root_dir, top_n, p_threshold, parse_path):
    """One row per (category, tool, method, dataset, unit, GO_ID, description)."""
    rows = []
    for f in Path(root_dir).rglob("*.gowinda.enrichment.tsv"):
        meta = parse_path(f.relative_to(root_dir).as_posix())
        if meta is None:
            continue
        for go_id, desc in top_terms(f, top_n, p_threshold):
            rows.append({**meta, "GO_ID": go_id, "description": desc})
    return pd.DataFrame(rows)
 
 
def write_csv(hits, category, tool, method, dataset_a, dataset_b, out_dir):
    sub = hits[
        (hits["category"] == category) & (hits["tool"] == tool) &
        (hits["method"] == method) & (hits["dataset"].isin([dataset_a, dataset_b]))
    ].copy()
    sub["population"] = sub["dataset"] + ":" + sub["unit"]
 
    wide = sub.pivot_table(index=["GO_ID", "description"], columns="population",
                            values="unit", aggfunc="size", fill_value=0)
    wide = wide.reindex(sorted(wide.columns), axis=1)
    count = wide.sum(axis=1)
    out = (wide > 0).replace({True: "Y", False: "N"})
    out.insert(0, "count", count)
    out = out.reset_index().sort_values("count", ascending=False)
 
    tool_us = tool.replace("-", "_")
    fname = f"{category}.{tool_us}.{method}.{dataset_a}_vs_{dataset_b}.go_terms.csv"
    out.to_csv(Path(out_dir) / fname, index=False)
 
 
def compare_pair(hits, dataset_a, dataset_b, out_dir):
    """Write one CSV per (category, tool, method) shared by both datasets."""
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    if hits.empty:
        return
    sub = hits[hits["dataset"].isin([dataset_a, dataset_b])]
    keys = ["category", "tool", "method"]
    groups_a = set(map(tuple, sub[sub["dataset"] == dataset_a][keys].drop_duplicates().values))
    groups_b = set(map(tuple, sub[sub["dataset"] == dataset_b][keys].drop_duplicates().values))
    for category, tool, method in sorted(groups_a & groups_b):
        write_csv(hits, category, tool, method, dataset_a, dataset_b, out_dir)

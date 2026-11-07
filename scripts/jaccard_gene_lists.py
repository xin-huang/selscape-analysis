import sys
from pathlib import Path
 
sys.path.insert(0, snakemake.scriptdir)
from jaccard_common import pooled_gene_sets, overlap_and_differences
from compare_common import parse_positive_path, parse_balancing_path
 
log_fh = open(snakemake.log[0], "w") if snakemake.log else sys.stdout
sys.stderr = log_fh
sys.stdout = log_fh
 
 
def write_gene_list(genes, path):
    with open(path, "w") as f:
        f.write("Gene\n")
        for gene in sorted(genes):
            f.write(f"{gene}\n")
 
 
def write_gene_lists(comparisons, outdir, log_fh):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
 
    sources = [
        ("results/positive_selection", parse_positive_path),
        ("results/balancing_selection", parse_balancing_path),
    ]
 
    for root_dir, parse_path in sources:
        pooled = pooled_gene_sets(Path(root_dir), parse_path)
        if not pooled:
            print(f"no outlier.genes files found under {root_dir}, skipping", file=log_fh)
            continue
 
        combos = set()
        for category, tool, method, dataset in pooled:
            combos.add((category, tool, method))
        combos = sorted(combos)
 
        for category, tool, method in combos:
            for dataset_a, dataset_b in comparisons:
                genes_a = pooled.get((category, tool, method, dataset_a), set())
                genes_b = pooled.get((category, tool, method, dataset_b), set())
                overlap, unique_a, unique_b = overlap_and_differences(genes_a, genes_b)
 
                stem = f"{category}.{tool}.{method}.{dataset_a}_vs_{dataset_b}"
                write_gene_list(overlap, outdir / f"{stem}.overlap.genes")
                write_gene_list(unique_a, outdir / f"{stem}.unique_{dataset_a}.genes")
                write_gene_list(unique_b, outdir / f"{stem}.unique_{dataset_b}.genes")
                print(
                    f"wrote {stem}: overlap={len(overlap)}, "
                    f"unique_{dataset_a}={len(unique_a)}, unique_{dataset_b}={len(unique_b)}",
                    file=log_fh,
                )
 
 
comparisons = snakemake.params.comparisons
assert comparisons, "config['gowinda_comparisons'] is empty -- nothing to compare."
 
outdir = Path(snakemake.output.outdir)
write_gene_lists(comparisons, outdir, log_fh)

import sys
from pathlib import Path
 
sys.path.insert(0, snakemake.scriptdir)
from compare_common import build_hits, compare_pair
 
sys.stderr = open(snakemake.log[0], "w")
 
 
def parse_path(rel_path):
    parts = rel_path.split("/")
 
    if len(parts) == 7 and parts[0] == "selscan" and parts[3] in ("1pop", "2pop"):
        method = parts[5].split("_", 1)[0]  # folder is "{method}_{maf}", e.g. "ihs_0.05"
        return dict(category="positive_selection", tool="selscan", method=method,
                    dataset=parts[2], unit=parts[4])
 
    if len(parts) == 8 and parts[0] == "scikit-allel" and parts[3] in ("1pop", "2pop"):
        return dict(category="positive_selection", tool="scikit-allel", method=parts[5],
                    dataset=parts[2], unit=parts[4])
 
    return None
 
 
results_dir = Path("results/positive_selection")
hits = build_hits(results_dir, snakemake.params.top_n, snakemake.params.get("p_adjusted_threshold"), parse_path)
compare_pair(hits, snakemake.wildcards.dataset_a, snakemake.wildcards.dataset_b, snakemake.output[0])
 

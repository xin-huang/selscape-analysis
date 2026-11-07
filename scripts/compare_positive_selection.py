import sys
from pathlib import Path
 
sys.path.insert(0, snakemake.scriptdir)
from compare_common import build_hits, compare_pair, parse_positive_path as parse_path
 
sys.stderr = open(snakemake.log[0], "w")
 
results_dir = Path("results/positive_selection")
hits = build_hits(results_dir, snakemake.params.top_n, snakemake.params.get("p_adjusted_threshold"), parse_path)
compare_pair(hits, snakemake.wildcards.dataset_a, snakemake.wildcards.dataset_b, snakemake.output[0])

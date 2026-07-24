import sys
from pathlib import Path
 
sys.path.insert(0, snakemake.scriptdir)
from compare_common import build_hits, compare_pair
 
sys.stderr = open(snakemake.log[0], "w")
 
 
def parse_path(rel_path):
    """rel_path relative to results/balancing_selection/, e.g.
    'scikit-allel/Human/1kg_high_hg38/windowed_tajima_d/YRI/100000_1/YRI.windowed_tajima_d.top_0.05.gowinda.enrichment.tsv'
    'betascan/Pan/pan/PPA/m_0.15/PPA.hg38.m_0.15.b1.top_0.0005.gowinda.enrichment.tsv'
    """
    parts = rel_path.split("/")
    if len(parts) == 7 and parts[0] == "scikit-allel":
        return dict(category="balancing_selection", tool="scikit-allel", method=parts[3],
                    dataset=parts[2], unit=parts[4])
    if len(parts) == 6 and parts[0] == "betascan":
        return dict(category="balancing_selection", tool="betascan", method="betascan_b1",
                    dataset=parts[2], unit=parts[3])
    return None
 
 
results_dir = Path("results/balancing_selection")
hits = build_hits(results_dir, snakemake.params.top_n, snakemake.params.get("p_adjusted_threshold"), parse_path)
compare_pair(hits, snakemake.wildcards.dataset_a, snakemake.wildcards.dataset_b, snakemake.output[0])
 

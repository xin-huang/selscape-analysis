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

import pandas as pd
import re

# get input files from Snakemake
populations = snakemake.params.populations
ci_files = snakemake.input.ci_files
tsv_files = snakemake.input.tsv_files
output_file = snakemake.output.merged

# parse .godambe.ci file to extract lower and upper bounds
# uses last set of bounds (most refined)
def parse_ci_file(ci_file):
    with open(ci_file, 'r') as f:
        lines = f.readlines()
    
    lower_bounds = None
    upper_bounds = None
    
    # read from b ottom up to get last bounds
    for line in reversed(lines):
        if line.startswith("Lower bounds"):
            match = re.search(r'\[(.*?)\]', line)
            if match:
                lower_bounds = [float(x) for x in match.group(1).split()]
        elif line.startswith("Upper bounds"):
            match = re.search(r'\[(.*?)\]', line)
            if match:
                upper_bounds = [float(x) for x in match.group(1).split()]
        # stop once we have both bounds
        if lower_bounds and upper_bounds:
            break
    
    return lower_bounds, upper_bounds

# parse .godambe.ci.tsv file to get point estimates
# resturns last row (most refined)
def parse_tsv_file(tsv_file):
    df = pd.read_csv(tsv_file, sep=r'\s+')
    last_row = df.iloc[-1]
    return last_row['log_mu'], last_row['log_sigma'], last_row['misid']

data = []

for pop, ci_file, tsv_file in zip(populations, ci_files, tsv_files):
    lower_bounds, upper_bounds = parse_ci_file(ci_file)
    log_mu, log_sigma, misid = parse_tsv_file(tsv_file)
    # store all parameters: point estimates + lower/upper bounds    
    data.append({
        'Pop': pop,
        'mu': log_mu,              # μ (log mean of DFE)
        'mu_lb': lower_bounds[0],  # μ lower bound
        'mu_ub': upper_bounds[0],  # μ upper bound
        'sigma': log_sigma,        # σ (log std dev of DFE)
        'sigma_lb': lower_bounds[1],  # σ lower bound
        'sigma_ub': upper_bounds[1],  # σ upper bound
        'misid': misid,            # Misidentification rate
        'misid_lb': lower_bounds[2],  # Misid lower bound
        'misid_ub': upper_bounds[2]   # Misid upper bound
    })

df = pd.DataFrame(data)
df.to_csv(output_file, sep='\t', index=False)

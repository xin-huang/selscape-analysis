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


rule all_visualization:
    input:
        expand("results/plots/circos/{species}/{ppl}/{ppl}_positive_selection_circos_scores.png", species=config["species"], ppl=config["populations"]),
        expand("results/plots/circos/{species}/{ppl}/{ppl}_balancing_selection_circos_scores.png", species=config["species"], ppl=config["populations"]),
        expand("results/plots/dfe/{species}/{species}.two_epoch.lognormal.dfe_params.png", species=config["species"])


rule make_positive_selection_circos_scores:
    input:
        ihs_scores="results/positive_selection/selscan/{species}/1pop/{ppl}/ihs_0.05/{ppl}.normalized.ihs.scores",
        nsl_scores="results/positive_selection/selscan/{species}/1pop/{ppl}/nsl_0.05/{ppl}.normalized.nsl.scores",
        mtjd_scores="results/positive_selection/scikit-allel/{species}/1pop/{ppl}/moving_tajima_d/100_1/{ppl}.moving_tajima_d.scores",
        wtjd_scores="results/positive_selection/scikit-allel/{species}/1pop/{ppl}/windowed_tajima_d/100000_1/{ppl}.windowed_tajima_d.scores",
    output:
        plot="results/plots/circos/{species}/{ppl}/{ppl}_positive_selection_circos_scores.png"
    params:
        population="{ppl}",
        plot_type="positive_selection"
    resources:
        mem_gb=32,
    script:
        "scripts/plot_circos_scores.py"


rule make_balancing_selection_circos_scores:
    input:
        b1_scores="results/balancing_selection/betascan/{species}/{ppl}/m_0.15/{ppl}.hg38.m_0.15.b1.scores",
        mtjd_bal_scores="results/balancing_selection/scikit-allel/{species}/moving_tajima_d/{ppl}/100_1/{ppl}.moving_tajima_d.merged.scores",
        wtjd_bal_scores="results/balancing_selection/scikit-allel/{species}/windowed_tajima_d/{ppl}/100000_1/{ppl}.windowed_tajima_d.merged.scores",
    output:
        plot="results/plots/circos/{species}/{ppl}/{ppl}_balancing_selection_circos_scores.png"
    params:
        population="{ppl}",
        plot_type="balancing_selection"
    script:
        "scripts/plot_circos_scores.py"


rule merge_dfe_confidence_intervals:
    input:
        bestfit_files=expand("results/dadi/{species}/dfe/{ppl}/InferDFE/{ppl}.hg38.two_epoch.lognormal.InferDFE.bestfits", ppl=config["populations"], allow_missing=True),
        ci_files=expand("results/dadi/{species}/dfe/{ppl}/StatDFE/{ppl}.hg38.two_epoch.lognormal.godambe.ci", ppl=config["populations"], allow_missing=True),
    output:
        merged="results/plots/dfe/{species}/{species}.two_epoch.lognormal.dfe_params.tsv"
    params:
        populations=config["populations"],
    script:
        "scripts/merge_dfe_ci.py"


rule plot_dfe_confidence_intervals:
    input:
        data=rules.merge_dfe_confidence_intervals.output.merged
    output:
        plot="results/plots/dfe/{species}/{species}.two_epoch.lognormal.dfe_params.png"
    script:
        "scripts/plot_dfe_params.py"

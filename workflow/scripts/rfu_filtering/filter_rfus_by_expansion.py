"""
Snakemake script for compiling RFU filtering results, whereby RFUs are filtered
by
- TCR size (in order for a TCR to be counted towards an RFU)
- RFU expansion in terms of minimum number of individuals with at least a given
  level of RFU expression
"""

import pandas as pd
import statsmodels.stats.multitest


def main(
    rfu_stats_fname,
    rfu_counts_fname,
    min_expansion,
    min_expanded_indivs,
    fdr_cutoff,
    output_tsv,
    output_rfu_stats,
):
    # Read in inputs.
    rfu_stats = pd.read_csv(rfu_stats_fname, header=0, sep="\t", index_col=0)
    rfu_counts = pd.read_csv(rfu_counts_fname, sep="\t", header=0, index_col=0)

    # Select RFUs that exceed the minimum expansion cutoff.
    expanded_rfus = rfu_counts.index[
        (rfu_counts >= min_expansion).sum(1) >= min_expanded_indivs
    ]
    rfu_stats_filtered = rfu_stats.loc[
        rfu_stats.index.isin(expanded_rfus)
    ].copy()

    # Compute RFUs tested and passing restricted testing.
    n_tested_rfus = rfu_stats_filtered.shape[0]
    tested_variables = rfu_stats_filtered.columns[
        rfu_stats_filtered.columns.str.endswith(".fdr")
    ].str.replace(".fdr$", "", regex=True)
    n_fdr_pass_rfus = pd.Series({"n_tested_rfus": n_tested_rfus})
    if n_tested_rfus > 0:
        for covariate_label in tested_variables:
            rfu_stats_filtered[covariate_label + ".fdr"] = (
                statsmodels.stats.multitest.multipletests(
                    rfu_stats_filtered[covariate_label + ".pvalue"].fillna(1.0),
                    method="fdr_bh",
                )[1]
            )
            n_fdr_pass_rfus["n_fdr_pass_rfus." + covariate_label] = sum(
                rfu_stats_filtered[covariate_label + ".fdr"] <= fdr_cutoff
            )

    # Output results
    pd.DataFrame([n_fdr_pass_rfus]).to_csv(
        output_tsv, header=True, index=False, sep="\t"
    )
    rfu_stats_filtered.index.name = "centroid"
    rfu_stats_filtered.to_csv(
        output_rfu_stats, header=True, index=True, sep="\t"
    )


if __name__ == "__main__":
    main(
        rfu_stats_fname=snakemake.input.rfu_stats,
        rfu_counts_fname=snakemake.input.rfu_counts,
        min_expansion=snakemake.params.min_expansion,
        min_expanded_indivs=snakemake.params.min_expanded_indivs,
        fdr_cutoff=snakemake.params.fdr_cutoff,
        output_tsv=snakemake.output.n_signif_hits,
        output_rfu_stats=snakemake.output.filtered_rfus,
    )

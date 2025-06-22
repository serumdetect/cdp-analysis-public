"""Snakemake script to filter RFUs for blinded covariates."""

import pandas as pd


def main(
    input_tsv,
    full_model_rfu_stats_tsv,
    fdr_cutoff,
    blinded_covariates,
    output_tsv,
):
    # Load data.
    rfu_df = pd.read_table(input_tsv, index_col=0)
    full_model_rfu_stats_tsv = pd.read_table(
        full_model_rfu_stats_tsv, index_col=0
    )

    # Sanity check: RFU filtering shouldn't be affected by covariate blinding.
    # Hence, all RFUs should be present in the RFU status table.
    assert rfu_df.index.equals(full_model_rfu_stats_tsv.index)

    dropped_covariates_fdrs = full_model_rfu_stats_tsv[
        [f"{covariate}.fdr" for covariate in blinded_covariates]
    ]
    keep_features = (dropped_covariates_fdrs >= fdr_cutoff).all(axis=1)
    rfu_df_filtered = rfu_df.loc[keep_features]

    # Output results.
    rfu_df_filtered.to_csv(output_tsv, sep="\t", index=True)


if __name__ == "__main__":
    main(
        input_tsv=snakemake.input.input_tsv,
        full_model_rfu_stats_tsv=snakemake.input.full_model_rfu_stats_tsv,
        fdr_cutoff=snakemake.params.fdr_cutoff,
        blinded_covariates=snakemake.params.blinded_covariates,
        output_tsv=snakemake.output.output_tsv,
    )

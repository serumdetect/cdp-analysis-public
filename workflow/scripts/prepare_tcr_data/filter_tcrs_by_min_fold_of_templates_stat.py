"""
Snakemake script for filtering TCR data by template count above a certain fold
of the reference quantile or the trimmed mean.
"""

import functools

import loguru
import numpy as np
import pandas as pd
import scipy.stats


def main(
    input_tcrs,
    min_fold_of_stat,
    reference_quantile,
    v_gene_specific,
    output_tcrs,
):
    # Read inputs.
    tcrs = pd.read_feather(input_tcrs)

    # Calculate the template count at reference quantile of each sample.
    if v_gene_specific:
        group_cols = ["sample_name", "v_gene"]
    else:
        group_cols = ["sample_name"]
    if reference_quantile is not None:
        reference_templates_of_sample = tcrs.groupby(group_cols, observed=True)[
            "templates"
        ].quantile(reference_quantile)
    else:
        # Use quartile-trimmed mean instead. Note, this is computed here even
        # though this information already exists in the metadata, in order to
        # be also compatible with V-gene specific filtering.
        reference_templates_of_sample = tcrs.groupby(group_cols, observed=True)[
            "templates"
        ].apply(functools.partial(scipy.stats.trim_mean, proportiontocut=0.25))

    # Filter TCRs that have `min_fold_of_stat` times the reference quantile
    # template count at minimum.
    tcrs = (
        tcrs.merge(
            np.floor(
                reference_templates_of_sample * min_fold_of_stat
            ).rename("min_templates"),
            left_on=group_cols,
            right_index=True,
        )
        .query("templates >= min_templates")
        .drop(columns="min_templates")
    )

    # Output results.
    tcrs.reset_index(drop=True).to_feather(output_tcrs)


if __name__ == "__main__":
    # Set loguru to output to snakemake.log.stderr.
    loguru.logger.add(snakemake.log.stderr)
    main(
        input_tcrs=snakemake.input[0],
        min_fold_of_stat=snakemake.params.min_fold_of_stat,
        reference_quantile=snakemake.params.reference_quantile,
        v_gene_specific=snakemake.params.v_gene_specific,
        output_tcrs=snakemake.output[0],
    )

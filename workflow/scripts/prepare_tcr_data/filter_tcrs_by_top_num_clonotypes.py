"""Snakemake script for filtering TCR data by top N clonotypes."""

import loguru
import pandas as pd


def filter_by_raw_ranks(tcrs, top_n):
    """Filter TCR data by the top N clonotypes per sample overall."""
    # Keep top N clonotypes per sample.
    tcrs = (
        tcrs.sort_values("templates", ascending=False)
        .groupby("sample_name")
        .head(top_n)
    )

    return tcrs


def filter_by_normalized_ranks(tcrs, top_n, reference_quantile):
    """
    Filter TCR data by the top N clonotypes computed by using normalized counts
    calculated within each V gene.
    """
    # Normalize template counts within each V gene to the reference rank.
    tcrs["v_normed_ranks"] = tcrs.groupby(["sample_name", "v_gene"])[
        "templates"
    ].transform(lambda x: x.rank(method="first", ascending=False, pct=True))

    # Keep top N clonotypes per sample.
    tcrs = tcrs.sort_values("v_normed_ranks").groupby("sample_name").head(top_n)

    return tcrs


def filter_by_normalized_templates(tcrs, top_n, reference_quantile):
    """
    Filter TCR data by the top N clonotypes computed by normalizing raw
    templates with the reference template quantile of that V gene.
    """
    # Normalize template counts within each V gene to the reference quantile.
    tcrs["v_normed_templates"] = tcrs.groupby(["sample_name", "v_gene"])[
        "templates"
    ].transform(lambda x: x / x.quantile(reference_quantile))

    # Keep top N clonotypes per sample.
    tcrs = (
        tcrs.sort_values("v_normed_templates", ascending=False)
        .groupby("sample_name")
        .head(top_n)
    )

    return tcrs


def main(
    input_tcrs, top_n, filtering_strategy, reference_quantile, output_tcrs
):
    # Read inputs.
    tcrs = pd.read_feather(input_tcrs)

    # Filter out samples with fewer than the top N clonotypes.
    num_clonotypes = tcrs["sample_name"].value_counts()
    samples_to_remove = num_clonotypes.index[num_clonotypes < top_n]
    loguru.logger.warning(
        f"Removing {len(samples_to_remove)}/{len(num_clonotypes)} samples with"
        f" fewer than {top_n} clonotypes."
    )
    tcrs = tcrs.loc[~tcrs["sample_name"].isin(samples_to_remove)]

    # Shuffle the data so that any ties in the ranking are broken randomly.
    tcrs = tcrs.sample(frac=1.0, random_state=0)
    if filtering_strategy == "raw_ranks":
        tcrs = filter_by_raw_ranks(tcrs, top_n)
    elif filtering_strategy == "normalized_ranks":
        tcrs = filter_by_normalized_ranks(tcrs, top_n, reference_quantile)
    elif filtering_strategy == "normalized_templates":
        tcrs = filter_by_normalized_templates(tcrs, top_n, reference_quantile)

    # Sanity check that the number of clonotypes per sample is equal to top_n.
    assert all(tcrs.groupby("sample_name", observed=True).size() == top_n)

    # Output results.
    tcrs.reset_index(drop=True).to_feather(output_tcrs)


if __name__ == "__main__":
    # Set loguru to output to snakemake.log.stderr.
    loguru.logger.add(snakemake.log.stderr)
    main(
        snakemake.input[0],
        snakemake.params.top_n,
        snakemake.wildcards.filtering_strategy,
        snakemake.params.reference_quantile,
        snakemake.output[0],
    )

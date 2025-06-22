"""
Compute query RFU counts taking into account the minimum size of each TCR.

An "expected sample substitution" scheme is used to compute the expected effect
of each query TCR being substituted with one of the reference TCRs contributing
to an RFU. This expected value is used as a weight to weigh each RFU count.
"""

import gzip

import loguru
import polars as pl  # noqa: E402
import tqdm


def _sanity_check_samples(query_samples, ref_samples):
    """Sanity check that there are no shared reference and query samples."""
    assert len(set(query_samples).intersection(ref_samples)) == 0


def _compute_query_tcr_weights(query_tcrs, ref_tcrs, min_size, num_ref_samples):
    """
    Compute the weight of each query TCR based on the number of reference
    samples that would bring the final TCR size below the cutoff.
    """
    # Prepare query centroid and sanity check that they are all amongst the
    # reference centroids.
    query_centroids = query_tcrs["centroid"].unique().to_list()
    assert all(
        [centroid in ref_tcrs["centroid"] for centroid in query_centroids]
    )
    query_tcrs = query_tcrs.lazy()
    ref_tcrs = ref_tcrs.lazy()

    # Simplify by throwing away irrelevant reference RFUs.
    ref_tcrs = ref_tcrs.filter(pl.col("centroid").is_in(query_centroids))

    # Compute the RFU count of each RFU including the new query sample. Only
    # keep RFUs that would exceed minimum size in the first place.
    combined_tcrs = pl.concat(
        [
            query_tcrs.with_columns(pl.lit(True).alias("is_query_sample")),
            ref_tcrs.with_columns(pl.lit(False).alias("is_query_sample")),
        ]
    )
    tcr_sizes = (
        combined_tcrs.group_by(["v_gene", "cdr3"])
        .len()
        .rename({"len": "total_size"})
    ).filter(pl.col("total_size") >= min_size)

    # Compute the number of clonotypes each reference sample contributes to the
    # overall TCR size.
    ref_tcr_sizes = (
        ref_tcrs.group_by(["sample_name", "v_gene", "cdr3"])
        .len()
        .rename({"len": "ref_sample_size"})
    ).join(tcr_sizes, on=["v_gene", "cdr3"], how="inner")

    # For each TCR, compute the proportion of reference samples that would bring
    # the final TCR size below the cutoff.
    #
    # We do this by first only keeping samples that are needed to keep TCR size
    # above the cutoff.
    num_critical_ref_samples = (
        ref_tcr_sizes.filter(
            (pl.col("total_size") - pl.col("ref_sample_size") < min_size)
        )
        .group_by(["v_gene", "cdr3"])
        .len()
        .rename({"len": "num_critical_ref_samples"})
    )

    # Assign a weight to each query TCR based on the number of reference samples
    # that would bring the final TCR size below the cutoff.
    #
    # We first perform an inner join with `tcr_sizes`, because this gets rid of
    # any TCRs that would not exceed the `min_size` cutoff even with the query
    # sample added in.
    query_tcrs_with_weights = (
        query_tcrs.join(tcr_sizes, on=["v_gene", "cdr3"], how="inner")
        .join(
            num_critical_ref_samples,
            on=["v_gene", "cdr3"],
            how="left",
            coalesce=True,
        )
        .with_columns(pl.col("num_critical_ref_samples").fill_null(0))
        .with_columns(
            (
                (num_ref_samples - pl.col("num_critical_ref_samples"))
                / num_ref_samples
            ).alias("tcr_weight")
        )
    )

    return query_tcrs_with_weights


def filtered_query_rfu_count(query_tcrs, ref_tcrs, min_size, num_ref_samples):
    """Compute RFU counts after filtering query TCRs by a given minimum size.

    This function implements the sample substitution whereby the given query
    sample is "substituted" with one of the reference samples. After each
    substitution, TCRs with sufficient size are tallied to give RFU counts.
    Finally, the RFU count is computed as the expected value of the RFU counts
    across all possible substitutions.
    """
    query_tcrs_with_weights = _compute_query_tcr_weights(
        query_tcrs, ref_tcrs, min_size, num_ref_samples
    )
    rfu_counts = (
        query_tcrs_with_weights.group_by("centroid")
        .agg(pl.col("tcr_weight").sum().alias("tcr_count"))
        .collect()
    )
    return rfu_counts


def main(query_tcrs_ftr, ref_tcrs_ftr, min_size, is_cv_test, rfu_counts_tsv):
    # Read the query and reference TCRs.
    loguru.logger.info("Reading inputs...")
    ref_tcrs = pl.read_ipc(ref_tcrs_ftr)
    query_tcrs = pl.read_ipc(query_tcrs_ftr)

    # Sanity check that there are no shared samples between the query and
    # reference TCRs.
    query_samples = query_tcrs["sample_name"].unique().to_list()
    ref_samples = ref_tcrs["sample_name"].unique().to_list()
    if is_cv_test:
        # If this is CV test split, then ensure that the query samples are
        # distinct from the reference samples.
        _sanity_check_samples(query_samples, ref_samples)

    # Calculate the RFU counts with each TCR's contribution weighted by the
    # probability that this TCR would fall short of the minimum size cutoff if
    # the query sample were replaced with one of the reference samples at
    # random.
    #
    # If min_size == 1, then do not bother with the above and simply tally the
    # number of TCRs in each centroid (row) in each sample (column).
    if min_size > 1:
        loguru.logger.info("Compute RFU counts sample by sample...")
        pbar = tqdm.tqdm(query_samples)
        rfu_counts_of_sample = []
        ref_tcrs_columns_to_use = ref_tcrs.select(
            ["sample_name", "v_gene", "cdr3", "centroid"]
        )
        for sample_name in pbar:
            pbar.set_description(sample_name)
            pbar.refresh()
            sample_query_tcrs = query_tcrs.filter(
                pl.col("sample_name") == sample_name
            ).select(["sample_name", "v_gene", "cdr3", "centroid"])
            sample_rfu_counts = filtered_query_rfu_count(
                sample_query_tcrs,
                ref_tcrs_columns_to_use,
                min_size,
                len(ref_samples),
            ).with_columns(pl.lit(sample_name).alias("sample_name"))
            rfu_counts_of_sample.append(sample_rfu_counts)

        # Concatenate RFU counts of all samples and output it as a matrix.
        loguru.logger.info(
            "Converting RFU counts into a matrix and outputting..."
        )
        rfu_counts = (
            pl.concat(rfu_counts_of_sample)
            .pivot(
                values="tcr_count",
                columns="sample_name",
                index="centroid",
            )
            .with_columns(pl.exclude("centroid").fill_null(0.0))
        )
    else:
        loguru.logger.info("Compute RFU counts...")
        rfu_counts = (
            query_tcrs.group_by(["centroid", "sample_name"])
            .len(name="tcr_count")
            .pivot(values="tcr_count", index="centroid", columns="sample_name")
            .with_columns(pl.exclude("centroid").fill_null(0.0))
        )

    with gzip.open(rfu_counts_tsv, "wb") as out_fh:
        rfu_counts.write_csv(out_fh, separator="\t", include_header=True)


if __name__ == "__main__":
    main(
        query_tcrs_ftr=snakemake.input.query_ftr,
        ref_tcrs_ftr=snakemake.input.ref_ftr,
        min_size=snakemake.params.min_size,
        is_cv_test=snakemake.params.is_cv_test,
        rfu_counts_tsv=snakemake.output[0],
    )

"""Annotate and compute statistics of clusters."""

import gzip

import polars as pl


def _basic_cluster_stats(df, predictor):
    """Compute basic statistics for each cluster.

    Argument `df` is a data frame of TCRs.
    """
    df = df.with_columns(
        pl.col(predictor).fill_null(False).alias("positive_response"),
        pl.col(predictor).fill_null(True).not_().alias("negative_response"),
    )
    positive_tcr_count = (
        df.filter(pl.col("positive_response"))
        .group_by("centroid")
        .agg(pl.len().alias("positive_tcr_count"))
    )
    negative_tcr_count = (
        df.filter(pl.col("negative_response"))
        .group_by("centroid")
        .agg(pl.len().alias("negative_tcr_count"))
    )
    positive_subject_count = (
        df.filter(pl.col("positive_response"))
        .unique(["centroid", "sample_name"])
        .group_by("centroid")
        .agg(pl.len().alias("positive_subject_count"))
    )
    negative_subject_count = (
        df.filter(pl.col("negative_response"))
        .unique(["centroid", "sample_name"])
        .group_by("centroid")
        .agg(pl.len().alias("negative_subject_count"))
    )
    df_centroid = df.filter(pl.col("centroid") == pl.col("idx")).unique(
        "centroid"
    )
    centroid_delta_density = df_centroid.select(
        ["centroid", "delta", "density"]
    ).rename({"delta": "centroid_delta", "density": "centroid_density"})
    centroid_tcr = df_centroid.select(["centroid", "v_gene", "cdr3"]).rename(
        {"v_gene": "centroid_v_gene", "cdr3": "centroid_cdr3"}
    )

    cluster_stats = (
        df.group_by("centroid")
        .agg(pl.len().alias("tcr_count"))
        .join(positive_tcr_count, on="centroid", how="left")
        .join(negative_tcr_count, on="centroid", how="left")
        .join(positive_subject_count, on="centroid", how="left")
        .join(negative_subject_count, on="centroid", how="left")
        .join(centroid_delta_density, on="centroid", how="left")
        .join(centroid_tcr, on="centroid", how="left")
        .with_columns(
            pl.col("positive_tcr_count").fill_null(0),
            pl.col("negative_tcr_count").fill_null(0),
            pl.col("positive_subject_count").fill_null(0),
            pl.col("negative_subject_count").fill_null(0),
        )
    )

    return cluster_stats


def main(
    annotated_tcrs,
    metadata_tsv,
    predictor,
    output_tsv,
):
    # Load all inputs.
    tcrs = pl.scan_ipc(annotated_tcrs)
    metadata = pl.scan_csv(
        metadata_tsv,
        has_header=True,
        separator="\t",
        schema_overrides={predictor: bool},
        null_values={predictor: "NA"},
    ).select(["sample_name", predictor])
    tcrs = tcrs.cast({"sample_name": pl.String}).join(
        metadata.cast({predictor: pl.Boolean}), on="sample_name", how="left"
    )

    cluster_stats_df = _basic_cluster_stats(tcrs, predictor).collect()

    with gzip.open(output_tsv, "wb") as output_tsv_gz:
        cluster_stats_df.write_csv(output_tsv_gz, separator="\t")


if __name__ == "__main__":
    main(
        annotated_tcrs=snakemake.input.annotated_tcrs,
        metadata_tsv=snakemake.input.metadata_tsv,
        predictor=snakemake.wildcards.predictor,
        output_tsv=snakemake.output[0],
    )

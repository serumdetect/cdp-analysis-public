"""
Snakemake script for computing a matrix of RFU (rows) counts by sample
(columns). Only RFUs with at least two TCRs are included.
"""

import gzip

import polars as pl
from loguru import logger


def main(
    tcrs_ftr, output_tsv, count_templates, min_subject_count, min_tcr_size=None
):
    logger.info("Reading TCR data...")
    columns = ["sample_name", "centroid"]
    if count_templates:
        columns.append("templates")
    if min_tcr_size is not None:
        columns.append("size")
    tcrs_df = pl.read_ipc(tcrs_ftr, columns=columns)

    if min_tcr_size is not None:
        logger.info("Filtering by TCR size...")
        tcrs_df = tcrs_df.filter(pl.col("size") >= min_tcr_size)

    logger.info("Computing RFU count matrix...")
    # Only keep RFUs with observations from more than one sample.
    rfu_counts_df = tcrs_df.join(
        (
            tcrs_df.group_by("centroid").agg(
                pl.col("sample_name").n_unique().alias("unique_samples")
            )
        ),
        on="centroid",
    ).filter(
        (pl.col("unique_samples") >= min_subject_count)
        & (pl.col("centroid") != -1)
    )
    if count_templates:
        rfu_counts_df = rfu_counts_df.group_by(["sample_name", "centroid"]).agg(
            pl.sum("templates").alias("n")
        )
    else:
        rfu_counts_df = (
            rfu_counts_df.group_by(["sample_name", "centroid"])
            .len()
            .rename({"len": "n"})
        )
    rfu_counts_mat = (
        rfu_counts_df.pivot(
            index="centroid",
            on="sample_name",
            values="n",
        )
        .fill_null(0)
        .with_columns(pl.exclude("centroid").cast(pl.Int32))
    )

    logger.info("Exporting the counts matrix...")
    with gzip.open(output_tsv, "wb") as output_tsv_gz:
        rfu_counts_mat.write_csv(output_tsv_gz, separator="\t")


if __name__ == "__main__":
    main(
        snakemake.input[0],
        snakemake.output[0],
        snakemake.params.count_templates,
        snakemake.params.min_subject_count,
        snakemake.params.get("min_tcr_size"),
    )

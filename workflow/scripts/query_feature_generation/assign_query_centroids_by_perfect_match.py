"""Snakemake script for assigning centroids by perfect amino acid match."""

import polars as pl


def main(query_ftr, ref_ftr, output_ftr):
    # Read the query and reference TCRs.
    query_df = pl.scan_ipc(query_ftr)
    ref_df = (
        pl.scan_ipc(ref_ftr)
        .unique(["v_gene", "cdr3"])
        .select(["v_gene", "cdr3", "centroid"])
    )

    # Join the `centroid` column back to the query dataframe.
    query_df = query_df.join(
        ref_df, on=["v_gene", "cdr3"], how="left"
    ).with_columns(pl.col("centroid").fill_null(-1))

    # Output results.
    query_df.collect().write_ipc(output_ftr)


if __name__ == "__main__":
    main(
        snakemake.input.query_ftr,
        snakemake.input.ref_clusters_ftr,
        snakemake.output[0],
    )

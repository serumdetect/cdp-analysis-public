"""Snakemake script for annotating TCRs with phenotype, CDP and VDJ DB data."""

import polars as pl


def main(tcr_data_ftr, cluster_data_ftr, output_ftr):
    # Load all the data.
    clustering_data = pl.scan_ipc(cluster_data_ftr)
    tcr_data = pl.scan_ipc(tcr_data_ftr)

    # Combine TCRs with clustering information.
    output_df = tcr_data.join(
        clustering_data, how="left", on=["v_gene", "cdr3"]
    )

    output_df.collect().write_ipc(output_ftr, compression="lz4")


if __name__ == "__main__":
    main(
        tcr_data_ftr=snakemake.input.tcr_data_ftr,
        cluster_data_ftr=snakemake.input.cluster_data_ftr,
        output_ftr=snakemake.output[0],
    )

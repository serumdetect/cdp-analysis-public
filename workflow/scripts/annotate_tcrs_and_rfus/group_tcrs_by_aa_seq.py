"""
Snakemake script for assigning a `centroid` column based on columns `v_gene`
and `cdr3`.
"""

import polars as pl


def main(tcrs_ftr, output_ftr):
    # Load input.
    tcr_df = pl.scan_ipc(tcrs_ftr)

    # Group by `v_gene` and `cdr3` and assign a `centroid` column.
    rfu_assignments = (
        tcr_df.group_by(["v_gene", "cdr3"])
        .len("size")
        .with_row_index(name="centroid")
    )

    # Join the `centroid` column back to the original dataframe.
    tcr_df = tcr_df.join(rfu_assignments, on=["v_gene", "cdr3"], how="left")

    # Output results.
    tcr_df.collect().write_ipc(output_ftr)


if __name__ == "__main__":
    main(snakemake.input[0], snakemake.output[0])

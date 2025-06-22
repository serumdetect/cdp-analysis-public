"""Snakemake script for normalizing RFU counts using size factors."""

import pandas as pd


def main(rfu_counts_tsv, size_factors_tsv, output_norm_rfu_counts_tsv):
    rfu_counts = pd.read_csv(rfu_counts_tsv, sep="\t", index_col=0, header=0)
    size_factors = pd.read_csv(
        size_factors_tsv, sep="\t", index_col=0, header=None
    ).squeeze()
    assert set(rfu_counts.columns) == set(size_factors.index)
    norm_rfu_counts = rfu_counts.div(size_factors, axis="columns")
    norm_rfu_counts.to_csv(output_norm_rfu_counts_tsv, sep="\t")


if __name__ == "__main__":
    main(
        snakemake.input.rfu_counts_tsv,
        snakemake.input.size_factors_tsv,
        snakemake.output[0],
    )

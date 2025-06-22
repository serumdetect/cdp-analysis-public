"""
Snakemake script for computing size factors naively as total template counts.
"""

import pandas as pd


def main(rfu_counts_tsv, output_tsv):
    rfu_counts = pd.read_table(rfu_counts_tsv, index_col=0)
    total_rfu_counts = rfu_counts.sum()
    size_factors = total_rfu_counts / total_rfu_counts.mean()
    size_factors.to_csv(output_tsv, sep="\t", header=False, index=True)


if __name__ == "__main__":
    main(
        rfu_counts_tsv=snakemake.input[0],
        output_tsv=snakemake.output[0],
    )

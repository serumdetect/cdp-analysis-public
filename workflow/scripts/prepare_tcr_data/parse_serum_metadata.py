"""Snakemake script for parsing the Serum metadata table."""

import pandas as pd


def main(input_tsv, output_tsv):
    metadata = pd.read_csv(input_tsv, header=0, sep="\t")

    # Add column for total UMI count in millions.
    metadata["umi_count_pre_M"] = metadata["umi_count_pre"] / 1e6

    metadata.replace({True: "TRUE", False: "FALSE"}).to_csv(
        output_tsv, header=True, index=False, sep="\t", na_rep="NA"
    )


if __name__ == "__main__":
    main(snakemake.input[0], snakemake.output[0])

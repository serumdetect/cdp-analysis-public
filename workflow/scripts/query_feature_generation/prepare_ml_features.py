"""Snakemake script for combining RFU features for ML."""

import pandas as pd


def read_features_of_dc(fname, dc):
    """Read a feature matrix and add the d_c as a prefix."""
    df = pd.read_csv(fname, sep="\t", index_col=0)
    df.index = [f"{dc}/{centroid}" for centroid in df.index]
    return df


def main(input_fnames, d_c_values, output_tsv):
    """Combine RFU features for ML."""
    df = pd.concat(
        [
            read_features_of_dc(fname, d_c)
            for fname, d_c in zip(input_fnames, d_c_values)
        ]
    ).T
    df.to_csv(output_tsv, sep="\t", index=True)


if __name__ == "__main__":
    main(
        input_fnames=snakemake.input,
        d_c_values=snakemake.params.d_c_values,
        output_tsv=snakemake.output[0],
    )

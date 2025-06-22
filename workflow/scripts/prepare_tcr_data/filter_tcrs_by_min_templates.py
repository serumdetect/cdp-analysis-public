"""Snakemake script for filtering TCRs by minimum number of templates."""


import pandas as pd


def main(input_ftr, min_templates, output_ftr):
    tcrs = pd.read_feather(input_ftr)
    tcrs = tcrs.query("templates >= @min_templates")
    tcrs.reset_index(drop=True).to_feather(output_ftr)


if __name__ == "__main__":
    main(
        snakemake.input[0],
        snakemake.params["min_templates"],
        snakemake.output[0],
    )

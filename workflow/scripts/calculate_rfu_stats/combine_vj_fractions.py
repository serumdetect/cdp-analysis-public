import pandas as pd


def main(
    input_filepaths,
    output_filepath,
    rfu_prefixes,
):
    rfu_stats = pd.concat(
        [
            (
                pd.read_csv(filepath, sep="\t")
                .pipe(
                    lambda df: df.assign(
                        rfu=f"{rfu_prefix}/" + df["centroid"].astype(str)
                    )
                )
                .drop(columns="centroid")
            )
            for filepath, rfu_prefix in zip(input_filepaths, rfu_prefixes)
        ]
    )[["rfu", "v_gene", "v_gene_fraction", "j_gene", "j_gene_fraction"]]
    rfu_stats.to_csv(output_filepath, sep="\t", index=False, compression="gzip")


if __name__ == "__main__":
    main(
        input_filepaths=snakemake.input,
        output_filepath=snakemake.output[0],
        rfu_prefixes=snakemake.params.rfu_prefixes,
    )

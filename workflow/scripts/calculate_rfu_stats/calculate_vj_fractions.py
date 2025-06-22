import pandas as pd


def main(
    output_filepath,
    tcrs_filepath,
):
    tcrs = pd.read_feather(
        tcrs_filepath, columns=["v_gene", "j_gene", "centroid"]
    ).pipe(lambda df: df[df["centroid"] != -1])
    centroid_totals = tcrs.value_counts("centroid").rename("centroid_total")
    chunks = []
    for vj in ["v", "j"]:
        chunks.append(
            tcrs.value_counts(["centroid", f"{vj}_gene"])
            .rename(f"{vj}_gene_count")
            .reset_index(f"{vj}_gene")
            .join(
                centroid_totals,
            )
            .pipe(
                lambda df: df.assign(
                    **{
                        f"{vj}_gene_fraction": df[f"{vj}_gene_count"]
                        / df["centroid_total"]
                    }
                )
            )[[f"{vj}_gene", f"{vj}_gene_fraction"]]
            .groupby("centroid")
            .apply(
                lambda df: df.iloc[df[f"{vj}_gene_fraction"].argmax()],
                include_groups=False,
            )
        )
    vj_fractions = pd.concat(chunks, axis=1)
    vj_fractions.to_csv(
        output_filepath, sep="\t", index=True, compression="gzip"
    )


if __name__ == "__main__":
    main(
        output_filepath=snakemake.output[0],
        tcrs_filepath=snakemake.input[0],
    )

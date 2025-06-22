"""Snakemake script for combining all cluster information including P values."""

import pandas as pd


def main(cluster_stats_tsv, glm_stats_tsv, output_tsv):
    # Load all inputs.
    cluster_stats_df = pd.read_csv(
        cluster_stats_tsv, header=0, sep="\t", index_col="centroid"
    )
    glm_pvals = pd.read_csv(
        glm_stats_tsv,
        header=0,
        sep="\t",
        index_col=0,
    )
    glm_pvals.columns = "glm_" + glm_pvals.columns

    output_df = cluster_stats_df.merge(
        glm_pvals, how="left", left_index=True, right_index=True
    )
    output_df.index.name = "centroid"
    output_df = output_df.reset_index()
    output_df["centroid"] = output_df["centroid"].astype(int)
    output_df.to_csv(output_tsv, header=True, sep="\t", index=False)


if __name__ == "__main__":
    main(
        snakemake.input.cluster_stats,
        snakemake.input.glm_stats,
        snakemake.output[0],
    )

"""
Snakemake script for joining filtered RFUs with RFU stats and combining then
across CDP parameters.
"""

import pandas as pd


def main(rfu_stats_tsvs, cluster_stats_tsvs, d_c_params, output_tsv):
    inputs = zip(d_c_params, rfu_stats_tsvs, cluster_stats_tsvs)
    cluster_stats_list = []
    for d_c_param, rfu_stats_tsv, cluster_stats_tsv in inputs:
        rfu_stats = pd.read_csv(rfu_stats_tsv, sep="\t")
        cluster_stats = pd.read_csv(cluster_stats_tsv, sep="\t")
        cluster_stats_joined = rfu_stats.merge(
            cluster_stats, how="left", on="centroid"
        )
        cluster_stats_joined.insert(0, "d_c", d_c_param)
        cluster_stats_list.append(cluster_stats_joined)
    cluster_stats_all = pd.concat(
        cluster_stats_list, ignore_index=True
    )
    cluster_stats_all["rfu"] = (
        cluster_stats_all["d_c"]
        + "/"
        + cluster_stats_all["centroid"].astype(str)
    )

    # Output results.
    cluster_stats_all.to_csv(output_tsv, index=False, header=True, sep="\t")


if __name__ == "__main__":
    main(
        rfu_stats_tsvs=snakemake.input.rfu_stats,
        cluster_stats_tsvs=snakemake.input.cluster_stats,
        d_c_params=snakemake.params.d_c_params,
        output_tsv=snakemake.output[0],
    )

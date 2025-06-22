"""Snakemake script for matching query TCRs with reference TCRs RFUs.

Matching is done by searching the TCR index and finding the nearest reference
TCR that is within the d_c cutoff. The neighboring TCRs are prioritized by
density.
"""

import numpy as np
import pandas as pd


def main(
    query_ftr,
    query_nearest_ref_tcrs_ftr,
    ref_clusters_ftr,
    d_c_cutoff,
    n_jobs,
    output_ftr,
):
    # Load the data.
    query_tcrs = pd.read_feather(query_ftr)
    nearest_ref_tcrs = pd.read_feather(query_nearest_ref_tcrs_ftr)

    # Compute query RFU centroids (when disregarding delta cutoff).
    ref_clusters = (
        pd.read_feather(
            ref_clusters_ftr, columns=["v_gene", "cdr3", "centroid"]
        )
        .rename(
            columns={"v_gene": "ref_v_gene", "cdr3": "ref_cdr3"},
        )
        .drop_duplicates()
    )
    nearest_ref_tcrs = nearest_ref_tcrs.merge(
        ref_clusters, how="left", on=["ref_v_gene", "ref_cdr3"]
    )

    query_tcr_centroids = query_tcrs.merge(
        nearest_ref_tcrs,
        how="left",
        on=["v_gene", "cdr3"],
    )

    # Only keep counts with a d_c less than the cutoff. Mark the remaining as
    # no centroid (-1).
    query_tcr_centroids.loc[
        (query_tcr_centroids["delta"] > d_c_cutoff)
        & (~np.isclose(query_tcr_centroids["delta"], d_c_cutoff)),
        "centroid",
    ] = -1

    # Output data.
    query_tcr_centroids.to_feather(output_ftr)


if __name__ == "__main__":
    main(
        query_ftr=snakemake.input.query_ftr,
        query_nearest_ref_tcrs_ftr=snakemake.input.query_nearest_ref_tcrs_ftr,
        ref_clusters_ftr=snakemake.input.ref_clusters_ftr,
        d_c_cutoff=snakemake.params.d_c_cutoff,
        n_jobs=snakemake.threads,
        output_ftr=snakemake.output.tcr_centroids,
    )

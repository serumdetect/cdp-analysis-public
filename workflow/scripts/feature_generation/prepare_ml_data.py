"""Snakemake script for preparing cross-validation RFU data."""

import os

import loguru
import pandas as pd
import tqdm


def read_inputs(
    counts_mat_tsv,
    glm_stats_tsv,
    cdp_settings_list,
    fdr_cutoff,
    positive_association_only,
    predictor,
):
    ml_data = dict()
    inputs = zip(cdp_settings_list, counts_mat_tsv, glm_stats_tsv)
    pbar = tqdm.tqdm(inputs)
    for cdp_setting, fname, glm_stats_fname in pbar:
        pbar.set_description(cdp_setting)

        if os.path.getsize(fname) <= 21:  # Size of empty gzip file
            continue

        rfu_counts = pd.read_csv(fname, header=0, index_col=0, sep="\t")

        # Only keep RFUs with cancer FDR <= fdr_cutoff and optionally only
        # RFUs positively associated with cancer.
        glm_stats = pd.read_csv(
            glm_stats_fname, header=0, index_col="centroid", sep="\t"
        )
        glm_stats = glm_stats.loc[glm_stats[f"{predictor}.fdr"] <= fdr_cutoff]
        if positive_association_only:
            # Filter out RFUs negatively associated with cancer.
            glm_stats = glm_stats.loc[
                glm_stats[f"{predictor}.log2FoldChange"] > 0
            ]

        centroids_to_keep = glm_stats.index
        is_centroid_to_keep = rfu_counts.index.isin(centroids_to_keep)
        if sum(is_centroid_to_keep) < len(centroids_to_keep):
            loguru.logger.warning(
                (
                    "Only {} out of {} centroids are present in the RFU counts"
                    " matrix."
                ).format(sum(is_centroid_to_keep), len(centroids_to_keep))
            )
        rfu_counts = rfu_counts.reindex(centroids_to_keep).fillna(0)

        # Combine current d_c settings with the RFU as the index.
        rfu_counts = rfu_counts.set_index(
            cdp_setting + "/" + rfu_counts.index.to_series().astype(str),
        )
        ml_data[cdp_setting] = rfu_counts

    concatenated_ml_data = pd.concat(ml_data.values()).T
    return concatenated_ml_data


def main(
    counts_mat_tsv,
    glm_stats_tsv,
    cdp_settings_list,
    fdr_cutoff,
    positive_association_only,
    predictor,
    output_tsv,
):
    ml_data = read_inputs(
        counts_mat_tsv,
        glm_stats_tsv,
        cdp_settings_list,
        fdr_cutoff,
        positive_association_only,
        predictor,
    )
    ml_data.to_csv(output_tsv, header=True, index=True, sep="\t")


if __name__ == "__main__":
    main(
        snakemake.input.rfu_counts,
        snakemake.input.glm_stats_tsv,
        snakemake.params.params,
        snakemake.params.fdr_cutoff,
        snakemake.params.positive_association_only,
        snakemake.wildcards.predictor,
        snakemake.output[0],
    )

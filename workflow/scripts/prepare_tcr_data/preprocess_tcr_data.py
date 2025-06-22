"""Snakemake script for preparing Serum TCR data."""

import functools

import loguru
import pandas as pd
import scipy.stats


def load_tcr_data(data_path, metadata_path):
    loguru.logger.info("Loading train data.")
    data = pd.read_feather(data_path)
    metadata = pd.read_csv(metadata_path, sep="\t")
    return data, metadata


def compute_tcr_stats(data, metadata):
    """Compute various statistics about the TCR data."""
    grouped_templates = data.groupby("sample_name", observed=True)["templates"]
    tcr_stats = pd.DataFrame(
        {
            "clonotype_count_pre": grouped_templates.size(),
            "umi_count_pre": grouped_templates.sum(),
            "median_umi_count": grouped_templates.median(),
            "mean_umi_count_iqr": grouped_templates.apply(
                functools.partial(scipy.stats.trim_mean, proportiontocut=0.25)
            ),
        }
    )
    metadata = metadata.drop(
        columns=[col for col in tcr_stats.columns if col != "sample_name"],
        errors="ignore",
    ).merge(tcr_stats, "left", on="sample_name")
    return metadata


def filter_by_blacklisted_tcrs(data, blacklisted_tcrs_table):
    """Remove TCRs matching blacklisted TCRs."""
    join_cols = ["v_gene", "j_gene", "cdr3_nt"]
    merged_data = pd.merge(
        data,
        blacklisted_tcrs_table[join_cols],
        on=join_cols,
        how="left",
        indicator=True,
        validate="many_to_one",
    )
    data = merged_data.loc[merged_data["_merge"] == "left_only"].drop(
        columns="_merge"
    )
    filtered_data = merged_data.loc[merged_data["_merge"] != "left_only"]
    n_tcrs_filtered = filtered_data.shape[0]
    n_samples_filtered = filtered_data["sample_name"].nunique()
    loguru.logger.warning(
        "Filtered out {}/{} blacklisted TCRs from {} samples".format(
            n_tcrs_filtered, len(data), n_samples_filtered
        )
    )
    return data


def filter_data(
    data,
    metadata,
    min_productive_rearrangements,
    min_productive_templates,
    cdr3_lens,
    blacklisted_tcrs_csv,
    is_external_data=False,
):
    """Filter the TCR data based on various criteria."""
    # Filter out TCRs by CDR3 length.
    data = data.loc[data["cdr3"].str.len().isin(cdr3_lens)]

    # Filter out samples by min_productive_rearrangements and
    # min_productive_templates.
    metadata = metadata.loc[
        (metadata["productive_rearrangements"] >= min_productive_rearrangements)
        & (metadata["productive_templates"] >= min_productive_templates)
    ]
    data = data.loc[data["sample_name"].isin(metadata["sample_name"].tolist())]

    # External test data is only filtered for blacklisted TCRs.
    if blacklisted_tcrs_csv is not None:
        blacklisted_tcrs = pd.read_csv(blacklisted_tcrs_csv)
        data = filter_by_blacklisted_tcrs(data, blacklisted_tcrs)
    if data["sample_name"].dtype.name == "category":
        data["sample_name"] = data["sample_name"].cat.remove_unused_categories()
    return data, metadata


def sanity_check_data(data, metadata):
    # Make sure the TCR data is consistent with the metadata.
    valid_sample_names = set(metadata["sample_name"])
    observed_sample_names = set(data["sample_name"])
    if not observed_sample_names.issubset(valid_sample_names):
        raise ValueError("Found sample names in data that are not in metadata.")


def main(
    blacklisted_tcrs_csv,
    min_productive_rearrangements,
    min_productive_templates,
    cdr3_lens,
    data_path,
    metadata_path,
    output_ftr,
    output_metadata_tsv,
):
    # Load the TCR data and metadata.
    data, metadata = load_tcr_data(
        data_path=data_path, metadata_path=metadata_path
    )

    # Assume that all V and J genes are from allele 01.
    data["v_gene"] = data["v_gene"].astype("string") + "*01"
    data["j_gene"] = data["j_gene"].astype("string") + "*01"

    # Compute various statistics about the TCR data.
    metadata = compute_tcr_stats(data, metadata)

    # Filter the TCR data based on various criteria.
    data, metadata = filter_data(
        data=data,
        metadata=metadata,
        min_productive_rearrangements=min_productive_rearrangements,
        min_productive_templates=min_productive_templates,
        cdr3_lens=cdr3_lens,
        blacklisted_tcrs_csv=blacklisted_tcrs_csv,
    )

    # Sanity check the data.
    sanity_check_data(data=data, metadata=metadata)

    data.to_feather(output_ftr)
    metadata.to_csv(output_metadata_tsv, header=True, sep="\t", index=False)


if __name__ == "__main__":
    main(
        data_path=snakemake.input.data_path,
        metadata_path=snakemake.input.metadata_path,
        # An empty list is equivalent to None.
        blacklisted_tcrs_csv=snakemake.input.blacklisted_tcrs or None,
        min_productive_rearrangements=snakemake.params[
            "min_productive_rearrangements"
        ],
        min_productive_templates=snakemake.params["min_productive_templates"],
        cdr3_lens=snakemake.params.cdr3_lens,
        output_ftr=snakemake.output.data,
        output_metadata_tsv=snakemake.output.metadata,
    )

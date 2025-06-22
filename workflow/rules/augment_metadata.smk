"""Rules to augment metadata with additional columns.

Some additional columns are calculated on unfiltered data and appended to the
parsed metadata file of the unfiltered data.

Columns determined from filtered data are appended to the corresponding
unfiltered augmented metadata file.
"""


rule calculate_unfiltered_sequencing_error_count:
    """Calculates number of clonotypes with UMI count at most 2.

    This count is a rough proxy to sequencing error counts.
    """
    output:
        "{all_data_splits}/tcr_stats/mRFU/unfiltered_sequencing_error_count.tsv",
    input:
        "{all_data_splits}/raw/mRFU.ftr",
    log:
        "{all_data_splits}/tcr_stats/mRFU/unfiltered_sequencing_error_count.tsv.log",
    conda:
        "../envs/augment_metadata/augment_metadata.yaml"
    resources:
        mem_mb=lambda _, input, attempt: size_mb(input[0])
        * 5
        * (0.5 + attempt / 2),
    script:
        "../scripts/augment_metadata/calculate_unfiltered_sequencing_error_count.py"


rule augment_unfiltered_metadata:
    """Augment parsed metadata of unfiltered data with additional columns.

    The additional columns should contain information that is calculated from
    the unfiltered data.
    """
    output:
        "{any_unfiltered_prefix}/augmented_metadata/mRFU.metadata.tsv",
    input:
        metadata="{any_unfiltered_prefix}/parsed_metadata/mRFU.metadata.tsv",
        sequencing_error_counts="data/tcr_stats/mRFU/unfiltered_sequencing_error_count.tsv",
    log:
        "{any_unfiltered_prefix}/augmented_metadata/mRFU.metadata.tsv.log",
    conda:
        "../envs/augment_metadata/augment_metadata.yaml"
    script:
        "../scripts/augment_metadata/augment_unfiltered_metadata.py"


rule calculate_filtered_counts:
    output:
        "{all_data_splits}/{tcr_filtering}/tcr_stats/mRFU/filtered_counts.tsv",
    input:
        "{all_data_splits}/{tcr_filtering}/raw/mRFU.ftr",
    log:
        "{all_data_splits}/{tcr_filtering}/tcr_stats/mRFU/filtered_counts.log",
    conda:
        "../envs/augment_metadata/augment_metadata.yaml"
    resources:
        mem_mb=lambda _, input, attempt: size_mb(input[0])
        * 5
        * (0.5 + attempt / 2),
    script:
        "../scripts/augment_metadata/calculate_filtered_counts.py"


def unfiltered_prefix_from_filtered_prefix(filtered_prefix):
    """Removes the filtering infix from the prefix."""
    fields = filtered_prefix.split("/")
    return "/".join([fields[0], *fields[2:]])


rule augment_filtered_metadata:
    """Augment augmented metadata of unfiltered data with additional columns.

    The additional columns should contain information that is calculated from
    the filtered data.
    """
    output:
        "{any_filtered_prefix}/augmented_metadata/mRFU.metadata.tsv",
    input:
        filtered_counts=lambda wildcards: (
            "data/"
            + "/".join(
                re.sub(
                    "/rs~\d+/cv~\d+", "", wildcards.any_filtered_prefix
                    ).split("/")[1:]
                )
            + "/tcr_stats/mRFU/filtered_counts.tsv"
        ),
        metadata=lambda wildcards: (
            f"{unfiltered_prefix_from_filtered_prefix(wildcards.any_filtered_prefix)}/augmented_metadata/mRFU.metadata.tsv"
        ),
    log:
        "{any_filtered_prefix}/augmented_metadata/mRFU.metadata.tsv.log",
    conda:
        "../envs/augment_metadata/augment_metadata.yaml"
    script:
        "../scripts/augment_metadata/augment_filtered_metadata.py"

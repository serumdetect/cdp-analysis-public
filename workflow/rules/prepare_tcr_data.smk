"""Snakemake rules for preparing raw TCR data and indexing the TCRs."""

import snakemake.remote.S3

S3 = snakemake.remote.S3.RemoteProvider()


def local_or_s3_path(path):
    """Determine whether the current file needs to be loaded from S3."""
    if path.startswith("s3://"):
        path = S3.remote(path)
    return path


rule preprocess_tcr_data:
    """Preprocess and validate "global" TCR data using."""
    output:
        data="data/raw/mRFU.ftr",
        metadata="data/metadata/mRFU.metadata.tsv",
    input:
        data_path=local_or_s3_path(config["raw_data_train"]),
        metadata_path=local_or_s3_path(config["metadata_train"]),
        # Snakemake requires a list or a string as an input path.
        blacklisted_tcrs=lambda _: (
            []
            if config["blacklisted_tcrs"] is None
            else local_or_s3_path(config["blacklisted_tcrs"])
        ),
    benchmark:
        "data/raw/mRFU.ftr.benchmark"
    params:
        min_productive_rearrangements=config["tcr_params"][
            "min_productive_rearrangements"
        ],
        min_productive_templates=config["tcr_params"][
            "min_productive_templates"
        ],
        cdr3_lens=config["tcr_params"]["cdr3_len"],
    resources:
        mem_mb=lambda wildcards: config["mem_mb"],
    conda:
        "../envs/prepare_tcr_data/preprocess_tcr_data.yaml"
    script:
        "../scripts/prepare_tcr_data/preprocess_tcr_data.py"


# Rule for loading external test data against global/CV fold train data.
use rule preprocess_tcr_data as preprocess_tcr_test_data with:
    output:
        data="{external_data}/raw/mRFU.ftr",
        metadata="{external_data}/metadata/mRFU.metadata.tsv",
    input:
        data_path=lambda wildcards: local_or_s3_path(
            config["test_data"][wildcards.external_data]["data"]
        ),
        metadata_path=lambda wildcards: local_or_s3_path(
            config["test_data"][wildcards.external_data]["metadata"]
        ),
        # Snakemake requires a list or a string as an input path.
        blacklisted_tcrs=lambda _: (
            []
            if config["blacklisted_tcrs"] is None
            else local_or_s3_path(config["blacklisted_tcrs"])
        ),
    benchmark:
        "{external_data}/raw/mRFU.ftr.benchmark"


rule filter_tcrs_by_min_templates:
    """Filter TCR data based on minimum template count."""
    output:
        "{data_or_test_data_splits}/min_templates~{min_templates}/raw/mRFU.ftr",
    input:
        "{data_or_test_data_splits}/raw/mRFU.ftr",
    benchmark:
        "{data_or_test_data_splits}/min_templates~{min_templates}/raw/mRFU.ftr.benchmark"
    params:
        min_templates=lambda wildcards: int(wildcards.min_templates),
    # According to https://github.com/serumdetect/snakemake-resource-usage/blob/f77fdcdf983e21bce2c97ede4dadeb4a34718e19/notebooks/analyze-input-sizes-vs-memory-use.2024-05-28.py.ipynb,
    # this rule uses roughly 10 * input file size.
    resources:
        mem_mb=lambda _, input, attempt: size_mb(input[0])
        * 10
        * (0.5 + attempt / 2),
    conda:
        "../envs/prepare_tcr_data/filter_tcrs.yaml"
    script:
        "../scripts/prepare_tcr_data/filter_tcrs_by_min_templates.py"


rule filter_tcrs_by_min_fold_of_templates_quantile:
    """
    Keep clonotypes with at least a certain fold of the median template count.
    """
    output:
        "{data_or_test_data_splits}/min_fold_of_{v_or_nothing}quantile~{min_fold_of_stat}/raw/mRFU.ftr",
    input:
        "{data_or_test_data_splits}/raw/mRFU.ftr",
    log:
        stderr="{data_or_test_data_splits}/min_fold_of_{v_or_nothing}quantile~{min_fold_of_stat}/raw/mRFU.ftr.stderr",
    benchmark:
        "{data_or_test_data_splits}/min_fold_of_{v_or_nothing}quantile~{min_fold_of_stat}/raw/mRFU.ftr.benchmark"
    params:
        min_fold_of_stat=lambda wildcards: float(wildcards.min_fold_of_stat),
        reference_quantile=0.5,
        v_gene_specific=lambda wildcards: wildcards.v_or_nothing == "v_",
    # According to https://github.com/serumdetect/snakemake-resource-usage/blob/339dfbc5eb336a9c8ddab1796a33f685a172e31e/notebooks/analyze-input-sizes-vs-memory-use.2024-05-28.py.ipynb,
    # this rule uses roughly 10 * input file size.
    wildcard_constraints:
        v_or_nothing="v_|",
    resources:
        mem_mb=lambda _, input, attempt: size_mb(input[0])
        * 10
        * (0.5 + attempt / 2),
    conda:
        "../envs/prepare_tcr_data/filter_tcrs.yaml"
    script:
        "../scripts/prepare_tcr_data/filter_tcrs_by_min_fold_of_templates_stat.py"


# Rule for filtering TCRs by trimmed mean instead of quantile.
use rule filter_tcrs_by_min_fold_of_templates_quantile as filter_tcrs_by_min_fold_of_templates_mean with:
    output:
        "{data_or_test_data_splits}/min_fold_of_{v_or_nothing}mean~{min_fold_of_stat}/raw/mRFU.ftr",
    input:
        "{data_or_test_data_splits}/raw/mRFU.ftr",
    log:
        stderr="{data_or_test_data_splits}/min_fold_of_{v_or_nothing}mean~{min_fold_of_stat}/raw/mRFU.ftr.stderr",
    benchmark:
        "{data_or_test_data_splits}/min_fold_of_{v_or_nothing}mean~{min_fold_of_stat}/raw/mRFU.ftr.benchmark"
    params:
        min_fold_of_stat=lambda wildcards: float(wildcards.min_fold_of_stat),
        reference_quantile=None,
        v_gene_specific=lambda wildcards: wildcards.v_or_nothing == "v_",


rule filter_tcrs_by_top_num_clonotypes:
    """Keep top N clonotypes. Remove samples with fewer than N clonotypes."""
    output:
        "{data_or_test_data_splits}/top_num_clonotypes_{filtering_strategy}~{top_n}/raw/mRFU.ftr",
    input:
        "{data_or_test_data_splits}/raw/mRFU.ftr",
    log:
        stderr="{data_or_test_data_splits}/top_num_clonotypes_{filtering_strategy}~{top_n}/raw/mRFU.ftr.stderr",
    benchmark:
        "{data_or_test_data_splits}/top_num_clonotypes_{filtering_strategy}~{top_n}/raw/mRFU.ftr.benchmark"
    wildcard_constraints:
        filtering_strategy="raw_ranks|norm_templates|norm_ranks",
    params:
        top_n=lambda wildcards: int(wildcards.top_n),
        # Quantile using which to normalize the clonotype template counts within
        # each V gene. Only used with filtering_strategy~norm_templates.
        reference_quantile=0.25,
    resources:
        mem_mb=lambda _, input, attempt: size_mb(input[0])
        * 10
        * (0.5 + attempt / 2),
    conda:
        "../envs/prepare_tcr_data/filter_tcrs.yaml"
    script:
        "../scripts/prepare_tcr_data/filter_tcrs_by_top_num_clonotypes.py"


rule copy_tcr_filtering_metadata:
    """
    Just copy the metadata file for a TCR filtering to keep the file paths
    consistent.
    """
    localrule: True
    output:
        "{data_or_test_data_splits}/{tcr_filtering}/metadata/mRFU.metadata.tsv",
    input:
        "{data_or_test_data_splits}/metadata/mRFU.metadata.tsv",
    shell:
        """
        cp {input} {output}
        """


if config["cv"]["num_cv_folds"] is not None:

    rule prepare_train_test_splits:
        """Prepare train-test splits for the TCR data."""
        output:
            data="data/rs~{rs}/cv~{cv}/raw/mRFU.ftr",
            metadata="data/rs~{rs}/cv~{cv}/metadata/mRFU.metadata.tsv",
            held_out_data="held_out_data/rs~{rs}/cv~{cv}/raw/mRFU.ftr",
            held_out_metadata="held_out_data/rs~{rs}/cv~{cv}/metadata/mRFU.metadata.tsv",
        input:
            data="data/raw/mRFU.ftr",
            metadata="data/metadata/mRFU.metadata.tsv",
        benchmark:
            "data/rs~{rs}/cv~{cv}/raw/mRFU.ftr.benchmark"
        params:
            cv_grouping_colname=config["cv"]["cv_grouping_colname"],
            num_cv_folds=config["cv"]["num_cv_folds"],
            rs=lambda wildcards: int(wildcards.rs),
            cv=lambda wildcards: int(wildcards.cv),
        conda:
            "../envs/prepare_tcr_data/prepare_train_test_splits.yaml"
        # According to https://github.com/serumdetect/snakemake-resource-usage/blob/f77fdcdf983e21bce2c97ede4dadeb4a34718e19/notebooks/analyze-input-sizes-vs-memory-use.2024-05-28.py.ipynb,
        # this rule uses ~10 * input file size.
        resources:
            mem_mb=lambda _, input, attempt: 10
            * size_mb(input.data)
            * (0.5 + attempt / 2),
        script:
            "../scripts/prepare_tcr_data/prepare_train_test_splits.py"

    use rule prepare_train_test_splits as prepare_train_test_splits_with_filtering with:
        output:
            data="data/{tcr_filtering}/rs~{rs}/cv~{cv}/raw/mRFU.ftr",
            metadata="data/{tcr_filtering}/rs~{rs}/cv~{cv}/metadata/mRFU.metadata.tsv",
            held_out_data="held_out_data/{tcr_filtering}/rs~{rs}/cv~{cv}/raw/mRFU.ftr",
            held_out_metadata="held_out_data/{tcr_filtering}/rs~{rs}/cv~{cv}/metadata/mRFU.metadata.tsv",
        input:
            data="data/{tcr_filtering}/raw/mRFU.ftr",
            metadata="data/metadata/mRFU.metadata.tsv",
        benchmark:
            "data/{tcr_filtering}/rs~{rs}/cv~{cv}/raw/mRFU.ftr.benchmark"

    rule prepare_external_test_splits:
        """External data test splits are simple copies of the external data."""
        output:
            data="{external_data}/rs~{rs}/cv~{cv}/raw/mRFU.ftr",
            metadata="{external_data}/rs~{rs}/cv~{cv}/metadata/mRFU.metadata.tsv",
        input:
            data="{external_data}/raw/mRFU.ftr",
            metadata="{external_data}/metadata/mRFU.metadata.tsv",
        shell:
            """
            cp {input.data} {output.data}
            cp {input.metadata} {output.metadata}
            """

    use rule prepare_external_test_splits as prepare_external_test_splits_with_filtering with:
        output:
            data="{external_data}/{tcr_filtering}/rs~{rs}/cv~{cv}/raw/mRFU.ftr",
            metadata="{external_data}/{tcr_filtering}/rs~{rs}/cv~{cv}/metadata/mRFU.metadata.tsv",
        input:
            data="{external_data}/{tcr_filtering}/raw/mRFU.ftr",
            metadata="{external_data}/metadata/mRFU.metadata.tsv",


rule parse_serum_metadata:
    """Preprocess clinical covariates for the Serum metadata."""
    localrule: True
    output:
        "{any_prefix}/parsed_metadata/mRFU.metadata.tsv",
    input:
        "{any_prefix}/metadata/mRFU.metadata.tsv",
    conda:
        "../envs/prepare_tcr_data/parse_serum_metadata.yaml"
    script:
        "../scripts/prepare_tcr_data/parse_serum_metadata.py"


rule index_tcr_data:
    """Create and index TCR repertoire."""
    output:
        f"{{prefix_no_drop_covariates}}/nn_indices/mRFU/{dist_normalization_paramspace.wildcard_pattern}.sav",
    input:
        "{prefix_no_drop_covariates}/raw/mRFU.ftr",
    benchmark:
        f"{{prefix_no_drop_covariates}}/nn_indices/mRFU/{dist_normalization_paramspace.wildcard_pattern}.sav.benchmark"
    threads: lambda wildcards: config["threads"]
    # According to https://github.com/serumdetect/snakemake-resource-usage/blob/339dfbc5eb336a9c8ddab1796a33f685a172e31e/notebooks/analyze-input-sizes-vs-memory-use.2024-05-28.py.ipynb,
    # this rule uses around 10GB + 60 * input file size.
    resources:
        mem_mb=lambda _, input, attempt: 10000
        + 80 * size_mb(input[0]) * (0.5 + attempt / 2),
    params:
        n_neighbors=30,
    # Requires internal package tcrnn.
    container:
        config["serum_rfu_formation_docker"]
    script:
        "../scripts/prepare_tcr_data/index_tcr_data.py"


rule compute_density:
    """Compute density of CFSFDP using ANN search"""
    group:
        "cdp"
    output:
        f"{{prefix_no_drop_covariates}}/densities/mRFU/{cdp_paramspace.wildcard_pattern}.ftr",
    input:
        f"{{prefix_no_drop_covariates}}/nn_indices/mRFU/{dist_normalization_paramspace.wildcard_pattern}.sav",
    benchmark:
        f"{{prefix_no_drop_covariates}}/densities/mRFU/{cdp_paramspace.wildcard_pattern}.ftr.benchmark"
    threads: lambda wildcards: config["threads"]
    # According to https://github.com/serumdetect/snakemake-resource-usage/blob/f77fdcdf983e21bce2c97ede4dadeb4a34718e19/notebooks/analyze-input-sizes-vs-memory-use.2024-05-28.py.ipynb,
    # this rule uses roughly 4 * input file size.
    resources:
        mem_mb=lambda _, input, attempt: size_mb(input[0])
        * 5
        * (0.5 + attempt / 2),
    # Requires internal package tcrnn.
    container:
        config["serum_rfu_formation_docker"]
    script:
        "../scripts/prepare_tcr_data/compute_density.py"


rule compute_decision_graph:
    """Compute density and delta of CFSFDP using ANN search"""
    group:
        "cdp"
    output:
        f"{{prefix_no_drop_covariates}}/decision_graph/mRFU/{cdp_paramspace.wildcard_pattern}.ftr",
    input:
        tcr_repertoire=f"{{prefix_no_drop_covariates}}/nn_indices/mRFU/{dist_normalization_paramspace.wildcard_pattern}.sav",
        density=f"{{prefix_no_drop_covariates}}/densities/mRFU/{cdp_paramspace.wildcard_pattern}.ftr",
    benchmark:
        f"{{prefix_no_drop_covariates}}/decision_graph/mRFU/{cdp_paramspace.wildcard_pattern}.ftr.benchmark"
    threads: lambda wildcards: config["threads"]
    # According to https://github.com/serumdetect/snakemake-resource-usage/blob/f77fdcdf983e21bce2c97ede4dadeb4a34718e19/notebooks/analyze-input-sizes-vs-memory-use.2024-05-28.py.ipynb,
    # this rule uses roughly 4 * input file size.
    resources:
        mem_mb=lambda _, input, attempt: size_mb(input.tcr_repertoire)
        * 5
        * (0.5 + attempt / 2),
    # Requires internal package tcrnn.
    container:
        config["serum_rfu_formation_docker"]
    script:
        "../scripts/prepare_tcr_data/compute_decision_graph.py"

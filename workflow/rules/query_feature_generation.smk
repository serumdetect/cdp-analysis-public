"""Snakemake workflow for generating features for the query (i.e. test) data.

These rules are derived from the predict-query-repertoire package.
"""

import re


def replace_root_with_data(path):
    """Performs the following replacement: <root>/... --> data/..."""
    return re.sub(r"^[^/]+", "data", path)


rule fuzzy_match_ref_tcrs:
    """
    Match query TCRs with reference TCRs using inexact matching of the TCR
    index.
    """
    output:
        query_nearest_ref_tcrs_ftr="{test_data_prefix}/nearest_ref_tcrs/mRFU/normalize_dist~{normalize_dist}.ftr",
    input:
        query_ftr="{test_data_prefix}/raw/mRFU.ftr",
        ref_repertoire_path=lambda wildcards: (
            replace_root_with_data(wildcards.test_data_prefix)
            + "/nn_indices/mRFU/normalize_dist~{normalize_dist}.sav"
        ),
    benchmark:
        "{test_data_prefix}/nearest_ref_tcrs/mRFU/normalize_dist~{normalize_dist}.ftr.benchmark"
    threads: lambda wildcards: config["threads"]
    # According to https://github.com/serumdetect/snakemake-resource-usage/blob/f77fdcdf983e21bce2c97ede4dadeb4a34718e19/notebooks/analyze-input-sizes-vs-memory-use.2024-05-28.py.ipynb,
    # this rule always ~2 * input file size.
    resources:
        mem_mb=lambda _, input, attempt: 2
        * size_mb(input.ref_repertoire_path)
        * (0.5 + attempt / 2),
    params:
        normalize_dist=lambda wildcards: int(wildcards.normalize_dist),
        monkey_patch_ref_repertoire=False,
    container:
        config["serum_rfu_formation_docker"]
    script:
        "../scripts/query_feature_generation/fuzzy_match_ref_tcrs.py"


rule assign_query_centroids:
    """Assign query TCRs to reference TCR centroids.

    The reference TCR centroid is taken as the nearest TCR to the query TCR if
    the distance between the two is less than the `d_c` cutoff.
    """
    output:
        tcr_centroids=f"{{test_data_prefix}}/annotated_tcrs/mRFU/{cdp_paramspace.wildcard_pattern}.ftr",
    input:
        query_ftr="{test_data_prefix}/raw/mRFU.ftr",
        query_nearest_ref_tcrs_ftr="{test_data_prefix}/nearest_ref_tcrs/mRFU/normalize_dist~{normalize_dist}.ftr",
        ref_clusters_ftr=lambda wildcards: (
            replace_root_with_data(wildcards.test_data_prefix)
            + f"/decision_graph/mRFU/{cdp_paramspace.wildcard_pattern}.ftr"
        ),
    benchmark:
        f"{{test_data_prefix}}/annotated_tcrs/mRFU/{cdp_paramspace.wildcard_pattern}.ftr.benchmark"
    # According to https://github.com/serumdetect/snakemake-resource-usage/blob/339dfbc5eb336a9c8ddab1796a33f685a172e31e/notebooks/analyze-input-sizes-vs-memory-use.2024-05-28.py.ipynb,
    # this rule always uses around 10x the input file size.
    resources:
        mem_mb=lambda _, input, attempt: 10
        * size_mb(input.ref_clusters_ftr)
        * (0.5 + attempt / 2),
    params:
        d_c_cutoff=lambda wildcards: float(wildcards.d_c),
    conda:
        "../envs/query_feature_generation/assign_query_centroids.yaml"
    script:
        "../scripts/query_feature_generation/assign_query_centroids.py"


rule assign_query_centroids_by_perfect_match:
    """
    Assign query TCRs to reference TCR centroids by perfect amino acid match.
    """
    output:
        tcr_centroids=f"{{test_data_prefix}}/annotated_tcrs/mRFU/normalize_dist~0/d_c~0.ftr",
    input:
        query_ftr="{test_data_prefix}/raw/mRFU.ftr",
        ref_clusters_ftr=lambda wildcards: (
            replace_root_with_data(wildcards.test_data_prefix)
            + f"/annotated_tcrs/mRFU/normalize_dist~0/d_c~0.ftr"
        ),
    benchmark:
        f"{{test_data_prefix}}/annotated_tcrs/mRFU/normalize_dist~0/d_c~0.ftr.benchmark"
    resources:
        mem_mb=lambda _, input, attempt: 10
        * size_mb(input.ref_clusters_ftr)
        * (0.5 + attempt / 2),
    conda:
        "../envs/query_feature_generation/assign_query_centroids.yaml"
    script:
        "../scripts/query_feature_generation/assign_query_centroids_by_perfect_match.py"


ruleorder: assign_query_centroids_by_perfect_match > assign_query_centroids


rule compute_query_rfu_counts:
    """
    Compute RFU counts for query samples while accounting for parameter `size`.

    This rule takes into account the fact that the size estimates are slightly
    biased in query samples since size is computed for num_reF_samples + 1 query
    sample.
    """
    output:
        rfu_counts=f"{{test_data_prefix}}/rfu_counts/mRFU/{cdp_paramspace.wildcard_pattern}/{rfu_size_paramspace.wildcard_pattern}.tsv.gz",
    input:
        query_ftr=f"{{test_data_prefix}}/annotated_tcrs/mRFU/{cdp_paramspace.wildcard_pattern}.ftr",
        ref_ftr=lambda wildcards: (
            replace_root_with_data(wildcards.test_data_prefix)
            + f"/annotated_tcrs/mRFU/{cdp_paramspace.wildcard_pattern}.ftr"
        ),
    benchmark:
        f"{{test_data_prefix}}/rfu_counts/mRFU/{cdp_paramspace.wildcard_pattern}/{rfu_size_paramspace.wildcard_pattern}.tsv.gz.benchmark"
    params:
        min_size=lambda wildcards: int(wildcards.min_size),
        is_cv_test=lambda wildcards: wildcards.test_data_prefix.startswith(
            "held_out_data/"
        ),
    # According to https://github.com/serumdetect/snakemake-resource-usage/blob/782943ec7fd1ba5fc448e1719c8e5a5da1fa577e/notebooks/analyze-input-sizes-vs-memory-use.2024-05-28.py.ipynb,
    # this rule uses below 5x input size.
    resources:
        mem_mb=lambda _, input, attempt: 5
        * size_mb(input.ref_ftr)
        * (0.5 + attempt / 2),
    conda:
        "../envs/query_feature_generation/compute_query_rfu_counts.yaml"
    script:
        "../scripts/query_feature_generation/compute_query_rfu_counts.py"


rule query_size_factors:
    """Compute size factors for query TCRs.

    Size factors are computed only once using the last `d_c` parameter value,
    which is conventionally the most lenient parameter.

    Currently size factors are computed with min_size = 1 RFU counts.
    """
    output:
        "{test_data_prefix}/size_factors/mRFU.tsv",
    input:
        query_rfu_counts_tsv="{{test_data_prefix}}/rfu_counts/mRFU/{d_c_param}/{min_size_param}.tsv.gz".format(
            d_c_param=list(cdp_paramspace.instance_patterns)[-1],
            min_size_param=NO_MIN_SIZE_FILTERING,
        ),
        ref_rfu_counts_tsv=lambda wildcards: (
            replace_root_with_data(wildcards.test_data_prefix)
            + "/rfu_counts/mRFU/{d_c_param}/{min_size_param}.tsv.gz"
        ).format(
            d_c_param=list(cdp_paramspace.instance_patterns)[-1],
            min_size_param=NO_MIN_SIZE_FILTERING,
        ),
        ref_size_factors_tsv=lambda wildcards: (
            replace_root_with_data(wildcards.test_data_prefix)
            + "/size_factors/mRFU.tsv"
        ),
    benchmark:
        "{test_data_prefix}/size_factors/mRFU.tsv.benchmark"
    # According to https://github.com/serumdetect/snakemake-resource-usage/blob/f77fdcdf983e21bce2c97ede4dadeb4a34718e19/notebooks/analyze-input-sizes-vs-memory-use.2024-05-28.py.ipynb,
    # this rule uses below 4 GB of memory.
    resources:
        mem_mb=4000,
    conda:
        "../envs/query_feature_generation/query_size_factors.yaml"
    script:
        "../scripts/query_feature_generation/query_size_factors.R"


rule compute_query_variance_stabilized_counts:
    """Compute variance stabilized counts using transformGamPoi."""
    output:
        output_rfu_counts_tsv=f"{{test_data_prefix}}/pearson_vst_rfu_counts/mRFU/{cdp_paramspace.wildcard_pattern}/{rfu_size_paramspace.wildcard_pattern}.tsv.gz",
    input:
        rfu_counts_tsv=f"{{test_data_prefix}}/rfu_counts/mRFU/{cdp_paramspace.wildcard_pattern}/{rfu_size_paramspace.wildcard_pattern}.tsv.gz",
        size_factors_tsv=f"{{test_data_prefix}}/size_factors/mRFU.tsv",
        reference_mu_dispersion_tsv=lambda wildcards: (
            replace_root_with_data(wildcards.test_data_prefix)
            + f"/pearson_vst_mu_dispersion/mRFU/{cdp_paramspace.wildcard_pattern}/{rfu_size_paramspace.wildcard_pattern}.tsv.gz"
        ),
    benchmark:
        f"{{test_data_prefix}}/pearson_vst_rfu_counts/mRFU/{cdp_paramspace.wildcard_pattern}/{rfu_size_paramspace.wildcard_pattern}.tsv.gz.benchmark"
    resources:
        mem_mb=lambda _, input, attempt: 2000
        + size_mb(input[0]) * 100 * (0.5 + attempt / 2),
    conda:
        "../envs/query_feature_generation/compute_query_variance_stabilized_counts.yaml"
    script:
        "../scripts/query_feature_generation/compute_query_variance_stabilized_counts.R"


# Rule to adjust query RFU counts for covariates.
use rule adjust_rfu_counts_for_covariates as adjust_query_rfu_counts_for_covariates with:
    output:
        f"{{test_data_prefix}}/adjusted_rfu_counts/{{stats_type}}~{{predictor}}/mRFU/glm_adjustment~{{glm_adjustment}}/{cdp_paramspace.wildcard_pattern}/{rfu_filtering_paramspace.wildcard_pattern}.tsv.gz",
    input:
        rfu_counts_tsv=lambda wildcards: (
            f"{{test_data_prefix}}/norm_rfu_counts/mRFU/{cdp_paramspace.wildcard_pattern}/{rfu_size_paramspace.wildcard_pattern}.tsv.gz"
            if wildcards.glm_adjustment == "1"
            else f"{{test_data_prefix}}/rfu_counts/mRFU/{cdp_paramspace.wildcard_pattern}/{rfu_size_paramspace.wildcard_pattern}.tsv.gz"
        ),
        size_factors_tsv=f"{{test_data_prefix}}/size_factors/mRFU.tsv",
        metadata_tsv=f"{{test_data_prefix}}/parsed_metadata/mRFU.metadata.tsv",
        deseq2_obj=lambda wildcards: (
            replace_root_with_data(wildcards.test_data_prefix)
            + f"/deseq2_objs~{{predictor}}/mRFU/{cdp_paramspace.wildcard_pattern}/{rfu_size_paramspace.wildcard_pattern}.rds"
        ),
        glm_stats_tsv=lambda wildcards: (
            replace_root_with_data(wildcards.test_data_prefix)
            + f"/rfu_filtering_results/{{stats_type}}~{{predictor}}/mRFU/{cdp_paramspace.wildcard_pattern}/{rfu_filtering_paramspace.wildcard_pattern}.tsv.gz"
        ),
    benchmark:
        f"{{test_data_prefix}}/adjusted_rfu_counts/{{stats_type}}~{{predictor}}/mRFU/glm_adjustment~{{glm_adjustment}}/{cdp_paramspace.wildcard_pattern}/{rfu_filtering_paramspace.wildcard_pattern}.tsv.gz.benchmark"
    # According to https://github.com/serumdetect/snakemake-resource-usage/blob/f77fdcdf983e21bce2c97ede4dadeb4a34718e19/notebooks/analyze-input-sizes-vs-memory-use.2024-05-28.py.ipynb,
    # this rule always uses ~3000x input file size.
    resources:
        mem_mb=lambda _, input, attempt: 3000
        * size_mb(input.rfu_counts_tsv)
        * (0.5 + attempt / 2),


# Rule to adjust query RFU counts for covariates with blinded covariates.
use rule adjust_query_rfu_counts_for_covariates as adjust_query_rfu_counts_for_covariates_drop_covariates with:
    output:
        f"{{test_data_prefix}}/drop_covariates~{{fdr_cutoff}}~{{drop_covariates}}/adjusted_rfu_counts/{{stats_type}}~{{predictor}}/mRFU/glm_adjustment~{{glm_adjustment}}/{cdp_paramspace.wildcard_pattern}/{rfu_filtering_paramspace.wildcard_pattern}.tsv.gz",
    input:
        rfu_counts_tsv=lambda wildcards: (
            f"{{test_data_prefix}}/norm_rfu_counts/mRFU/{cdp_paramspace.wildcard_pattern}/{rfu_size_paramspace.wildcard_pattern}.tsv.gz"
            if wildcards.glm_adjustment == "1"
            else f"{{test_data_prefix}}/rfu_counts/mRFU/{cdp_paramspace.wildcard_pattern}/{rfu_size_paramspace.wildcard_pattern}.tsv.gz"
        ),
        size_factors_tsv=f"{{test_data_prefix}}/size_factors/mRFU.tsv",
        metadata_tsv=f"{{test_data_prefix}}/parsed_metadata/mRFU.metadata.tsv",
        deseq2_obj=lambda wildcards: (
            replace_root_with_data(wildcards.test_data_prefix)
            + f"/drop_covariates~{{drop_covariates}}/deseq2_objs~{{predictor}}/mRFU/{cdp_paramspace.wildcard_pattern}/{rfu_size_paramspace.wildcard_pattern}.rds"
        ),
        glm_stats_tsv=lambda wildcards: (
            replace_root_with_data(wildcards.test_data_prefix)
            + f"/drop_covariates~{{fdr_cutoff}}~{{drop_covariates}}/rfu_filtering_results/{{stats_type}}~{{predictor}}/mRFU/{cdp_paramspace.wildcard_pattern}/{rfu_filtering_paramspace.wildcard_pattern}.tsv.gz"
        ),
    benchmark:
        f"{{test_data_prefix}}/drop_covariates~{{fdr_cutoff}}~{{drop_covariates}}/adjusted_rfu_counts/{{stats_type}}~{{predictor}}/mRFU/glm_adjustment~{{glm_adjustment}}/{cdp_paramspace.wildcard_pattern}/{rfu_filtering_paramspace.wildcard_pattern}.tsv.gz.benchmark"
    # According to https://github.com/serumdetect/snakemake-resource-usage/blob/f77fdcdf983e21bce2c97ede4dadeb4a34718e19/notebooks/analyze-input-sizes-vs-memory-use.2024-05-28.py.ipynb,
    # this rule always uses ~3000x input file size.
    resources:
        mem_mb=lambda _, input, attempt: 3000
        * size_mb(input.rfu_counts_tsv)
        * (0.5 + attempt / 2),


if config["batch_jobs"]:

    # Batched rule to adjust query RFU counts for covariates. Only implemented
    # for log_lr features.
    use rule adjust_rfu_counts_for_covariates_batch as adjust_query_rfu_counts_for_covariates_batch with:
        output:
            expand(
                "{{test_data_prefix}}/adjusted_rfu_counts/{{stats_type}}~{{predictor}}/mRFU/glm_adjustment~log_lr/{cdp_paramspace_wildcard_pattern}/{rfu_filtering_paramspace_instance_patterns}.tsv.gz",
                cdp_paramspace_wildcard_pattern=[
                    cdp_paramspace.wildcard_pattern
                ],
                rfu_filtering_paramspace_instance_patterns=rfu_filtering_paramspace.instance_patterns,
            ),
        input:
            rfu_counts_tsv=lambda wildcards: _adjust_rfu_counts_for_covariates_batch_files(
                "{{{{test_data_prefix}}}}/rfu_counts/mRFU/{cdp_paramspace_wildcard_pattern}/{{rfu_size_paramspace_instance_patterns}}.tsv.gz".format(
                    cdp_paramspace_wildcard_pattern=_escape_curly_braces(
                        cdp_paramspace.wildcard_pattern
                    ),
                )
            ),
            size_factors_tsv="{test_data_prefix}/size_factors/mRFU.tsv",
            metadata_tsv="{test_data_prefix}/parsed_metadata/mRFU.metadata.tsv",
            deseq2_obj=lambda wildcards: _adjust_rfu_counts_for_covariates_batch_files(
                replace_root_with_data(wildcards.test_data_prefix)
                + "/deseq2_objs~{{{{predictor}}}}/mRFU/{cdp_paramspace_wildcard_pattern}/{{rfu_size_paramspace_instance_patterns}}.rds".format(
                    cdp_paramspace_wildcard_pattern=_escape_curly_braces(
                        cdp_paramspace.wildcard_pattern
                    ),
                )
            ),
            glm_stats_tsv=lambda wildcards: _adjust_rfu_counts_for_covariates_batch_files(
                replace_root_with_data(wildcards.test_data_prefix)
                + "/rfu_filtering_results/{{{{stats_type}}}}~{{{{predictor}}}}/mRFU/{cdp_paramspace_wildcard_pattern}/{{rfu_filtering_paramspace_instance_patterns}}.tsv.gz".format(
                    cdp_paramspace_wildcard_pattern=_escape_curly_braces(
                        cdp_paramspace.wildcard_pattern
                    ),
                )
            ),
        benchmark:
            f"{{test_data_prefix}}/adjusted_rfu_counts/{{stats_type}}~{{predictor}}/mRFU/glm_adjustment~log_lr/{cdp_paramspace.wildcard_pattern}.tsv.gz.benchmark"
        # According to https://github.com/serumdetect/snakemake-resource-usage/blob/339dfbc5eb336a9c8ddab1796a33f685a172e31e/notebooks/analyze-input-sizes-vs-memory-use.2024-05-28.py.ipynb,
        # this rule always uses ~3000x input file size. Give it a bit more
        # memory since we have to iterate through all samples.
        resources:
            mem_mb=lambda _, input, attempt: 4500
            * max(size_mb(fname) for fname in input.rfu_counts_tsv)
            * (0.5 + attempt / 2),

    ruleorder: adjust_query_rfu_counts_for_covariates_batch > adjust_query_rfu_counts_for_covariates

    # Batched rule to adjust query RFU counts for covariates with blinded
    # covariates. Only implemented for log_lr features.
    use rule adjust_query_rfu_counts_for_covariates_batch as adjust_query_rfu_counts_for_covariates_batch_drop_covariates with:
        output:
            expand(
                "{{test_data_prefix}}/drop_covariates~{{drop_covariates}}/adjusted_rfu_counts/{{stats_type}}~{{predictor}}/mRFU/glm_adjustment~log_lr/{cdp_paramspace_wildcard_pattern}/{rfu_filtering_paramspace_instance_patterns}.tsv.gz",
                cdp_paramspace_wildcard_pattern=[
                    cdp_paramspace.wildcard_pattern
                ],
                rfu_filtering_paramspace_instance_patterns=rfu_filtering_paramspace.instance_patterns,
            ),
        input:
            rfu_counts_tsv=lambda wildcards: _adjust_rfu_counts_for_covariates_batch_files(
                "{{{{test_data_prefix}}}}/rfu_counts/mRFU/{cdp_paramspace_wildcard_pattern}/{{rfu_size_paramspace_instance_patterns}}.tsv.gz".format(
                    cdp_paramspace_wildcard_pattern=_escape_curly_braces(
                        cdp_paramspace.wildcard_pattern
                    ),
                )
            ),
            size_factors_tsv="{test_data_prefix}/size_factors/mRFU.tsv",
            metadata_tsv="{test_data_prefix}/parsed_metadata/mRFU.metadata.tsv",
            deseq2_obj=lambda wildcards: _adjust_rfu_counts_for_covariates_batch_files(
                replace_root_with_data(wildcards.test_data_prefix)
                + "/drop_covariates~{{{{drop_covariates}}}}/deseq2_objs~{{{{predictor}}}}/mRFU/{cdp_paramspace_wildcard_pattern}/{{rfu_size_paramspace_instance_patterns}}.rds".format(
                    cdp_paramspace_wildcard_pattern=_escape_curly_braces(
                        cdp_paramspace.wildcard_pattern
                    ),
                )
            ),
            glm_stats_tsv=lambda wildcards: _adjust_rfu_counts_for_covariates_batch_files(
                replace_root_with_data(wildcards.test_data_prefix)
                + "/drop_covariates~{{{{drop_covariates}}}}/rfu_filtering_results/{{{{stats_type}}}}~{{{{predictor}}}}/mRFU/{cdp_paramspace_wildcard_pattern}/{{rfu_filtering_paramspace_instance_patterns}}.tsv.gz".format(
                    cdp_paramspace_wildcard_pattern=_escape_curly_braces(
                        cdp_paramspace.wildcard_pattern
                    ),
                )
            ),
        benchmark:
            f"{{test_data_prefix}}/drop_covariates~{{drop_covariates}}/adjusted_rfu_counts/{{stats_type}}~{{predictor}}/mRFU/glm_adjustment~log_lr/{cdp_paramspace.wildcard_pattern}.tsv.gz.benchmark"

    ruleorder: adjust_query_rfu_counts_for_covariates_batch_drop_covariates > adjust_query_rfu_counts_for_covariates_drop_covariates


use rule prepare_ml_data as prepare_query_ml_data with:
    output:
        f"{{test_data_prefix}}/adjusted_rfu_counts/{{stats_type}}~{{predictor}}/mRFU/glm_adjustment~{{glm_adjustment}}/{rfu_filtering_paramspace.wildcard_pattern}/ml_features.tsv.gz",
    input:
        rfu_counts=lambda wildcards: [
            fname.replace("{prefix}", "{test_data_prefix}")
            for fname in ml_features_input(wildcards)
        ],
        glm_stats_tsv=lambda wildcards: expand(
            replace_root_with_data(wildcards.test_data_prefix)
            + "/rfu_filtering_results/{{stats_type}}~{{predictor}}/mRFU/{params}/{rfu_filtering_wildcard_pattern}.tsv.gz",
            params=list(cdp_paramspace.instance_patterns),
            rfu_filtering_wildcard_pattern=rfu_filtering_paramspace.wildcard_pattern,
        ),
    benchmark:
        f"{{test_data_prefix}}/adjusted_rfu_counts/{{stats_type}}~{{predictor}}/mRFU/glm_adjustment~{{glm_adjustment}}/{rfu_filtering_paramspace.wildcard_pattern}/ml_features.tsv.gz.benchmark"


# Compile ML input features with dropped covariates.
use rule prepare_query_ml_data as prepare_query_ml_data_drop_covariates with:
    output:
        f"{{test_data_prefix}}/drop_covariates~{{fdr_cutoff}}~{{drop_covariates}}/adjusted_rfu_counts/{{stats_type}}~{{predictor}}/mRFU/glm_adjustment~{{glm_adjustment}}/{rfu_filtering_paramspace.wildcard_pattern}/ml_features.tsv.gz",
    input:
        rfu_counts=expand(
            "{{test_data_prefix}}/drop_covariates~{{fdr_cutoff}}~{{drop_covariates}}/adjusted_rfu_counts/{{stats_type}}~{{predictor}}/mRFU/glm_adjustment~{{glm_adjustment}}/{params}/{rfu_filtering_wildcard_pattern}.tsv.gz",
            params=list(cdp_paramspace.instance_patterns),
            rfu_filtering_wildcard_pattern=rfu_filtering_paramspace.wildcard_pattern,
        ),
        glm_stats_tsv=lambda wildcards: expand(
            replace_root_with_data(wildcards.test_data_prefix)
            + "/drop_covariates~{{fdr_cutoff}}~{{drop_covariates}}/rfu_filtering_results/{{stats_type}}~{{predictor}}/mRFU/{params}/{rfu_filtering_wildcard_pattern}.tsv.gz",
            params=list(cdp_paramspace.instance_patterns),
            rfu_filtering_wildcard_pattern=rfu_filtering_paramspace.wildcard_pattern,
        ),
    benchmark:
        f"{{test_data_prefix}}/drop_covariates~{{fdr_cutoff}}~{{drop_covariates}}/adjusted_rfu_counts/{{stats_type}}~{{predictor}}/mRFU/glm_adjustment~{{glm_adjustment}}/{rfu_filtering_paramspace.wildcard_pattern}/ml_features.tsv.gz.benchmark"

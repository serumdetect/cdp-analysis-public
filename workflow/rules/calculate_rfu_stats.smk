"""Calculates stats pertaining to RFUs"""


rule calculate_vj_fractions:
    """Calculates V-/J-gene fractions of clonotypes per RFU."""
    output:
        f"{{any_prefix}}/rfu_stats/mRFU/{cdp_paramspace.wildcard_pattern}/{rfu_size_paramspace.wildcard_pattern}/vj_fractions.tsv.gz",
    input:
        f"{{any_prefix}}/annotated_tcrs/mRFU/{cdp_paramspace.wildcard_pattern}.ftr",
    benchmark:
        f"{{any_prefix}}/rfu_stats/mRFU/{cdp_paramspace.wildcard_pattern}/{rfu_size_paramspace.wildcard_pattern}/vj_fractions.tsv.gz.benchmark"
    resources:
        mem_mb=lambda _, input, attempt: size_mb(input[0])
        * 3
        * (0.5 + attempt / 2),
    conda:
        "../envs/calculate_rfu_stats/rfu_stats.yaml"
    script:
        "../scripts/calculate_rfu_stats/calculate_vj_fractions.py"


rule combine_vj_fractions:
    """Combines V-/J-gene fractions across cdp parameter values."""
    output:
        f"{{any_prefix}}/rfu_stats/mRFU/{rfu_size_paramspace.wildcard_pattern}/vj_fractions.tsv.gz",
    input:
        lambda wildcards: expand(
            "{{any_prefix}}/rfu_stats/mRFU/{cdp_paramspace_instance}/{rfu_size_paramspace_wildcard_pattern}/vj_fractions.tsv.gz",
            cdp_paramspace_instance=cdp_paramspace.instance_patterns,
            rfu_size_paramspace_wildcard_pattern=[
                rfu_size_paramspace.wildcard_pattern
            ],
        ),
    params:
        rfu_prefixes=list(cdp_paramspace.instance_patterns),
    conda:
        "../envs/calculate_rfu_stats/rfu_stats.yaml"
    script:
        "../scripts/calculate_rfu_stats/combine_vj_fractions.py"


rule gather_annotated_clusters_with_rfu_filtering:
    """Gather cluster stats for RFUs that pass filtering into an .xlsx file."""
    group:
        "rfu_stats"
    output:
        f"{{prefix_no_drop_covariates}}/rfu_filtering_results/{{stats_type}}~{{predictor}}/mRFU/{rfu_filtering_paramspace.wildcard_pattern}/annotated_clusters_with_rfu_filtering.tsv.gz",
    input:
        rfu_stats=expand(
            "{{prefix_no_drop_covariates}}/rfu_filtering_results/{{stats_type}}~{{predictor}}/mRFU/{cdp_paramspace}/{rfu_filtering_paramspace_wildcard_pattern}.tsv.gz",
            cdp_paramspace=cdp_paramspace.instance_patterns,
            rfu_filtering_paramspace_wildcard_pattern=[
                rfu_filtering_paramspace.wildcard_pattern
            ],
        ),
        cluster_stats=expand(
            "{{prefix_no_drop_covariates}}/cluster_stats~{{predictor}}/mRFU/{cdp_paramspace}.tsv.gz",
            cdp_paramspace=cdp_paramspace.instance_patterns,
        ),
    benchmark:
        f"{{prefix_no_drop_covariates}}/rfu_filtering_results/{{stats_type}}~{{predictor}}/mRFU/{rfu_filtering_paramspace.wildcard_pattern}/annotated_clusters_with_rfu_filtering.xlsx.benchmark"
    params:
        d_c_params=list(cdp_paramspace.instance_patterns),
    resources:
        mem_mb=lambda _, input, attempt: sum(
            [size_mb(input_file) for input_file in input.cluster_stats]
        )
        * 10
        * (0.5 + attempt / 2),
    conda:
        "../envs/calculate_rfu_stats/rfu_stats.yaml"
    script:
        "../scripts/calculate_rfu_stats/gather_annotated_clusters_with_rfu_filtering.py"


# Gather results for RFU filtering with covariates dropped.
use rule gather_annotated_clusters_with_rfu_filtering as gather_annotated_clusters_with_rfu_filtering_drop_covariates with:
    group:
        "rfu_stats"
    output:
        f"{{prefix_no_drop_covariates}}/drop_covariates~{{fdr_cutoff}}~{{drop_covariates}}/rfu_filtering_results/{{stats_type}}~{{predictor}}/mRFU/{rfu_filtering_paramspace.wildcard_pattern}/annotated_clusters_with_rfu_filtering.tsv.gz",
    input:
        rfu_stats=expand(
            "{{prefix_no_drop_covariates}}/drop_covariates~{{fdr_cutoff}}~{{drop_covariates}}/rfu_filtering_results/{{stats_type}}~{{predictor}}/mRFU/{cdp_paramspace}/{rfu_filtering_paramspace_wildcard_pattern}.tsv.gz",
            cdp_paramspace=cdp_paramspace.instance_patterns,
            rfu_filtering_paramspace_wildcard_pattern=[
                rfu_filtering_paramspace.wildcard_pattern
            ],
        ),
        cluster_stats=expand(
            "{{prefix_no_drop_covariates}}/cluster_stats~{{predictor}}/mRFU/{cdp_paramspace}.tsv.gz",
            cdp_paramspace=cdp_paramspace.instance_patterns,
        ),
    benchmark:
        f"{{prefix_no_drop_covariates}}/drop_covariates~{{fdr_cutoff}}~{{drop_covariates}}/rfu_filtering_results/{{stats_type}}~{{predictor}}/mRFU/{rfu_filtering_paramspace.wildcard_pattern}/annotated_clusters_with_rfu_filtering.xlsx.benchmark"
    params:
        d_c_params=list(cdp_paramspace.instance_patterns),


rule join_rfu_stats:
    """Joins additional RFU stats to GLM stats."""
    group:
        "rfu_stats"
    output:
        f"{{any_prefix}}/{{stats_type}}~{{predictor}}/mRFU/{rfu_filtering_paramspace.wildcard_pattern}/rfu_stats.tsv.gz",
    input:
        glm_stats=f"{{any_prefix}}/rfu_filtering_results/{{stats_type}}~{{predictor}}/mRFU/{rfu_filtering_paramspace.wildcard_pattern}/annotated_clusters_with_rfu_filtering.tsv.gz",
        vj_fractions=f"{{any_prefix}}/rfu_stats/mRFU/{rfu_size_paramspace.wildcard_pattern}/vj_fractions.tsv.gz",
    conda:
        "../envs/calculate_rfu_stats/rfu_stats.yaml"
    script:
        "../scripts/calculate_rfu_stats/join_rfu_stats.py"


use rule join_rfu_stats as join_rfu_stats_drop_covariates with:
    group:
        "rfu_stats"
    output:
        f"{{any_prefix}}/drop_covariates~{{fdr_cutoff_pat}}~{{drop_covariates}}/{{stats_type}}~{{predictor}}/mRFU/{rfu_filtering_paramspace.wildcard_pattern}/rfu_stats.tsv.gz",
    input:
        glm_stats=f"{{any_prefix}}/drop_covariates~{{fdr_cutoff_pat}}~{{drop_covariates}}/rfu_filtering_results/{{stats_type}}~{{predictor}}/mRFU/{rfu_filtering_paramspace.wildcard_pattern}/annotated_clusters_with_rfu_filtering.tsv.gz",
        vj_fractions=f"{{any_prefix}}/rfu_stats/mRFU/{rfu_size_paramspace.wildcard_pattern}/vj_fractions.tsv.gz",

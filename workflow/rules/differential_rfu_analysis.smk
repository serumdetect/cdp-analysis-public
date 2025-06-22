"""
Various Snakemake rules for running RFU statistics with various settings.

These rules for major branching points in the main project workflow. The base
rule `compute_deseq2_stats` is used to run standard DESeq2 with no TCR filters.
The inherited rules are as follows.
"""

import math
import re


################################################################################
# Computation of RFU counts matrices and size factors.
################################################################################


rule compute_rfu_counts_matrix:
    """Compute a matrix of RFU counts"""
    group:
        "rfu_counts"
    output:
        f"{{prefix}}/rfu_counts/mRFU/{cdp_paramspace.wildcard_pattern}/{rfu_size_paramspace.wildcard_pattern}.tsv.gz",
    input:
        f"{{prefix}}/annotated_tcrs/mRFU/{cdp_paramspace.wildcard_pattern}.ftr",
    benchmark:
        f"{{prefix}}/rfu_counts/mRFU/{cdp_paramspace.wildcard_pattern}/{rfu_size_paramspace.wildcard_pattern}.tsv.gz.benchmark"
    # According to https://github.com/serumdetect/snakemake-resource-usage/blob/f77fdcdf983e21bce2c97ede4dadeb4a34718e19/notebooks/analyze-input-sizes-vs-memory-use.2024-05-28.py.ipynb,
    # this rule uses 10 * input file size GB of memory regardless of input size.
    # A fraction of these jobs will fail but upon retrying with 15 * input file
    # size GB of memory, they will succeed.
    #
    # 2024-10-31: Reducing memory requirement after transitioning to Polars.
    resources:
        mem_mb=lambda _, input, attempt: size_mb(input[0]) * (0.5 + attempt / 2),
    retries: 3
    params:
        count_templates=config["glm"]["count_templates"],
        min_subject_count=config["glm"]["min_subject_count"],
        min_tcr_size=lambda wildcards: int(wildcards.min_size),
    conda:
        "../envs/differential_rfu_analysis/compute_rfu_counts_matrix.yaml"
    script:
        "../scripts/differential_rfu_analysis/compute_rfu_counts_matrix.py"


if not config["glm"]["count_templates"]:

    rule compute_size_factors:
        """Compute size factors using a big clustering setting.

        Size factors are computed using the last CDP setting, which conventionally
        is the most lenient clustering setting generating the largest RFUs.
        """
        group:
            "rfu_counts"
        output:
            "{prefix}/size_factors/mRFU.tsv",
        input:
            # Default to using min_size~1 for computing size factors.
            "{{prefix}}/rfu_counts/mRFU/{d_c_param}/{rfu_size_param}.tsv.gz".format(
                d_c_param=list(cdp_paramspace.instance_patterns)[-1],
                rfu_size_param=NO_MIN_SIZE_FILTERING,
            ),
        benchmark:
            "{prefix}/size_factors/mRFU.tsv.benchmark"
        resources:
            mem_mb=5000,
        conda:
            "../envs/differential_rfu_analysis/compute_size_factors.yaml"
        script:
            "../scripts/differential_rfu_analysis/compute_size_factors.R"

else:

    rule compute_size_factors:
        """Compute size factors naively as the total template count."""
        group:
            "rfu_counts"
        output:
            "{prefix}/size_factors/mRFU.tsv",
        input:
            # Default to using min_size~1 for computing size factors.
            "{{prefix}}/rfu_counts/mRFU/{d_c_param}/{rfu_size_param}.tsv.gz".format(
                d_c_param=list(cdp_paramspace.instance_patterns)[-1],
                rfu_size_param=NO_MIN_SIZE_FILTERING,
            ),
        benchmark:
            "{prefix}/size_factors/mRFU.tsv.benchmark"
        resources:
            mem_mb=5000,
        conda:
            "../envs/differential_rfu_analysis/compute_naive_size_factors.yaml"
        script:
            "../scripts/differential_rfu_analysis/compute_naive_size_factors.py"


rule normalize_rfu_counts:
    """Normalize RFU counts by size factors."""
    group:
        "rfu_counts"
    output:
        f"{{any_prefix}}/norm_rfu_counts/mRFU/{cdp_paramspace.wildcard_pattern}/{rfu_size_paramspace.wildcard_pattern}.tsv.gz",
    input:
        rfu_counts_tsv=f"{{any_prefix}}/rfu_counts/mRFU/{cdp_paramspace.wildcard_pattern}/{rfu_size_paramspace.wildcard_pattern}.tsv.gz",
        size_factors_tsv=f"{{any_prefix}}/size_factors/mRFU.tsv",
    benchmark:
        f"{{any_prefix}}/norm_rfu_counts/mRFU/{cdp_paramspace.wildcard_pattern}/{rfu_size_paramspace.wildcard_pattern}.tsv.gz.benchmark"
    resources:
        mem_mb=lambda _, input, attempt: 2000
        + size_mb(input[0]) * 5 * (0.5 + attempt / 2),
    conda:
        "../envs/differential_rfu_analysis/normalize_rfu_counts.yaml"
    script:
        "../scripts/differential_rfu_analysis/normalize_rfu_counts.py"


rule compute_variance_stabilized_counts:
    """Compute variance stabilized counts using transformGamPoi."""
    output:
        output_rfu_counts_tsv=f"{{prefix_no_drop_covariates}}/pearson_vst_rfu_counts/mRFU/{cdp_paramspace.wildcard_pattern}/{rfu_size_paramspace.wildcard_pattern}.tsv.gz",
        output_mu_dispersion_tsv=f"{{prefix_no_drop_covariates}}/pearson_vst_mu_dispersion/mRFU/{cdp_paramspace.wildcard_pattern}/{rfu_size_paramspace.wildcard_pattern}.tsv.gz",
    input:
        rfu_counts_tsv=f"{{prefix_no_drop_covariates}}/rfu_counts/mRFU/{cdp_paramspace.wildcard_pattern}/{rfu_size_paramspace.wildcard_pattern}.tsv.gz",
        size_factors_tsv=f"{{prefix_no_drop_covariates}}/size_factors/mRFU.tsv",
    benchmark:
        f"{{prefix_no_drop_covariates}}/pearson_vst_rfu_counts/mRFU/{cdp_paramspace.wildcard_pattern}/{rfu_size_paramspace.wildcard_pattern}.tsv.gz.benchmark"
    params:
        min_subject_count=config["glm"]["min_subject_count"],
    resources:
        mem_mb=lambda _, input, attempt: 2000
        + size_mb(input[0]) * 600 * (0.5 + attempt / 2),
    conda:
        "../envs/differential_rfu_analysis/compute_variance_stabilized_counts.yaml"
    script:
        "../scripts/differential_rfu_analysis/compute_variance_stabilized_counts.R"


################################################################################
# Computing cancer-associated RFUs using DESeq2.
################################################################################


rule compute_deseq2_stats:
    """Compute differential RFU analysis using DESeq2"""
    output:
        deseq2_obj=f"{{prefix_no_drop_covariates}}/deseq2_objs~{{predictor}}/mRFU/{cdp_paramspace.wildcard_pattern}/{rfu_size_paramspace.wildcard_pattern}.rds",
        glm_stats_tsv=f"{{prefix_no_drop_covariates}}/glm_stats~{{predictor}}/mRFU/{cdp_paramspace.wildcard_pattern}/{rfu_size_paramspace.wildcard_pattern}.tsv.gz",
        coef_mat_tsv=f"{{prefix_no_drop_covariates}}/glm_coef_mats/glm_stats~{{predictor}}/mRFU/{cdp_paramspace.wildcard_pattern}/{rfu_size_paramspace.wildcard_pattern}.tsv.gz",
    input:
        counts_mat=f"{{prefix_no_drop_covariates}}/rfu_counts/mRFU/{cdp_paramspace.wildcard_pattern}/{rfu_size_paramspace.wildcard_pattern}.tsv.gz",
        size_factors="{prefix_no_drop_covariates}/size_factors/mRFU.tsv",
        metadata="{prefix_no_drop_covariates}/parsed_metadata/mRFU.metadata.tsv",
    benchmark:
        f"{{prefix_no_drop_covariates}}/deseq2_objs~{{predictor}}/mRFU/{cdp_paramspace.wildcard_pattern}/{rfu_size_paramspace.wildcard_pattern}.rds.benchmark"
    # According to https://github.com/serumdetect/snakemake-resource-usage/blob/f77fdcdf983e21bce2c97ede4dadeb4a34718e19/notebooks/analyze-input-sizes-vs-memory-use.2024-05-28.py.ipynb,
    # this rule uses roughly 1500 * input file size.
    resources:
        mem_mb=lambda _, input, attempt: size_mb(input.counts_mat)
        * 2250
        * (0.5 + attempt / 2),
    params:
        min_subject_count=config["glm"]["min_subject_count"],
        depth_colnames=config["glm"]["depth_colnames"],
        depth_poly_order=config["glm"]["depth_poly_order"],
        covariates=get_glm_covariates,
    conda:
        "../envs/differential_rfu_analysis/compute_deseq2_stats.yaml"
    script:
        "../scripts/differential_rfu_analysis/compute_deseq2_stats.R"


# Compute DESeq2 stats with some covariates removed from the GLM model.
use rule compute_deseq2_stats as compute_deseq2_stats_drop_covariates with:
    output:
        deseq2_obj=f"{{prefix_no_drop_covariates}}/drop_covariates~{{drop_covariates}}/deseq2_objs~{{predictor}}/mRFU/{cdp_paramspace.wildcard_pattern}/{rfu_size_paramspace.wildcard_pattern}.rds",
        glm_stats_tsv=f"{{prefix_no_drop_covariates}}/drop_covariates~{{drop_covariates}}/glm_stats~{{predictor}}/mRFU/{cdp_paramspace.wildcard_pattern}/{rfu_size_paramspace.wildcard_pattern}.tsv.gz",
        coef_mat_tsv=f"{{prefix_no_drop_covariates}}/drop_covariates~{{drop_covariates}}/glm_coef_mats/glm_stats~{{predictor}}/mRFU/{cdp_paramspace.wildcard_pattern}/{rfu_size_paramspace.wildcard_pattern}.tsv.gz",
    benchmark:
        f"{{prefix_no_drop_covariates}}/drop_covariates~{{drop_covariates}}/deseq2_objs~{{predictor}}/mRFU/{cdp_paramspace.wildcard_pattern}/{rfu_size_paramspace.wildcard_pattern}.rds.benchmark"
    params:
        min_subject_count=config["glm"]["min_subject_count"],
        depth_colnames=config["glm"]["depth_colnames"],
        depth_poly_order=config["glm"]["depth_poly_order"],
        disease_subset=None,
        covariates=lambda wildcards: [
            covariate
            for covariate in get_glm_covariates(wildcards)
            if covariate not in wildcards.drop_covariates.split("__")
        ],


ruleorder: compute_deseq2_stats_drop_covariates > compute_deseq2_stats


rule compute_nonparametric_stats:
    """Compute differential RFU analysis using Wilcoxon rank-sum test."""
    group:
        "rfu_counts"
    output:
        f"{{prefix_no_drop_covariates}}/{{stats_type}}_stats~{{predictor}}/mRFU/{cdp_paramspace.wildcard_pattern}/{rfu_size_paramspace.wildcard_pattern}.tsv.gz",
    input:
        norm_rfu_counts_tsv=f"{{prefix_no_drop_covariates}}/norm_rfu_counts/mRFU/{cdp_paramspace.wildcard_pattern}/{rfu_size_paramspace.wildcard_pattern}.tsv.gz",
        metadata="{prefix_no_drop_covariates}/parsed_metadata/mRFU.metadata.tsv",
    benchmark:
        f"{{prefix_no_drop_covariates}}/{{stats_type}}_stats~{{predictor}}/mRFU/{cdp_paramspace.wildcard_pattern}/{rfu_size_paramspace.wildcard_pattern}.rds.benchmark"
    # According to https://github.com/serumdetect/snakemake-resource-usage/blob/782943ec7fd1ba5fc448e1719c8e5a5da1fa577e/notebooks/analyze-input-sizes-vs-memory-use.2024-05-28.py.ipynb,
    # this rule uses roughly 100 * input file size. File size got reduced by a
    # smaller input matrix.
    resources:
        mem_mb=lambda _, input, attempt: 1000
        + size_mb(input.norm_rfu_counts_tsv) * 5 * (0.5 + attempt / 2),
    wildcard_constraints:
        stats_type="wilcox|fisher|bm",
    conda:
        "../envs/differential_rfu_analysis/compute_nonparametric_stats.yaml"
    script:
        "../scripts/differential_rfu_analysis/compute_nonparametric_stats.py"

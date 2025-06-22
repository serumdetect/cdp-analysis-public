"""
Snakemake rules generating basic RFU features for ML by adjusting RFU counts for
fitted GLM models' covariates.
"""

import re

import itertools
import snakemake.io


def adjust_rfu_counts_for_covariates_threads(wildcards):
    """Only glm_adjustment~1 benefits from multithreading."""
    return config["threads"] if wildcards.glm_adjustment == "1" else 1


rule adjust_rfu_counts_for_covariates:
    """Compute normalized RFU counts adjusted for nuisance variables.

    "Adjusted count" is a now a misnomer due to a relic of the original purpose
    of this rule, which was to adjust RFU counts for covariates associated with
    it (except for cancer status).

    With `glm_adjustment~log_lr`, the "adjusted" counts are actually
    log-likelihood ratios. With `glm_adjustment~0`, raw counts are returned.
    This is included to simplify the usage of unadjusted counts in ML testing.
    """
    output:
        f"{{prefix_no_drop_covariates}}/adjusted_rfu_counts/{{stats_type}}~{{predictor}}/mRFU/glm_adjustment~{{glm_adjustment}}/{cdp_paramspace.wildcard_pattern}/{rfu_filtering_paramspace.wildcard_pattern}.tsv.gz",
    input:
        rfu_counts_tsv=lambda wildcards: (
            f"{{prefix_no_drop_covariates}}/norm_rfu_counts/mRFU/{cdp_paramspace.wildcard_pattern}/{rfu_size_paramspace.wildcard_pattern}.tsv.gz"
            if wildcards.glm_adjustment == "1"
            else f"{{prefix_no_drop_covariates}}/rfu_counts/mRFU/{cdp_paramspace.wildcard_pattern}/{rfu_size_paramspace.wildcard_pattern}.tsv.gz"
        ),
        size_factors_tsv="{prefix_no_drop_covariates}/size_factors/mRFU.tsv",
        metadata_tsv="{prefix_no_drop_covariates}/parsed_metadata/mRFU.metadata.tsv",
        deseq2_obj=f"{{prefix_no_drop_covariates}}/deseq2_objs~{{predictor}}/mRFU/{cdp_paramspace.wildcard_pattern}/{rfu_size_paramspace.wildcard_pattern}.rds",
        glm_stats_tsv=f"{{prefix_no_drop_covariates}}/rfu_filtering_results/{{stats_type}}~{{predictor}}/mRFU/{cdp_paramspace.wildcard_pattern}/{rfu_filtering_paramspace.wildcard_pattern}.tsv.gz",
    benchmark:
        f"{{prefix_no_drop_covariates}}/adjusted_rfu_counts/{{stats_type}}~{{predictor}}/mRFU/glm_adjustment~{{glm_adjustment}}/{cdp_paramspace.wildcard_pattern}/{rfu_filtering_paramspace.wildcard_pattern}.tsv.gz.benchmark"
    params:
        glm_adjustment=lambda wildcards: wildcards.glm_adjustment,
        fdr_cutoff=config["fdr_cutoff"],
        depth_colnames=config["glm"]["depth_colnames"],
        depth_poly_order=config["glm"]["depth_poly_order"],
    threads: adjust_rfu_counts_for_covariates_threads
    # According to https://github.com/serumdetect/snakemake-resource-usage/blob/f77fdcdf983e21bce2c97ede4dadeb4a34718e19/notebooks/analyze-input-sizes-vs-memory-use.2024-05-28.py.ipynb,
    # this rule uses below 1250x of input file size (sometimes more).
    resources:
        mem_mb=lambda _, input, attempt: 2500
        * size_mb(input.rfu_counts_tsv)
        * (0.5 + attempt / 2),
    conda:
        "../envs/feature_generation/adjust_rfu_counts.yaml"
    script:
        "../scripts/feature_generation/adjust_rfu_counts.R"


# Compute adjusted RFU counts some covariates removed from the GLM model.
use rule adjust_rfu_counts_for_covariates as adjust_rfu_counts_for_covariates_drop_covariates with:
    output:
        f"{{prefix_no_drop_covariates}}/drop_covariates~{{drop_covariates}}/adjusted_rfu_counts/{{stats_type}}~{{predictor}}/mRFU/glm_adjustment~{{glm_adjustment}}/{cdp_paramspace.wildcard_pattern}/{rfu_filtering_paramspace.wildcard_pattern}.tsv.gz",
    input:
        rfu_counts_tsv=lambda wildcards: (
            f"{{prefix_no_drop_covariates}}/norm_rfu_counts/mRFU/{cdp_paramspace.wildcard_pattern}/{rfu_size_paramspace.wildcard_pattern}.tsv.gz"
            if wildcards.glm_adjustment == "1"
            else f"{{prefix_no_drop_covariates}}/rfu_counts/mRFU/{cdp_paramspace.wildcard_pattern}/{rfu_size_paramspace.wildcard_pattern}.tsv.gz"
        ),
        size_factors_tsv="{prefix_no_drop_covariates}/size_factors/mRFU.tsv",
        metadata_tsv="{prefix_no_drop_covariates}/parsed_metadata/mRFU.metadata.tsv",
        deseq2_obj=f"{{prefix_no_drop_covariates}}/drop_covariates~{{drop_covariates}}/deseq2_objs~{{predictor}}/mRFU/{cdp_paramspace.wildcard_pattern}/{rfu_size_paramspace.wildcard_pattern}.rds",
        glm_stats_tsv=f"{{prefix_no_drop_covariates}}/drop_covariates~{{drop_covariates}}/rfu_filtering_results/{{stats_type}}~{{predictor}}/mRFU/{cdp_paramspace.wildcard_pattern}/{rfu_filtering_paramspace.wildcard_pattern}.tsv.gz",
    benchmark:
        f"{{prefix_no_drop_covariates}}/drop_covariates~{{drop_covariates}}/adjusted_rfu_counts/{{stats_type}}~{{predictor}}/mRFU/glm_adjustment~{{glm_adjustment}}/{cdp_paramspace.wildcard_pattern}/{rfu_filtering_paramspace.wildcard_pattern}.tsv.gz.benchmark"


if config["batch_jobs"]:
    # Use a rule batching over all cdp_paramspace and rfu_filtering_paramspace
    # values.

    def _escape_curly_braces(s):
        """Convert single curly braces to double curly braces."""
        return s.replace("{", "{{").replace("}", "}}")

    def _adjust_rfu_counts_for_covariates_batch_files(input_file_pattern):
        """Expand input/output file names based on input file pattern."""
        # Compile the RFU size wildcard pattern regex.
        rfu_size_wildcard_name = list(
            snakemake.io.get_wildcard_names(
                rfu_size_paramspace.wildcard_pattern
            )
        )[0]
        rfu_size_pat = re.compile(
            re.sub(
                r"\$$",
                "",
                snakemake.io.regex(
                    snakemake.io.update_wildcard_constraints(
                        rfu_size_paramspace.wildcard_pattern,
                        workflow.wildcard_constraints,
                        {},
                    )
                ),
            )
        )

        # Compile a list of rfu_size_paramspace instance patterns.
        rfu_size_paramspace_instance_patterns = [
            rfu_size_paramspace.wildcard_pattern.replace(
                "{" + rfu_size_wildcard_name + "}",
                re.match(rfu_size_pat, instance).group(1),
            )
            for instance in rfu_filtering_paramspace.instance_patterns
        ]

        return expand(
            input_file_pattern,
            rfu_filtering_paramspace_instance_patterns=list(
                rfu_filtering_paramspace.instance_patterns
            ),
            # Produce a list of size parameters corresponding to each
            # rfu_filtering_paramspace instance.
            rfu_size_paramspace_instance_patterns=rfu_size_paramspace_instance_patterns,
        )

    rule adjust_rfu_counts_for_covariates_batch:
        """Compute normalized RFU counts adjusted for nuisance variables.

        "Adjusted count" is a now a misnomer due to a relic of the original
        purpose of this rule, which was to adjust RFU counts for covariates
        associated with it (except for cancer status).

        With `glm_adjustment~log_lr`, the "adjusted" counts are actually
        log-likelihood ratios. With `glm_adjustment~0`, raw counts are returned.
        This is included to simplify the usage of unadjusted counts in ML
        testing.

        This is a batched version of the adjust_rfu_counts_for_covariates rule
        and is only implemented for glm_adjustment~log_lr, since this adjustment
        is fastest and hence benefits from batching the most.
        """
        output:
            expand(
                "{{prefix_no_drop_covariates}}/adjusted_rfu_counts/{{stats_type}}~{{predictor}}/mRFU/glm_adjustment~log_lr/{cdp_paramspace_wildcard_pattern}/{rfu_filtering_paramspace_instance_patterns}.tsv.gz",
                cdp_paramspace_wildcard_pattern=[
                    cdp_paramspace.wildcard_pattern
                ],
                rfu_filtering_paramspace_instance_patterns=rfu_filtering_paramspace.instance_patterns,
            ),
        input:
            rfu_counts_tsv=lambda wildcards: _adjust_rfu_counts_for_covariates_batch_files(
                "{{{{prefix_no_drop_covariates}}}}/rfu_counts/mRFU/{cdp_paramspace_wildcard_pattern}/{{rfu_size_paramspace_instance_patterns}}.tsv.gz".format(
                    cdp_paramspace_wildcard_pattern=_escape_curly_braces(
                        cdp_paramspace.wildcard_pattern
                    )
                )
            ),
            size_factors_tsv="{prefix_no_drop_covariates}/size_factors/mRFU.tsv",
            metadata_tsv="{prefix_no_drop_covariates}/parsed_metadata/mRFU.metadata.tsv",
            deseq2_obj=_adjust_rfu_counts_for_covariates_batch_files(
                "{{{{prefix_no_drop_covariates}}}}/deseq2_objs~{{{{predictor}}}}/mRFU/{cdp_paramspace_wildcard_pattern}/{{rfu_size_paramspace_instance_patterns}}.rds".format(
                    cdp_paramspace_wildcard_pattern=_escape_curly_braces(
                        cdp_paramspace.wildcard_pattern
                    )
                )
            ),
            glm_stats_tsv=_adjust_rfu_counts_for_covariates_batch_files(
                "{{{{prefix_no_drop_covariates}}}}/rfu_filtering_results/{{{{stats_type}}}}~{{{{predictor}}}}/mRFU/{cdp_paramspace_wildcard_pattern}/{{rfu_filtering_paramspace_instance_patterns}}.tsv.gz".format(
                    cdp_paramspace_wildcard_pattern=_escape_curly_braces(
                        cdp_paramspace.wildcard_pattern
                    )
                )
            ),
        benchmark:
            f"{{prefix_no_drop_covariates}}/adjusted_rfu_counts/{{stats_type}}~{{predictor}}/mRFU/glm_adjustment~log_lr/{cdp_paramspace.wildcard_pattern}.tsv.gz.benchmark"
        params:
            glm_adjustment="log_lr",
            fdr_cutoff=config["fdr_cutoff"],
            depth_colnames=config["glm"]["depth_colnames"],
            depth_poly_order=config["glm"]["depth_poly_order"],
        # According to https://github.com/serumdetect/snakemake-resource-usage/blob/339dfbc5eb336a9c8ddab1796a33f685a172e31e/notebooks/analyze-input-sizes-vs-memory-use.2024-05-28.py.ipynb,
        # this rule uses around 1000x of input file size (sometimes more). Give
        # it a bit more memory since we have to iterate through all samples.
        resources:
            mem_mb=lambda _, input, attempt: 1500
            * max(size_mb(fname) for fname in input.rfu_counts_tsv)
            * (0.5 + attempt / 2),
        conda:
            "../envs/feature_generation/adjust_rfu_counts.yaml"
        script:
            "../scripts/feature_generation/adjust_rfu_counts.R"

    ruleorder: adjust_rfu_counts_for_covariates_batch > adjust_rfu_counts_for_covariates

    # Compute adjusted RFU counts some covariates removed from the GLM model - batch version.
    use rule adjust_rfu_counts_for_covariates_batch as adjust_rfu_counts_for_covariates_batch_drop_covariates with:
        output:
            expand(
                "{{prefix_no_drop_covariates}}/drop_covariates~{{drop_covariates}}/adjusted_rfu_counts/{{stats_type}}~{{predictor}}/mRFU/glm_adjustment~log_lr/{cdp_paramspace_wildcard_pattern}/{rfu_filtering_paramspace_instance_patterns}.tsv.gz",
                cdp_paramspace_wildcard_pattern=[
                    cdp_paramspace.wildcard_pattern
                ],
                rfu_filtering_paramspace_instance_patterns=rfu_filtering_paramspace.instance_patterns,
            ),
        input:
            rfu_counts_tsv=lambda wildcards: _adjust_rfu_counts_for_covariates_batch_files(
                "{{{{prefix_no_drop_covariates}}}}/rfu_counts/mRFU/{cdp_paramspace_wildcard_pattern}/{{rfu_size_paramspace_instance_patterns}}.tsv.gz".format(
                    cdp_paramspace_wildcard_pattern=_escape_curly_braces(
                        cdp_paramspace.wildcard_pattern
                    )
                )
            ),
            size_factors_tsv="{prefix_no_drop_covariates}/size_factors/mRFU.tsv",
            metadata_tsv="{prefix_no_drop_covariates}/parsed_metadata/mRFU.metadata.tsv",
            deseq2_obj=_adjust_rfu_counts_for_covariates_batch_files(
                "{{{{prefix_no_drop_covariates}}}}/drop_covariates~{{{{drop_covariates}}}}/deseq2_objs~{{{{predictor}}}}/mRFU/{cdp_paramspace_wildcard_pattern}/{{rfu_size_paramspace_instance_patterns}}.rds".format(
                    cdp_paramspace_wildcard_pattern=_escape_curly_braces(
                        cdp_paramspace.wildcard_pattern
                    )
                )
            ),
            glm_stats_tsv=_adjust_rfu_counts_for_covariates_batch_files(
                "{{{{prefix_no_drop_covariates}}}}/drop_covariates~{{{{drop_covariates}}}}/rfu_filtering_results/{{{{stats_type}}}}~{{{{predictor}}}}/mRFU/{cdp_paramspace_wildcard_pattern}/{{rfu_filtering_paramspace_instance_patterns}}.tsv.gz".format(
                    cdp_paramspace_wildcard_pattern=_escape_curly_braces(
                        cdp_paramspace.wildcard_pattern
                    )
                )
            ),
        benchmark:
            f"{{prefix_no_drop_covariates}}/drop_covariates~{{drop_covariates}}/adjusted_rfu_counts/{{stats_type}}~{{predictor}}/mRFU/glm_adjustment~log_lr/{cdp_paramspace.wildcard_pattern}.tsv.gz.benchmark"

    ruleorder: adjust_rfu_counts_for_covariates_batch_drop_covariates > adjust_rfu_counts_for_covariates_drop_covariates


rule filter_rfus_for_blinded_covariates:
    """
    From the filtered RFUs, filter out those that are associated with blinded
    covariates in the main GLM model.
    """
    localrule: True
    output:
        output_tsv=f"{{prefix_no_drop_covariates}}/drop_covariates~{{fdr_cutoff}}~{{drop_covariates}}/adjusted_rfu_counts/{{stats_type}}~{{predictor}}/mRFU/glm_adjustment~{{glm_adjustment}}/{cdp_paramspace.wildcard_pattern}/{rfu_filtering_paramspace.wildcard_pattern}.tsv.gz",
    input:
        input_tsv=f"{{prefix_no_drop_covariates}}/drop_covariates~{{drop_covariates}}/adjusted_rfu_counts/{{stats_type}}~{{predictor}}/mRFU/glm_adjustment~{{glm_adjustment}}/{cdp_paramspace.wildcard_pattern}/{rfu_filtering_paramspace.wildcard_pattern}.tsv.gz",
        full_model_rfu_stats_tsv=f"{{prefix_no_drop_covariates}}/rfu_filtering_results/{{stats_type}}~{{predictor}}/mRFU/{cdp_paramspace.wildcard_pattern}/{rfu_filtering_paramspace.wildcard_pattern}.tsv.gz",
    conda:
        "../envs/feature_generation/filter_rfus_for_blinded_covariates.yaml"
    params:
        fdr_cutoff=lambda wildcards: float(wildcards.fdr_cutoff),
        blinded_covariates=lambda wildcards: wildcards.drop_covariates.split(
            "__"
        ),
    script:
        "../scripts/feature_generation/filter_rfus_for_blinded_covariates.py"


ruleorder: filter_rfus_for_blinded_covariates > adjust_rfu_counts_for_covariates_drop_covariates > adjust_rfu_counts_for_covariates


if config["batch_jobs"]:

    ruleorder: filter_rfus_for_blinded_covariates > adjust_rfu_counts_for_covariates_batch_drop_covariates


# Same but for GLM stats.
use rule filter_rfus_for_blinded_covariates as filter_glm_stats_for_blinded_covariates with:
    localrule: True
    output:
        output_tsv=f"{{prefix_no_drop_covariates}}/drop_covariates~{{fdr_cutoff}}~{{drop_covariates}}/rfu_filtering_results/{{stats_type}}~{{predictor}}/mRFU/{cdp_paramspace.wildcard_pattern}/{rfu_filtering_paramspace.wildcard_pattern}.tsv.gz",
    input:
        input_tsv=f"{{prefix_no_drop_covariates}}/drop_covariates~{{drop_covariates}}/rfu_filtering_results/{{stats_type}}~{{predictor}}/mRFU/{cdp_paramspace.wildcard_pattern}/{rfu_filtering_paramspace.wildcard_pattern}.tsv.gz",
        full_model_rfu_stats_tsv=f"{{prefix_no_drop_covariates}}/rfu_filtering_results/{{stats_type}}~{{predictor}}/mRFU/{cdp_paramspace.wildcard_pattern}/{rfu_filtering_paramspace.wildcard_pattern}.tsv.gz",


ruleorder: filter_glm_stats_for_blinded_covariates > filter_rfus_by_expansion


def ml_features_input(wildcards):
    """Get ML features input path depending on GLM stats."""
    if wildcards.glm_adjustment == "0":
        return expand(
            "{{prefix}}/norm_rfu_counts/mRFU/{cdp_params}"
            + "/{rfu_size_wildcard_pattern}.tsv.gz",
            cdp_params=list(cdp_paramspace.instance_patterns),
            rfu_size_wildcard_pattern=rfu_size_paramspace.wildcard_pattern,
        )
    else:
        return expand(
            "{{prefix}}/adjusted_rfu_counts/{{stats_type}}~{{predictor}}/mRFU"
            + "/glm_adjustment~{{glm_adjustment}}/{cdp_params}"
            + "/{rfu_filtering_wildcard_pattern}.tsv.gz",
            cdp_params=list(cdp_paramspace.instance_patterns),
            rfu_filtering_wildcard_pattern=rfu_filtering_paramspace.wildcard_pattern,
        )


rule prepare_ml_data:
    """Compile RFU data for ML using all samples."""
    output:
        f"{{prefix}}/adjusted_rfu_counts/{{stats_type}}~{{predictor}}/mRFU/glm_adjustment~{{glm_adjustment}}/{rfu_filtering_paramspace.wildcard_pattern}/ml_features.tsv.gz",
    input:
        rfu_counts=ml_features_input,
        glm_stats_tsv=expand(
            "{{prefix}}/rfu_filtering_results/{{stats_type}}~{{predictor}}/mRFU/{params}/{rfu_filtering_wildcard_pattern}.tsv.gz",
            params=list(cdp_paramspace.instance_patterns),
            rfu_filtering_wildcard_pattern=rfu_filtering_paramspace.wildcard_pattern,
        ),
    benchmark:
        f"{{prefix}}/adjusted_rfu_counts/{{stats_type}}~{{predictor}}/mRFU/glm_adjustment~{{glm_adjustment}}/{rfu_filtering_paramspace.wildcard_pattern}/ml_features.tsv.gz.benchmark"
    params:
        params=list(cdp_paramspace.instance_patterns),
        fdr_cutoff=config["fdr_cutoff"],
        positive_association_only=False,
    conda:
        "../envs/feature_generation/prepare_ml_data.yaml"
    script:
        "../scripts/feature_generation/prepare_ml_data.py"

"""
Snakemake rules for filtering RFUs through various cutoffs before performing
multiple testing correction.

There are three parameters to adjust (as listed in `rfu_filtering_paramspace`
below):

1. TCR size to be used for computing the "expression" of an RFU
2. Minimum RFU "expression" for it to be considered expanded in a given
   individual
3. Minimum number of individuals in which an RFU is expanded for it to be
   included for testing

While TCR size filtering has to be performed prior to GLM, the latter two
filters in combination restricted the RFUs to be tested. This can be done post-
hoc after GLM P values have been computed.
"""

import re


rule filter_rfus_by_expansion:
    """
    Count the number of significant RFUs by performing FDR restricted to
    RFUs of a certain level of expansion.
    """
    localrule: True
    group:
        "rfu_stats" if config["batch_jobs"] else None
    output:
        n_signif_hits=f"{{prefix}}/rfu_filtering_results/{{stats_type}}~{{predictor}}/mRFU/{cdp_paramspace.wildcard_pattern}/{rfu_filtering_paramspace.wildcard_pattern}.num_signif_filtered_rfus.tsv",
        filtered_rfus=f"{{prefix}}/rfu_filtering_results/{{stats_type}}~{{predictor}}/mRFU/{cdp_paramspace.wildcard_pattern}/{rfu_filtering_paramspace.wildcard_pattern}.tsv.gz",
    input:
        rfu_stats=f"{{prefix}}/{{stats_type}}~{{predictor}}/mRFU/{cdp_paramspace.wildcard_pattern}/{rfu_size_paramspace.wildcard_pattern}.tsv.gz",
        rfu_counts=lambda wildcards:
            f"{re.sub('/drop_covariates~[^/]+', '', wildcards.prefix)}/norm_rfu_counts/mRFU/{cdp_paramspace.wildcard_pattern}/{rfu_size_paramspace.wildcard_pattern}.tsv.gz",
    benchmark:
        f"{{prefix}}/rfu_filtering_results/{{stats_type}}~{{predictor}}/mRFU/{cdp_paramspace.wildcard_pattern}/{rfu_filtering_paramspace.wildcard_pattern}.tsv.gz.benchmark"
    params:
        min_expansion=lambda wildcards: float(wildcards.min_expansion),
        min_expanded_indivs=lambda wildcards: int(wildcards.min_expanded_indivs),
        fdr_cutoff=config["fdr_cutoff"],
    conda:
        "../envs/rfu_filtering/filter_rfus_by_expansion.yaml"
    resources:
        mem_mb=lambda _, attempt: 1000 + 2000 * attempt,
    script:
        "../scripts/rfu_filtering/filter_rfus_by_expansion.py"


# For slow d_cs, send the job to the cluster.
use rule filter_rfus_by_expansion as filter_rfus_by_expansion_cluster with:
    localrule: False
    group:
        "rfu_stats"
    output:
        n_signif_hits=f"{{prefix}}/rfu_filtering_results/{{stats_type}}~{{predictor}}/mRFU/{cdp_paramspace.wildcard_pattern}/{rfu_filtering_paramspace.wildcard_pattern}.num_signif_filtered_rfus.tsv",
        filtered_rfus=f"{{prefix}}/rfu_filtering_results/{{stats_type}}~{{predictor}}/mRFU/{cdp_paramspace.wildcard_pattern}/{rfu_filtering_paramspace.wildcard_pattern}.tsv.gz",
    input:
        rfu_stats=f"{{prefix}}/{{stats_type}}~{{predictor}}/mRFU/{cdp_paramspace.wildcard_pattern}/{rfu_size_paramspace.wildcard_pattern}.tsv.gz",
        rfu_counts=lambda wildcards:
            f"{re.sub('drop_covariates~[^/]+/', '', wildcards.prefix)}/norm_rfu_counts/mRFU/{cdp_paramspace.wildcard_pattern}/{rfu_size_paramspace.wildcard_pattern}.tsv.gz",
    benchmark:
        f"{{prefix}}/rfu_filtering_results/{{stats_type}}~{{predictor}}/mRFU/{cdp_paramspace.wildcard_pattern}/{rfu_filtering_paramspace.wildcard_pattern}.tsv.gz.benchmark"
    wildcard_constraints:
        # Specify the large d_c settings.
        d_c="|".join(
            snakemake.utils.Paramspace(
                pd.read_csv("config/cdp_paramspace.csv", dtype={"d_c": str})
                .query("is_slow_d_c")
                .drop(columns="is_slow_d_c")
            ).instance_patterns
        ),
    resources:
        mem_mb=lambda _, attempt: 1000 + 2000 * attempt,


ruleorder: filter_rfus_by_expansion_cluster > filter_rfus_by_expansion

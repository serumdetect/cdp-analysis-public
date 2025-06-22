"""Snakemake rules for annotating TCRs and RFUs."""

import snakemake.remote.S3


rule annotate_tcrs:
    """Annotate all TCRs with CDP data."""
    group:
        "annotate_tcrs"
    output:
        f"{{prefix_no_drop_covariates}}/annotated_tcrs/mRFU/{cdp_paramspace.wildcard_pattern}.ftr",
    input:
        tcr_data_ftr="{prefix_no_drop_covariates}/raw/mRFU.ftr",
        cluster_data_ftr=f"{{prefix_no_drop_covariates}}/decision_graph/mRFU/{cdp_paramspace.wildcard_pattern}.ftr",
    benchmark:
        f"{{prefix_no_drop_covariates}}/annotated_tcrs/mRFU/{cdp_paramspace.wildcard_pattern}.ftr.benchmark"
    # According to https://github.com/serumdetect/snakemake-resource-usage/blob/f77fdcdf983e21bce2c97ede4dadeb4a34718e19/notebooks/analyze-input-sizes-vs-memory-use.2024-05-28.py.ipynb,
    # this rule uses roughly 60 * input file size.
    resources:
        mem_mb=lambda _, input, attempt: size_mb(input.cluster_data_ftr)
        * 20
        * (0.5 + attempt / 2),
    conda:
        "../envs/annotate_tcrs_and_rfus/annotate_tcrs.yaml"
    script:
        "../scripts/annotate_tcrs_and_rfus/annotate_tcrs.py"


rule group_tcrs_by_aa_seq:
    """Group TCRs by amino acid sequence without RFU formation.

    This uses a special parameter setting of `normalize_dist~0` and `d_c~0`,
    which corresponds to RFU clustering requiring perfect amino acid match.
    """
    output:
        f"{{prefix_no_drop_covariates}}/annotated_tcrs/mRFU/normalize_dist~0/d_c~0.ftr",
    input:
        "{prefix_no_drop_covariates}/raw/mRFU.ftr",
    benchmark:
        f"{{prefix_no_drop_covariates}}/annotated_tcrs/mRFU/normalize_dist~0/d_c~0.ftr.benchmark"
    resources:
        mem_mb=lambda _, input, attempt: size_mb(input[0])
        * 20
        * (0.5 + attempt / 2),
    conda:
        "../envs/annotate_tcrs_and_rfus/annotate_tcrs.yaml"
    script:
        "../scripts/annotate_tcrs_and_rfus/group_tcrs_by_aa_seq.py"


ruleorder: group_tcrs_by_aa_seq > annotate_tcrs


rule compute_cluster_stats:
    """Compute basic cluster stats."""
    group:
        "annotate_tcrs"
    output:
        f"{{prefix_no_drop_covariates}}/cluster_stats~{{predictor}}/mRFU/{cdp_paramspace.wildcard_pattern}.tsv.gz",
    input:
        annotated_tcrs=f"{{prefix_no_drop_covariates}}/annotated_tcrs/mRFU/{cdp_paramspace.wildcard_pattern}.ftr",
        metadata_tsv="{prefix_no_drop_covariates}/parsed_metadata/mRFU.metadata.tsv",
    benchmark:
        f"{{prefix_no_drop_covariates}}/cluster_stats~{{predictor}}/mRFU/{cdp_paramspace.wildcard_pattern}.tsv.gz.benchmark"
    # According to https://github.com/serumdetect/snakemake-resource-usage/blob/f77fdcdf983e21bce2c97ede4dadeb4a34718e19/notebooks/analyze-input-sizes-vs-memory-use.2024-05-28.py.ipynb,
    # this rule uses 6 times input file size in MB.
    resources:
        mem_mb=lambda _, input, attempt: size_mb(input.annotated_tcrs)
        * 1.5
        * (0.5 + attempt / 2),
    conda:
        "../envs/annotate_tcrs_and_rfus/annotate_tcrs.yaml"
    script:
        "../scripts/annotate_tcrs_and_rfus/cluster_stats.py"

"""Wildcard constraints of the project."""

top_num_clonotypes_pat = "top_num_clonotypes_({})~\d+".format(
    "|".join(x for x in ["raw_ranks", "norm_templates", "norm_ranks"])
)
tcr_filtering_pattern = rf"(min_templates~\d+|{top_num_clonotypes_pat}|min_fold_of_(v_)?(quantile|mean)~\d+(\.\d+)?)"
cv_infix = r"(/rs~\d+/cv~\d+)?"
drop_covariates = r"[^~/]+"
fdr_cutoff_pat = r"\d+(\.\d+)?"
prefix_pattern_no_drop_covariates = f"(/{tcr_filtering_pattern})?" + cv_infix
prefix_pattern = (
    prefix_pattern_no_drop_covariates
    + f"(/drop_covariates(~{fdr_cutoff_pat})?~{drop_covariates})?"
)

# Prepare all splits, which are "data", "held_out_data" and all the possible
# keys in config["test_data"].
external_data_splits = "({})".format("|".join(list(config["test_data"])))
data_or_external_data_splits = "({})".format(
    "|".join(["data"] + list(config["test_data"]))
)
test_data_splits = "({})".format(
    "|".join(["held_out_data"] + list(config["test_data"]))
)
data_or_test_data_splits = "({})".format(
    "|".join(["data"] + list(config["test_data"]))
)
all_data_splits = "({})".format(
    "|".join(["data", "held_out_data"] + list(config["test_data"]))
)


wildcard_constraints:
    # Wildcards for TCR filtering.
    tcr_filtering=tcr_filtering_pattern,
    min_templates="\d+",
    top_n="\d+",
    min_fold_of_stat="\d+(\.\d+)?",
    #
    # The `normalize_dist` and `d_c` parameters are used for selecting the d_c
    # parameter of each output file.
    normalize_dist="[^~/]+",
    d_c="[^~/]+",
    #
    # Wildcards `min_dist`, `min_expansion` and `min_expanded_indivs` are the
    # RFU filtering wildards.
    min_size="\d+",
    min_expansion="\d+",
    min_expanded_indivs="\d+",
    #
    # If the `permutation` wildcard is present in the output file, it indicates
    # that the file is generated original from permuted metadata.
    permutation="\d+",
    #
    # Wildcard `prefix` starts with "data" followed by two optional components
    # for random state and CV fold.
    prefix="data" + prefix_pattern,
    prefix_no_drop_covariates="data" + prefix_pattern_no_drop_covariates,
    any_prefix=all_data_splits + prefix_pattern_no_drop_covariates,
    any_filtered_prefix=all_data_splits + "/" + tcr_filtering_pattern + cv_infix,
    any_unfiltered_prefix=all_data_splits + cv_infix,
    external_data=external_data_splits,
    external_data_prefix=external_data_splits
    + prefix_pattern_no_drop_covariates,
    test_data_prefix=test_data_splits + prefix_pattern_no_drop_covariates,
    rs="\d+",
    cv="\d+",
    cv_infix=cv_infix,
    drop_covariates=drop_covariates,
    fdr_cutoff=fdr_cutoff_pat,
    #
    # Wildcard `stats_type` indicates the type of GLM stats used to generate
    # the RFU statistical tests (including those used for GLM adjustment).
    # The values are:
    # - "glm_stats": GLM using DESeq2
    # - "rra_glm_stats": Robust rank aggregation using RFU counts adjusted with
    #   the GLM model for significant non-cancer status covariates.
    # - "robmixglm_stats": robust GLM using Robmixglm
    # - "rra_robmixglm_stats": Robust rank aggregation using RFU counts adjusted with
    #   the Robmixglm model for significant non-cancer status covariates.
    stats_type="(glm|rra_glm|rra_robmixglm|robmixglm|wilcox|fisher|bm)_stats",
    #
    # Wildcard `glm_adjustment` indicates whether the RFU counts were adjusted
    # for non-cancer status covariates based on the fitted GLM model. The values
    # are:
    # - 0: no adjustment
    # - 1: with adjustment
    # - log_lr: use log-likelihoods for cancer-not cancer states.
    glm_adjustment="0|1|log_lr",

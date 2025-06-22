"""
Snakemake script for computing Wilcoxon rank-sum test on normalized RFU counts.
"""

import functools

import pandas as pd
import scipy.stats
import sklearn.metrics
import statsmodels.stats.multitest
from tqdm.auto import tqdm

tqdm.pandas()


def compute_wilcox_stats(rfu, predictor):
    """
    Compute Mann-Whitney U rank-sum test for a given RFU and predictor.
    """
    mwu_res = scipy.stats.mannwhitneyu(
        rfu.loc[predictor], rfu.loc[~predictor], method="asymptotic"
    )
    auc = sklearn.metrics.roc_auc_score(predictor, rfu)
    return pd.Series({"auc": auc, "pvalue": mwu_res.pvalue})


def compute_fisher_stats(rfu, predictor):
    """
    Compute Fisher's exact test for a given RFU and predictor.
    """
    contingency_table = pd.crosstab(
        (rfu > 0).astype(int), predictor.astype(int)
    )
    categories = [1, 0]
    for category in categories:
        if category not in contingency_table.index:
            contingency_table.loc[category] = 0
        if category not in contingency_table.columns:
            contingency_table[category] = 0
    contingency_table = contingency_table.loc[categories, categories]
    fisher_res = scipy.stats.fisher_exact(contingency_table)
    return pd.Series({"odds_ratio": fisher_res[0], "pvalue": fisher_res[1]})


def compute_brunner_munzel_stats(rfu, predictor):
    """
    Compute Brunner-Munzel test for a given RFU and predictor.
    """
    bm_res = scipy.stats.brunnermunzel(rfu.loc[predictor], rfu.loc[~predictor])
    auc = sklearn.metrics.roc_auc_score(predictor, rfu)
    return pd.Series({"auc": auc, "pvalue": bm_res.pvalue})


def main(norm_rfu_counts_tsv, metadata, predictor, stats_type, output_tsv):
    # Read inputs.
    metadata = (
        pd.read_table(metadata, index_col="sample_name")
        .dropna(subset=predictor)
        .astype({predictor: bool})
    )
    norm_rfu_counts = pd.read_table(norm_rfu_counts_tsv, index_col=0)[
        metadata.index
    ]

    # Compute Wilcoxon rank-sum test for each RFU.
    predictor_vec = metadata[predictor]
    try:
        stats_fn = {
            "wilcox": compute_wilcox_stats,
            "fisher": compute_fisher_stats,
            "bm": compute_brunner_munzel_stats,
        }[stats_type]
    except KeyError:
        raise ValueError(f"Unknown stats type: {stats_type}")
    stats_df = norm_rfu_counts.progress_apply(
        functools.partial(stats_fn, predictor=predictor_vec), axis=1
    ).rename(columns=lambda col: f"{predictor}.{col}")

    # Create the output table and save it.
    stats_df[f"{predictor}.fdr"] = statsmodels.stats.multitest.fdrcorrection(
        stats_df[f"{predictor}.pvalue"]
    )[1]
    stats_df.sort_values(f"{predictor}.fdr").to_csv(
        output_tsv, sep="\t", index=True
    )


if __name__ == "__main__":
    main(
        norm_rfu_counts_tsv=snakemake.input.norm_rfu_counts_tsv,
        metadata=snakemake.input.metadata,
        predictor=snakemake.wildcards.predictor,
        stats_type=snakemake.wildcards.stats_type,
        output_tsv=snakemake.output[0],
    )

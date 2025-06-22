"""Snakemake script for splitting the TCR data into training and test sets."""

import loguru
import pandas as pd
import sklearn.model_selection


def validate_cv_params(num_cv_folds, cv_fold):
    """Validate the cross-validation parameters.

    One of the following must hold.
    - Both num_cv_folds and cv_fold are None.
    - `num_cv_folds` is an integer and `cv_fold` is an integer
    """
    if num_cv_folds is None and cv_fold is None:
        return
    else:
        assert isinstance(num_cv_folds, int)
        assert isinstance(cv_fold, int)


def cross_validate_tcr_data(
    num_cv_folds, cv_fold, tcrs, metadata, stratifier_colname, random_state
):
    """Stratified splitting of the TCR data into training and test sets."""
    # For the sake of reproducibility, we need to set the random state.
    if random_state is None:
        loguru.logger.warning("Setting random state to 0 for reproducibility.")
        random_state = 0

    if stratifier_colname is not None:
        splitter = sklearn.model_selection.StratifiedKFold(
            n_splits=num_cv_folds, shuffle=True, random_state=random_state
        )
        splits = list(
            splitter.split(
                X=metadata[[stratifier_colname]],
                y=metadata[stratifier_colname],
            )
        )
        train_idx, test_idx = splits[cv_fold]
    else:
        splitter = sklearn.model_selection.KFold(
            n_splits=num_cv_folds, shuffle=True, random_state=random_state
        )
        splits = list(splitter.split(X=metadata))
        train_idx, test_idx = splits[cv_fold]

    # Extract the test samples.
    metadata_test = metadata.iloc[test_idx]
    tcrs_test = tcrs.loc[
        tcrs["sample_name"].isin(metadata_test["sample_name"])
    ].reset_index(drop=True)
    if tcrs_test["sample_name"].dtype.name == "category":
        tcrs_test["sample_name"] = tcrs_test[
            "sample_name"
        ].cat.remove_unused_categories()

    # Extract the train samples. Overwrite the original data to save memory.
    metadata = metadata.iloc[train_idx]
    tcrs = tcrs.loc[
        tcrs["sample_name"].isin(metadata["sample_name"])
    ].reset_index(drop=True)
    if tcrs_test["sample_name"].dtype.name == "category":
        tcrs["sample_name"] = tcrs["sample_name"].cat.remove_unused_categories()

    return tcrs, metadata, tcrs_test, metadata_test


def main(
    tcrs_ftr,
    metadata_tsv,
    cv_grouping_colname,
    num_cv_folds,
    random_state,
    cv_fold,
    output_ftr,
    output_metadata_tsv,
    output_held_out_ftr,
    output_metadata_held_out_tsv,
):
    # Load inputs.
    tcrs = pd.read_feather(tcrs_ftr)
    metadata = pd.read_csv(metadata_tsv, sep="\t")

    # Issue 36: split to exclude the current test samples.
    validate_cv_params(num_cv_folds, cv_fold)
    tcrs, metadata, tcrs_test, metadata_test = cross_validate_tcr_data(
        num_cv_folds,
        cv_fold,
        tcrs,
        metadata,
        cv_grouping_colname,
        random_state,
    )

    tcrs_test.to_feather(output_held_out_ftr)
    metadata_test.to_csv(
        output_metadata_held_out_tsv, header=True, sep="\t", index=False
    )
    tcrs.to_feather(output_ftr)
    metadata.to_csv(output_metadata_tsv, header=True, sep="\t", index=False)


if __name__ == "__main__":
    main(
        tcrs_ftr=snakemake.input.data,
        metadata_tsv=snakemake.input.metadata,
        cv_grouping_colname=snakemake.params.cv_grouping_colname,
        num_cv_folds=snakemake.params.num_cv_folds,
        random_state=snakemake.params.rs,
        cv_fold=snakemake.params.cv,
        output_ftr=snakemake.output.data,
        output_metadata_tsv=snakemake.output.metadata,
        output_held_out_ftr=snakemake.output.held_out_data,
        output_metadata_held_out_tsv=snakemake.output.held_out_metadata,
    )

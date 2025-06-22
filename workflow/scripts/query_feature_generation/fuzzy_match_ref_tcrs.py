"""Snakemake script for matching query TCRs with reference TCRs RFUs.

Matching is done by searching the TCR index and finding the nearest reference
TCR.
"""

import joblib
import pandas as pd
import tcrnn


def monkey_patch_tcr_repertoire(tcr_repertoire, n_jobs):
    """Monkey patch a TCR repertoire object to its newest class version."""
    new_tcr_repertoire = tcrnn.TcrRepertoire(
        data=tcr_repertoire.data,
        tcr_metric=tcr_repertoire.tcr_metric,
        n_jobs=n_jobs,
        progress=tcr_repertoire.progress,
    )
    new_tcr_repertoire._dedup_data = tcr_repertoire._dedup_data
    new_tcr_repertoire._encoded_tcrs = tcr_repertoire._encoded_tcrs
    new_tcr_repertoire.nn_index = tcr_repertoire.nn_index
    return new_tcr_repertoire


def main(
    query_ftr,
    ref_repertoire_path,
    monkey_patch_ref,
    n_jobs,
    output_ftr,
):
    # Load the data.
    query_tcrs = pd.read_feather(query_ftr)
    ref_repertoire = joblib.load(ref_repertoire_path)

    if monkey_patch_ref:
        # If requested, monkey patch the reference TCR repertoire to the current
        # code version for searching for the nearest reference TCR without
        # having to repeat the indexing.
        ref_repertoire = monkey_patch_tcr_repertoire(
            ref_repertoire, n_jobs
        )

    # Find the nearest reference TCRs.
    query_tcrs = tcrnn.TcrRepertoire(
        query_tcrs,
        ref_repertoire.tcr_metric,
        n_jobs=n_jobs,
    )
    nearest_ref_tcrs, delta = ref_repertoire.nearest_tcr(query_tcrs)
    nearest_ref_tcrs = nearest_ref_tcrs.to_frame(index=False)
    nearest_ref_tcrs.columns = ["ref_v_gene", "ref_cdr3"]
    nearest_ref_tcrs["delta"] = delta
    nearest_ref_tcrs = pd.concat(
        [query_tcrs.dedup_data.index.to_frame(index=False), nearest_ref_tcrs],
        axis=1,
    )

    # Output data.
    nearest_ref_tcrs.to_feather(output_ftr)


if __name__ == "__main__":
    main(
        query_ftr=snakemake.input.query_ftr,
        ref_repertoire_path=snakemake.input.ref_repertoire_path,
        monkey_patch_ref=snakemake.params.monkey_patch_ref_repertoire,
        n_jobs=snakemake.threads,
        output_ftr=snakemake.output.query_nearest_ref_tcrs_ftr,
    )

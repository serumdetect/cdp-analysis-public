"""Snakemake script for computing CDP decision graph."""

import joblib
import numpy as np
import pandas as pd
import tcrnn.cfsfdp


def main(
    tcr_repertoire_sav,
    densities_ftr,
    d_c,
    output_ftr,
):
    # Load the TCR repertoire index.
    tcr_index = joblib.load(tcr_repertoire_sav)

    # Load the densities. For backwards compatiblity (for the `serum_rfu`)
    # dataset, ensure that a column named "v_gene" exists.
    density_df = pd.read_feather(densities_ftr)
    if "v_gene" not in density_df:
        density_df["v_gene"] = np.nan
    density = density_df.set_index(["v_gene", "cdr3"])["density"]

    # Create TCR CDP object. Use cutoff of 5.0 for normalized distances and 50
    # for unnormalized distances.
    delta_cutoff = 5.0 if d_c <= 5.0 else 50
    tcr_cdp = tcrnn.cfsfdp.NnCfsfdp(
        tcr_index, d_c=d_c, density=density, delta_cutoff=delta_cutoff
    )

    # Compute delta and centroid.
    centroid = tcr_cdp.assign_centroid()

    # Save the CDP results.
    cdp_results = pd.DataFrame(
        {
            "idx": np.arange(len(centroid)),
            "size": tcr_cdp.size.values,
            "density": tcr_cdp.density,
            "delta": tcr_cdp.delta,
            "parent_idx": tcr_cdp.parent_idx,
            "centroid": centroid,
        }
    )

    cdp_results.reset_index().to_feather(output_ftr)


if __name__ == "__main__":
    main(
        tcr_repertoire_sav=snakemake.input.tcr_repertoire,
        densities_ftr=snakemake.input.density,
        d_c=float(snakemake.wildcards.d_c),
        output_ftr=snakemake.output[0],
    )

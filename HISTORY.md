# Changelog


## [Unreleased]

- Added support for using HLA alleles as standalone terms.
- Bug fix: remove test samples with warning if they have factor levels that are
  missing in train data. This is relevant for instance when test samples have VJ
  VJ lots that are not in the train data.


## [1.2.0]

- Update Dockerfile s3data version to 2.2.0.
- Implement filtering rule `filter_by_fold_v_dependent_quantile`.
- Update Conda environment files to include hidden dependencies.
- Add VJ lot as a default covariate for GLM.
- Update `plot_top_rfus` background color (to white) and color gradient.
- Bug fix: ensure there are no nans in `race_parsed` column.
- Update config.yaml to the settings that have been used for the past few
  months.


## [1.1.0]

- Bug fixes on both the Snakemake workflow, code and tests.
- Minor black formatting changes.
- Add TCR repertoire depth (both median and total UMI count) as a covariate to
  be corrected on in the GLM.
- Add helper function for computing S3 file sizes.
- Make main rules use dynamic memory request depending on input file size.
- Remove excluded `sample_name` categories from the TCR feather files after TCR
  filtering.
- Compute pre-filtering UMI and clonotype counts per sample.
- Add batched jobs for RFU counts adjustment rules, combining multiple d_c and
  RFU filtering settings into a single rule.
- Make compute_query_rfu_counts use a single thread instead. Polars did not seem
  to be using more than one thread in any case.
- Efficiency update for compute_query_rfu_counts: when TCR min_size = 1, do not
  bother with the complicated size-aware RFU counting algorithm. Instead, simply
  tally the total number of TCR clonotypes for each sample and RFU.
- Allow the use of external test data.
- Conform with s3data v2.
- Add functionality to run the CDP pipeline with specific train GLM covariates
  unavailable in the test data.
- Enable local Snakemake rules execution by adjusting size_mb function.
- Add rules for calculating and visualizing proportions of V(J)-genes among
  clonotypes/tcrs in sample repertoires.
- Add rules for augmenting metadata with pre- and post-filtering
  TCR-table-derived stats.
- Add rules to produce RFU stats tables which augment GLM stats tables with
  VJ-gene annotations.
- Set explicit package versions in conda-env.yaml to speed up the Conda
  enviroment creation process.


## [1.0.0]

- Update s3data version to v1.2.3 in Dockerfile.
- Add rule for annotating dataset TCRs with inhouse TILs.
- Major refactoring to include full CV capability that includes TCR CDP.


## 0.1.0

- Start tracking change log.

# Cluster TCRs into RFUs and compute RFU features for ML (with cross-validation)

This repo implements a workflow for clustering TCRs into RFUs, finding cancer-
associated RFUs, and generating RFU counts-based features for ML downstream.

The workflow can be run in a cross-validation setting where only the training
samples are used for TCR clustering and cancer-association discovery. In
addition, the training data can be used to generate test data TCR/RFU features
that are adjusted using the fitted GLM models.

## Installation

Create a new Conda environment with **Python 3.11** and install the repo
requirements.

```bash
# Install all public dependencies.
conda create -n <env-name> "python=3.11"
conda install --file conda-requirements.txt

# Load and install tcrnn.
git clone https://github.com/serumdetect/tcrnn-public
cd tcrnn-public
git submodule update --init --recursive
conda install --file conda-requirements.txt
pip install ./pynndescent
pip install .
```


## Usage:

### Inputs

Inputs are TCR repertoires in a `feather` format table and an associated
metadata table in TSV format. See [config/config.yaml](./config/config.yaml) for
details.

Other workflow parameters are listed and documented in the same configuration
file.

Request the workflow outputs using Snakemake as follows with the target files
described below.

```bash
snakemake -n -c 1 --use-conda -- <targets>
```


### Main results and Snakemake target files

The log-LR ML features can be generated using target `data/adjusted_rfu_counts/glm_stats~is_cancer/mRFU/glm_adjustment~log_lr/min_size~1/min_expansion~8/min_expanded_indivs~3/ml_features.tsv.gz`. These are log-likelihood values for the fitted GLM model with predictor value being `TRUE` vs `FALSE` (or 1 vs 0) given the observed RFU counts and repertoire depths.

This file will combine RFU models fitted across all the d_c settings specified in configuration file [config/cdp_paramspace.csv](./config/cdp_paramspace.csv). The three columns in this file have the following meaning.

Column name | Description
:-----------|:-----------
normalize_dist | Whether the CDR3 dissimilarity metric is normalized by the length of the longer CDR3 sequence (roughly `#mismatches` / `CDR3 length`)
d_c | Maximum dissimilarity for two TCRs to be grouped together. See the d_c parameter in the original CFSFDP algorithm
is_slow_d_c | Whether this `normalize_dist`/`d_c` setting will result in a lot of RFUs. When running with Kubernetes, some rules dynamically keep jobs from fast `d_c` settings local while sending the slow jobs to Kubernetes.

Parameter value `normalize_dist` = 0, `d_c` = 0 has a special meaning of bypassing RFU formation and grouping TCRs directly by V gene and CDR3 sequence identity instead. If [config/cdp_paramspace.csv](./config/cdp_paramspace.csv) only contains one row with these values, then no RFU formation will be performed at all.


### Statistical inference (DESeq2) results

The GLM statistics pre-filtering can be generated using target `data/glm_stats~{predictor}/mRFU/normalize_dist~{normalize_dist}/d_c~{d_c}/min_size~{min_size}.tsv.gz`. Post-TCR/RFU filtering (see below) FDRs can be generated using target `data/glm_stats~{predictor}/mRFU/min_size~{min_size}/min_expansion~{min_expansion}/rfu_stats.tsv.gz`. This file contains the P and FDR values for each TCR/RFU after they have been filtered by *size* and the level of *expansion*.

The target file `data/glm_coef_mats/glm_stats~{predictor}/mRFU/normalize_dist~{normalize_dist}/d_c~{d_c}/min_size~{min_size}.tsv.gz` contains the actual fitted GLM coefficients.

See rule `compute_deseq2_stats` for additional details on DESeq2 statistics outputs that are produced.


### Meaning of the RFU filtering parameters

`min_size` (deprecated), `min_expansion` and `min_expanded_indivs` are RFU/TCR level filters that allow excluding RFUs/TCRs from statistical testing to increase statistical power. The meaning of these parameters are as follows.


Parameter label | Parameter meaning
----------------|------------------
*min_size* | **Deprecated**. Only TCRs with the total number of nucleotide level distinct clonotypes across all samples in the *training* data greater than this are used to tally the RFU count. This filter is experimental and only relevant for RFUs.
*min_expansion* | Depth-normalized TCR/RFU count threshold to use.
*min_expanded_indivs* | Minimum number of samples with depth-normalized RFU count above `min_expansion` for the RFU to be included.


### TCR filtering

Raw TCRs can be filtered before they are entered to TCR/RFU counting. This can be used to exclude low read/UMI count TCRs that are more likely to correspond to naive TCRs or PCR errors. In order to request results with TCR filtering, a TCR filtering tag can be added right after the dataset tag in the path (`data` for the global dataset, `held_out_data` for CV test splits, etc.). Currently available TCR filtering schemes are:

- `data/min_templates~<int:min_tcr_templates>/`
- `data/min_fold_of_{quantile|mean}~<float:min_fold>/`  # Min fold of median (hard-coded quantile) or trimmed mean.
- `data/top_num_clonotypes~{raw_ranks|norm_templates|norm_ranks}/~<int:top_n>`

For instance, target `data/min_fold_of_mean~0.5/glm_stats~{predictor}/mRFU/normalize_dist~{normalize_dist}/d_c~{d_c}/min_size~{min_size}.tsv.gz` generates DESeq2 statistics using a TCR set that only includes TCRs with `templates` greater than the 0.5 * mean `templates` count of each sample.


### Cross-validation

To run cross-validation, set `cv.num_cv_folds` to e.g. 10 in [config.yaml](config/config.yaml). Then request the following target files. The only difference between the two paths below is the `data`/`held_out_data` prefix.

* `data/rs~{random_state}/cv~{cv_fold}/adjusted_rfu_counts/glm_stats~{predictor}/mRFU/glm_adjustment~log_lr/min_size~{min_size}/min_expansion~{min_expansion}/min_expanded_indivs~{min_expanded_indivs}/ml_features.tsv.gz"`
* `held_out_data/rs~{random_state}/cv~{cv_fold}/adjusted_rfu_counts/glm_stats~{predictor}/mRFU/glm_adjustment~log_lr/min_size~{min_size}/min_expansion~{min_expansion}/min_expanded_indivs~{min_expanded_indivs}/ml_features.tsv.gz"`

To generate ML features for test data based on the GLMs fitted on the training data, set the test datasets in [config.yaml](./config/config.yaml) and generate their targets by replacing the `data` in target paths with the corresponding test dataset key used in configuration `test_data`.

To use cross-validation with TCR filtering, the pattern is e.g. `data/{tcr_filtering}/rs~{random_state}/cv~{cv_fold}...`.


### Generating features with specific GLM covariates being blinded in test data

These results can be requested by appending to the prefix (e.g. `data/`) a pattern of `drop_covariates~<fdr>~<blinded_covariates>`. A GLM is fitted against all the covariates specified in the [config.yaml](./config/config.yaml) minus the covariates specified in `<blinded_covariates>`, which is a double-underscore (`__`) delimited list of blinded covariate names as they appear in the metadata columns, e.g. `race_parsed__age_years__Gender`. `<fdr>` is the FDR cutoff used to exclude covariate-dropped-GLM cancer-associated RFUs that are associated with FDR ≤ `<fdr>` in the full train data GLM.

For example, target `data/{tcr_filtering}/rs~{random_state}/cv~{cv_fold}/drop_covariates~0.1~age_years/adjusted_rfu_counts/glm_stats~is_cancer/mRFU/glm_adjustment~log_lr/min_size~{min_size}/min_expansion~{min_expansion}/min_expanded_indivs~{min_expanded_indivs}/ml_features.tsv.gz` would generate log-LR features TCRs/RFUs significant for metadata column `is_cancer` for which the significance for covariate `age_years` is less than FDR of 0.1


### Batch mode for rule `adjust_rfu_counts_for_covariates`

Jobs from rule `adjust_rfu_counts_for_covariates` are brief and repetitive in that the input files are identical given `min_size`. In a CV analysis where a large number of RFU filtering parameter combinations are searched for, the deployment of a large number of small Snakemake jobs on the Kubernetes cluster can cause significant overhead. By batching, all `adjust_rfu_counts_for_covariates` jobs from the same `d_c` to a Kubernetes single job (across specified RFU filtering parameters).

Batching only applies to RFU counts adjustment with `glm_adjustment~log_lr`.

To enable this, set configuration value `batch_jobs` in `config/config.yaml` to `true`. When enabled, Snakemake will read all RFU filtering parameter combinations from [config/rfu_filtering_params.csv](./config/rfu_filtering_params.csv) and run them in batch mode.

Note that the batch rule will always recompute the outputs for all parameter values in [config/rfu_filtering_params.csv](./config/rfu_filtering_params.csv). To compute `log_lr` features for additional RFU filtering parameters, one should remove already completed values from [config/rfu_filtering_params.csv](./config/rfu_filtering_params.csv) before adding new values to that file.


## Testing

To test the Python code, run `pytest tests/`.

To run the R code, run `testthat::test_dir("./tests")` from an R console.

---

THIS REPOSITORY IS COPYRIGHT OF SERUM DETECT, INC

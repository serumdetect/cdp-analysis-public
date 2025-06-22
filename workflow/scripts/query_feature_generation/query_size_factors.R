#' Snakemake script for computing query size factors.


library(scuttle)
library(pbapply)


#' Compute query size factors for one sample.
#'
#' Compute size factors for test samples that are brought to the scale of the
#' provided reference samples.
#'
#' The function assumes that the reference size factors were computed using
#' function `scuttle::pooledSizeFactors`.
#'
#' @param ref_mat Matrix (dimensions of m * n_ref) of reference samples.
#' @param ref_size_factors Vector (length n) of size factors for reference
#'     samples.
#' @param query_vec Vector (length m) of query samples counts.
#' @param ... Additional arguments to be passed to `scuttle::pooledSizeFactors`.
#'
#' @importFrom scuttle pooledSizeFactors
#'
#' @return Size factor for the query sample.
#'
#' @export
test_sample_size_factor = function(query_vec, ref_mat, ref_size_factors, ...) {
    # Sanity checks.
    if (ncol(ref_mat) == 0) {
        stop("`ref_mat`` must have at least one column.")
    }
    if (ncol(ref_mat) != length(ref_size_factors)) {
        msg = paste0(
            "Number of columns in `ref_mat` and length of `ref_size_factors`",
            " must be the same."
        )
        stop(msg)
    }
    if (nrow(ref_mat) != length(query_vec)) {
        msg = "Number of rows in `ref_mat` must equal length of `query_vec`."
        stop(msg)
    }

    # The reference RFU counts matrix must have column (sample) names and they
    # must match the names of the reference size factors if available.
    if (is.null(colnames(ref_mat))) {
        stop("Column names of `ref_mat` must be present.")
    }
    else if (is.null(names(ref_size_factors))) {
        names(ref_size_factors) = colnames(ref_mat)
    } else {
        if (
            is.null(colnames(ref_mat))
            || !all(colnames(ref_mat) == names(ref_size_factors))
        ) {
            msg = "Column names of `ref_mat` and `ref_size_factors` must match."
            stop(msg)
        }
    }

    # Combine the query sample with the reference samples.
    combined_mat = cbind(ref_mat, .query = query_vec)
    sf_all = scuttle::pooledSizeFactors(combined_mat, ...)
    sf_ratio_train = median(
        sf_all[-ncol(combined_mat)] / ref_size_factors[colnames(ref_mat)]
    )
    sf_query = sf_all[ncol(combined_mat)] / sf_ratio_train

    return(sf_query)
}


#' Compute test data size factors for a group of test samples individually.
#'
#' Compute size factors for test samples that are brought to the scale of the
#' provided reference samples.
#'
#' @param ref_mat Matrix (dimensions of m * n_ref) of reference samples.
#' @param ref_size_factors Vector (length n) of size factors for reference
#'     samples.
#' @param query_mat Matrix (dimensions of m * n_query) of query samples.
#' @param ... Additional arguments to be passed to `test_sample_size_factor`.
#'
#' @importFrom scuttle pooledSizeFactors
#' @importFrom pbapply pbapply
#'
#' @return A vector of size factors for the query samples.
#'
#' @export
test_sample_size_factors = function(
    query_mat, ref_mat, ref_size_factors, n_jobs = NULL, ...
) {
    # Sanity checks.
    if (ncol(query_mat) == 0) {
        stop("`query_mat` must have at least one column.")
    }
    if (nrow(ref_mat) != nrow(query_mat)) {
        msg = "Number of rows in `ref_mat` and `query_mat` must be the same."
        stop(msg)
    }

    # Sanity check: make sure there are no common samples between the reference
    # and query matrices.
    common_samples = intersect(colnames(ref_mat), colnames(query_mat))
    if (length(common_samples)) {
        msg = sprintf(
            paste0(
                "The following samples (N = %s) are present in both `ref_mat`",
                " and `query_mat`: %s."
            ),
            length(common_samples),
            paste(common_samples, collapse = ", ")
        )
        warning(msg)
    }

    # Concatenate the reference and query matrices and compute new size factors.
    sf_query = pbapply(
        query_mat,
        2,
        test_sample_size_factor,
        ref_mat = ref_mat,
        ref_size_factors = ref_size_factors,
        ...,
        cl = n_jobs
    )
    names(sf_query) = colnames(query_mat)

    return(sf_query)
}


main = function(
    query_rfu_counts_tsv,
    ref_rfu_counts_tsv,
    ref_size_factors_tsv,
    output_size_factors
) {
    # Load input RFU counts and make sure the RFUs match between the query and
    # reference matrices.
    query_rfu_counts = read.table(
        query_rfu_counts_tsv,
        header = TRUE,
        sep = "\t",
        row.names = 1,
        check.names = FALSE
    )
    ref_rfu_counts = read.table(
        ref_rfu_counts_tsv,
        header = TRUE,
        sep = "\t",
        row.names = 1,
        check.names = FALSE
    )
    query_rfu_counts = query_rfu_counts[rownames(ref_rfu_counts), ]
    rownames(query_rfu_counts) = rownames(ref_rfu_counts)
    query_rfu_counts[is.na(query_rfu_counts)] = 0

    # Load reference size factors and compute query size factors.
    temp = read.table(ref_size_factors_tsv, header = FALSE, sep = "\t")
    ref_size_factors = setNames(temp[, 2], temp[, 1])

    query_size_factors = test_sample_size_factors(
        query_rfu_counts, ref_rfu_counts, ref_size_factors
    )
    stopifnot(all(names(query_size_factors) == colnames(query_rfu_counts)))

    # Output the query size factors.
    write.table(
        data.frame(query_size_factors),
        output_size_factors,
        sep = "\t",
        quote = FALSE,
        row.names = TRUE,
        col.names = FALSE
    )
}


if (!interactive()) {
    main(
        query_rfu_counts_tsv = snakemake@input[["query_rfu_counts_tsv"]],
        ref_rfu_counts_tsv = snakemake@input[["ref_rfu_counts_tsv"]],
        ref_size_factors_tsv = snakemake@input[["ref_size_factors_tsv"]],
        output_size_factors = snakemake@output[[1]]
    )
}

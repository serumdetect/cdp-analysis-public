#' Snakemake script for adjusting query RFU counts based on the fitted reference
#' GLM model.

suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(doFuture))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(future))
suppressPackageStartupMessages(library(logger))
suppressPackageStartupMessages(library(progressr))
suppressPackageStartupMessages(library(SummarizedExperiment))


#' Align matrix rows to provided row names and fill missing values with 0.
#'
#' @param mat The input matrix.
#' @param rownames The desired row names.
#' @param fill_value The value to fill for missing values.
#'
#' @return The aligned matrix.
align_matrix_to_rows = function(mat, rownames, fill_value = 0) {
    aligned_mat = matrix(
        fill_value,
        nrow = length(rownames),
        ncol = ncol(mat),
        dimnames = list(rownames, colnames(mat))
    )
    shared_rownames = intersect(rownames, rownames(mat))
    aligned_mat[shared_rownames, ] = mat[shared_rownames, ]

    return(aligned_mat)
}


#' Extract the reference model matrix.
extract_ref_model_matrix = function(ref_deseq2_dds) {
    ref_col_data = extract_ref_col_data(ref_deseq2_dds)

    # "(Intercept)" must be renamed to "Intercept" for the model.matrix() call
    # to match the column name in the coefficients matrix.
    ref_model_matrix = model.matrix(
        design(ref_deseq2_dds), data = ref_col_data
    ) %>%
        as.data.frame() %>%
        dplyr::rename(Intercept = "(Intercept)") %>%
        rename_with(make.names) %>%
        as.matrix()

    return(ref_model_matrix)
}


#' Create a Cached File Reader
#'
#' Creates a cached file reader using closures to read and cache data from a
#' file.
#'
#' @return A list of functions for reading and accessing cached file data.
#' @export
create_cached_file_reader = function(
    reader_fn, post_processing_fn = NULL, ...
) {
    current_file = NULL
    cached_data = NULL

    return(
        function(file_path) {
            if (!identical(file_path, current_file)) {
                logger::log_info(sprintf("Read new file %s...", file_path))
                current_file <<- file_path
                cached_data <<- reader_fn(file_path, ...)
                if (!is.null(post_processing_fn)) {
                    cached_data <<- post_processing_fn(cached_data)
                }
            } else {
                logger::log_info(sprintf("Return cached file %s...", file_path))
            }
            return(cached_data)
        }
    )
}


#' Create an input loader, which consists of cached file readers for each input.
create_cached_input_reader = function() {
    rfu_counts_reader = create_cached_file_reader(
        reader_fn = read.table,
        post_processing_fn = as.matrix,
        sep = "\t",
        header = TRUE,
        row.names = 1,
        check.names = FALSE
    )
    size_factors_reader = create_cached_file_reader(
        reader_fn = read.table,
        function(df) setNames(df[[2]], df[[1]]),
        sep = "\t",
        header = FALSE
    )
    metadata_reader = create_cached_file_reader(
        reader_fn = read.table,
        header = TRUE,
        sep = "\t",
        row.names = "sample_name",
        comment = "",
        quote = '"'
    )
    ref_deseq2_dds_reader = create_cached_file_reader(reader_fn = readRDS)
    ref_glm_stats_reader = create_cached_file_reader(
        reader_fn = read.table,
        sep = "\t",
        header = TRUE,
        row.names = 1,
        check.names = FALSE
    )

    list(
        rfu_counts_reader = rfu_counts_reader,
        size_factors_reader = size_factors_reader,
        metadata_reader = metadata_reader,
        ref_deseq2_dds_reader = ref_deseq2_dds_reader,
        ref_glm_stats_reader = ref_glm_stats_reader
    )
}


#' Load all inputs.
load_inputs = function(
    cached_input_reader,
    rfu_counts_tsv,
    size_factors_tsv,
    metadata_tsv,
    ref_deseq2_dds,
    ref_glm_stats_tsv,
    predictor
) {
    # Load all the inputs with the cached input reader.
    ref_glm_stats = cached_input_reader$ref_glm_stats_reader(ref_glm_stats_tsv)
    rfu_counts = cached_input_reader$rfu_counts_reader(rfu_counts_tsv) %>%
        align_matrix_to_rows(rownames(ref_glm_stats), 0)
    metadata = cached_input_reader$metadata_reader(metadata_tsv)[
        colnames(rfu_counts),
    ]
    ref_deseq2_dds = cached_input_reader$ref_deseq2_dds_reader(ref_deseq2_dds)[
        rownames(ref_glm_stats), , drop = FALSE
    ]
    size_factors = cached_input_reader$size_factors_reader(size_factors_tsv)[
        colnames(rfu_counts)
    ]

    # Remove samples with missing metadata that is necessary for the GLM model.
    ref_covariates = setdiff(
        colnames(colData(ref_deseq2_dds)),
        c("sizeFactor", predictor)
    )
    complete_cases = complete.cases(metadata[, ref_covariates, drop = FALSE])
    if (any(!complete_cases)) {
        msg = sprintf(
            paste0(
                "Removing %s samples with missing metadata necessary for the",
                " model."
            ),
            sum(!complete_cases)
        )
        logger::log_warn(msg)
        metadata = metadata[complete_cases, , drop = FALSE]
        rfu_counts = rfu_counts[, complete_cases, drop = FALSE]
        size_factors = size_factors[complete_cases]
    }

    logger::log_info("Loading reference matrix...")
    ref_model_matrix = extract_ref_model_matrix(ref_deseq2_dds)
    beta_mat = coef(ref_deseq2_dds)

    # Load the query model matrix.
    logger::log_info("Making query matrix...")
    query_model_matrix = make_query_model_matrix(
        metadata, ref_deseq2_dds, predictor
    )
    non_complete_cases = sum(
        !(colnames(rfu_counts) %in% rownames(query_model_matrix))
    )
    if (any(non_complete_cases)) {
        msg = sprintf(
            paste0(
                "Removing %s query samples with failed model matrix",
                " generation likely due to unseen factor levels."
            ),
            sum(non_complete_cases)
        )
        logger::log_warn(msg)
        rfu_counts = rfu_counts[, rownames(query_model_matrix), drop = FALSE]
        size_factors = size_factors[rownames(query_model_matrix)]
    }

    list(
        rfu_counts = rfu_counts,
        metadata = metadata,
        ref_deseq2_dds = ref_deseq2_dds,
        ref_glm_stats = ref_glm_stats,
        ref_model_matrix = ref_model_matrix,
        query_model_matrix = query_model_matrix,
        beta_mat = beta_mat,
        size_factors = size_factors
    )
}


#' Ensure that `target_df` has the same column order and factor levels as
#' `source_df`.
.unify_factor_levels = function(target_df, source_df) {
    target_df = target_df[, colnames(source_df), drop = FALSE]
    for (colname in colnames(source_df)) {
        if (is.factor(source_df[[colname]])) {
            target_df[[colname]] = factor(
                target_df[[colname]], levels = levels(source_df[[colname]])
            )
        }
    }

    return(target_df)
}


#' Infer coefficient name from predictor type.
.coefficient_name = function(predictor, metadata) {
    if (is.logical(metadata[[predictor]])) {
        coefficient = sprintf("%sTRUE", predictor)
    } else {
        coefficient = predictor
    }
    return(coefficient)
}


#' Make the query data model matrix.
make_query_model_matrix = function(
    metadata, ref_deseq2_dds, predictor
) {
    ref_col_data = extract_ref_col_data(ref_deseq2_dds)
    variables = colnames(ref_col_data)

    # Even if the covariate of interest is already in the query metadata, we
    # still mask it. We first set it to FALSE and make sure this column exists,
    # in order to allow the model matrix to be created. Then we mask the
    # covariate of interest to NA in the model matrix.
    if (is.logical(metadata[[predictor]])) {
        metadata[[predictor]] = FALSE
    } else {
        metadata[[predictor]] = 0
    }
    query_col_data = metadata %>%
        dplyr::select(all_of(variables)) %>%
        .unify_factor_levels(ref_col_data)
    query_model_matrix = model.matrix(
        design(ref_deseq2_dds), data = query_col_data
    ) %>%
        as.data.frame() %>%
        dplyr::rename(Intercept = "(Intercept)") %>%
        rename_with(make.names) %>%
        as.matrix()

    coefficient_of_interest = .coefficient_name(predictor, metadata)
    stopifnot(coefficient_of_interest %in% colnames(query_model_matrix))
    query_model_matrix[, coefficient_of_interest] = NA

    return(query_model_matrix)
}


#' Sanity check inputs.
sanity_check_inputs = function(
    rfu_counts, query_model_matrix, ref_model_matrix, beta_mat, size_factors
) {
    # Check that the column names of the reference model matrix == query model
    # matrix == the beta matrix.
    stopifnot(all(colnames(ref_model_matrix) == colnames(beta_mat)))
    stopifnot(all(colnames(query_model_matrix) == colnames(ref_model_matrix)))

    # Check that the RFUs match between RFU counts and beta matrix.
    if (!all(rownames(rfu_counts) == rownames(beta_mat))) {
        stop("RFUs do not match between RFU counts and beta matrix.")
    }

    # Check that the sample names match between RFU counts and query model
    # matrix.
    if (!all(colnames(rfu_counts) == rownames(query_model_matrix))) {
        stop(
            paste0(
                "Sample names do not match between RFU counts and query model",
                " matrix."
            )
        )
    }

    # Check that size factors match samples in RFU counts.
    if (length(size_factors) != ncol(rfu_counts)) {
        stop("Size factors do not match the number of samples.")
    }
    if (!all(names(size_factors) == colnames(rfu_counts))) {
        print(size_factors)
        print(colnames(rfu_counts))
        stop("Sample names do not match between size factors and RFU counts.")
    }
}


#' Extract the reference data covariates table.
#' Column `sizeFactor` is added by DESeq2 and should be removed.
extract_ref_col_data = function(ref_deseq2_dds) {
    ref_col_data = SummarizedExperiment::colData(ref_deseq2_dds) %>%
        as.data.frame() %>%
        dplyr::select(-sizeFactor)
    return(ref_col_data)
}


#' Get a list of variables to be adjusted for, which is the set of all variables
#' in the reference data that are not the covariate of interest.
.get_variables_to_adjust = function(covariates, predictor) {
    stopifnot(predictor %in% covariates)
    variables_to_adjust = setdiff(covariates, predictor)
    return(variables_to_adjust)
}


#' Adjust depth-normalized RFU counts according to a fitted GLM model's
#' coefficients.
#'
#' The purpose of this function is to adjust the depth-normalized RFU counts
#' according to the fitted coefficients of a GLM model, in order to remove the
#' fitted effect of undesired confounding covariates.
#'
#' @param normalized_counts_vec A vector of depth-normalized RFU counts with
#'   dimensions of length n, where n is the number of samples.
#' @param beta A vector of fitted coefficients.
#' @param model_matrix A matrix of covariates used to fit the model.
#' @param covariates_to_not_adjust A vector of covariates that should not be
#'   adjusted for.
#'
#' @return A vector of adjusted counts with length n.
covariate_adjusted_counts = function(
    normalized_counts_vec, beta, model_matrix, covariates_to_not_adjust
) {
    # Zero out all the covariates in the model matrix that are not to be
    # adjusted.
    for (covariate in covariates_to_not_adjust) {
        # Sanity check.
        if (!any(startsWith(colnames(model_matrix), covariate))) {
            msg = sprintf(
                "Covariate '%s' does not match model matrix coefficients %s.",
                covariate, colnames(model_matrix)
            )
            stop(msg)
        }

        model_matrix[
            ,
            colnames(model_matrix)[
                startsWith(colnames(model_matrix), covariate)
            ]
        ] = 0
    }

    # Calculate the adjusted counts.
    log2_mean = (model_matrix %*% matrix(beta, ncol = 1))[, 1]
    mu = 2^log2_mean
    adjusted_counts = normalized_counts_vec / mu
    return(adjusted_counts)
}


#' Compute GLM covariate-adjusted RFU counts.
#'
#' @param adjust_covariate A boolean matrix with dimensions R * C, where R is
#'   is the number of RFUs and C is the number of covariates. Values
#'   indicate whether a feature should be adjusted for a covariate. Must have
#'   row and column names.
#' @param beta_mat The GLM coefficients matrix with dimensions R * C. Must have
#'   row and column names.
#' @param norm_rfu_counts_mat A matrix of depth-normalized RFU counts with
#'   dimensions R * N, where R is the number of RFUs and N is the number of
#'   samples. Must have row and column names.
#' @param model_matrix The model matrix of dimensions N * C used to fit the GLM.
#'   Must have row and column names.
#'
#' @return A matrix of adjusted RFU counts with dimensions R * N.
glm_rfu_counts_adjustment = function(
    adjust_covariate,
    beta_mat,
    norm_rfu_counts_mat,
    model_matrix,
    use_pbar = TRUE
) {
    # Find the subset of RFUs that actually have something to be adjusted.
    # These RFUs have at least one covariate whose FDR is significant. Only
    # these RFUs will be adjusted in the forearch loop.
    rfu_is_adjustable = apply(adjust_covariate, 1, any)
    adjust_covariate = adjust_covariate[rfu_is_adjustable, , drop = F]
    covariate_names = colnames(adjust_covariate)

    final_counts = norm_rfu_counts_mat
    if (nrow(adjust_covariate) > 0) {
        # Adjust RFU counts one RFU at a time.
        logger::log_info("Adjusting RFU counts...")
        if (use_pbar) {
            pbar = progressr::progressor(along = rownames(adjust_covariate))
        }
        adjusted_counts_list = foreach::foreach(
            centroid = rownames(adjust_covariate),
            .options.future = list(scheduling = 10.0)  # Ten chunks per worker.
        ) %dofuture% {
            if (use_pbar) {
                pbar(message = sprintf("Adjusting RFU %s...", centroid))
            }

            # Adjust for all covariates except the ones in
            # `covariates_to_not_adjust`.
            adjusted_counts = covariate_adjusted_counts(
                normalized_counts_vec = norm_rfu_counts_mat[centroid, ],
                beta = beta_mat[centroid, ],
                model_matrix = model_matrix,
                covariates_to_not_adjust = covariate_names[
                    !adjust_covariate[centroid, ]
                ]
            )
        }

        logger::log_info("Combining adjusted counts into a matrix...")
        adjusted_counts = matrix(
            unlist(adjusted_counts_list),
            nrow = length(adjusted_counts_list),
            byrow = T
        )
        rownames(adjusted_counts) = rownames(adjust_covariate)

        # Combine the adjusted counts with the original counts of the unadjusted
        # RFUs.
        stopifnot(colnames(adjusted_counts) == colnames(final_counts))
        final_counts[rownames(adjusted_counts), ] = adjusted_counts
    }

    return(final_counts)
}


#' Compute GLM covariate-adjusted RFU counts for covariates with significant
#' association based on FDR cutoff.
compute_adjusted_rfu_counts = function(
    norm_rfu_counts,
    query_model_matrix,
    ref_glm_stats,
    beta_mat,
    variables_to_adjust,
    coefficient_of_interest,
    fdr_cutoff,
    depth_colnames,
    depth_poly_order
) {
    # Create a model matrix for the query data.
    stopifnot(all(colnames(query_model_matrix) == colnames(beta_mat)))
    query_model_matrix[, coefficient_of_interest] = 0
    stopifnot("Intercept" %in% colnames(query_model_matrix))
    query_model_matrix[, "Intercept"] = 0

    # Compute which covariates to adjust for for each RFU.
    adjust_covariate = ref_glm_stats[
        sprintf("%s.fdr", variables_to_adjust)
    ] <= fdr_cutoff
    colnames(adjust_covariate) = variables_to_adjust
    stopifnot(all(rownames(norm_rfu_counts) == rownames(adjust_covariate)))

    # The TCR depth terms come from a polynomial. Rename the columns to match
    # the polynomial coefficients.
    for (depth_colname in depth_colnames) {
        # Figure out what coefficient names this depth column name would make.
        depth_term = formula(
            sprintf(
                "~ 0 + poly(%s, %s, raw = TRUE)",
                depth_colname,
                depth_poly_order
            )
        )
        dummy_data = data.frame(1)
        colnames(dummy_data) = depth_colname
        coef_names = make.names(
            colnames(model.matrix(depth_term, data = dummy_data))
        )

        depth_covariates = sapply(
            coef_names, function(coef_name) adjust_covariate[, depth_colname]
        )

        adjust_covariate = cbind(
            adjust_covariate[
                , colnames(adjust_covariate) != depth_colname, drop = FALSE
            ],
            depth_covariates
        )
    }

    adjusted_counts = glm_rfu_counts_adjustment(
        adjust_covariate, beta_mat, norm_rfu_counts, query_model_matrix
    )

    return(adjusted_counts)
}


#' Compute log-likelihood ratio for the covariate of interest in the GLM model.
#'
#' @param raw_counts_mat Raw RFU counts matrix of RFUs * samples.
#' @param size_factors Size factors for each sample.
#' @param size Dispersion parameter for the negative binomial distribution. Note
#'   that this is the parametrization of \link{dnbinom}. The dispersion of
#'   DESeq2 alpha = 1 / size.
#' @param beta_mat Estimated coefficients matrix (RFUs * variables) for the GLM
#'   model.
#' @param model_matrix Model matrix for the GLM model.
#' @param coefficient_of_interest Variable of interest for which the
#'   log-likelihood ratio is computed.
#'
#' @return A matrix of log-likelihood ratios with dimensions RFUs * samples.
glm_log_likelihood_ratio = function(
    raw_counts_mat,
    size_factors,
    size,
    beta_mat,
    model_matrix,
    coefficient_of_interest
) {
    # Raise an error if the variable of interest is already in the model matrix.
    if (any(!is.na(model_matrix[, coefficient_of_interest]))) {
        stop("Variable of interest is already in the model matrix.")
    }

    # Raise an error if the coefficient of interest are not in the model matrix.
    if (!all(coefficient_of_interest %in% colnames(model_matrix))) {
        stop("Variable of interest is not in the model matrix.")
    }

    # Compute two versions of the model matrix for each class values.
    model_matrix_0 = model_matrix
    model_matrix_0[, coefficient_of_interest] = 0
    mu_0 = t(
        (2 ^ (model_matrix_0 %*% t(beta_mat)))
        * size_factors
    )

    model_matrix_1 = model_matrix
    model_matrix_1[, coefficient_of_interest] = 1
    mu_1 = t(
        (2 ^ (model_matrix_1 %*% t(beta_mat)))
        * size_factors
    )

    # Compute the log-likelihood ratio.
    loglik_0 = dnbinom(raw_counts_mat, size, mu = mu_0, log = T)
    loglik_1 = dnbinom(raw_counts_mat, size, mu = mu_1, log = T)
    log_lik_ratio = loglik_1 - loglik_0

    return(log_lik_ratio)
}


#' Compute adjusted RFU counts.
adjust_rfu_counts = function(
    rfu_counts,
    size_factors,
    query_model_matrix,
    ref_deseq2_dds,
    ref_glm_stats,
    beta_mat,
    glm_adjustment,
    predictor,
    fdr_cutoff,
    depth_colnames,
    depth_poly_order
) {
    coefficient_of_interest = .coefficient_name(
        predictor, colData(ref_deseq2_dds)
    )

    if (glm_adjustment == "1") {
        logger::log_info("Computing adjusted RFU counts...")

        # Get the list of variables to be adjusted for.
        ref_col_data = extract_ref_col_data(ref_deseq2_dds)
        variables_to_adjust = .get_variables_to_adjust(
            colnames(ref_col_data), predictor
        )

        # Compute the adjusted RFU counts.
        adjusted_counts = compute_adjusted_rfu_counts(
            norm_rfu_counts = rfu_counts,
            query_model_matrix = query_model_matrix,
            ref_glm_stats = ref_glm_stats,
            beta_mat = beta_mat,
            variables_to_adjust = variables_to_adjust,
            coefficient_of_interest = coefficient_of_interest,
            fdr_cutoff = fdr_cutoff,
            depth_colnames = depth_colnames,
            depth_poly_order = depth_poly_order
        )
    } else if (glm_adjustment == "log_lr") {
        logger::log_info("Computing log-likelihood ratios...")
        adjusted_counts = glm_log_likelihood_ratio(
            # We need to round the values in case the counts are not integers,
            # which currently is the case with query TCRs. See rule
            # `compute_query_rfu_counts`.
            raw_counts_mat = round(rfu_counts),
            size_factors = size_factors,
            size = 1 / dispersions(ref_deseq2_dds),
            beta_mat = beta_mat,
            model_matrix = query_model_matrix,
            coefficient_of_interest = coefficient_of_interest
        )
    } else {
        stop(
            sprintf(
                "Invalid GLM adjustment method: %s", glm_adjustment
            )
        )
    }

    return(adjusted_counts)
}


#' Function for processing a single set of files.
process_output = function(
    rfu_counts,
    size_factors,
    metadata,
    ref_deseq2_dds,
    ref_glm_stats,
    ref_model_matrix,
    query_model_matrix,
    beta_mat,
    glm_adjustment,
    predictor,
    fdr_cutoff,
    depth_colnames,
    depth_poly_order,
    output_tsv
) {
    # Sanity checks.
    sanity_check_inputs(
        rfu_counts, query_model_matrix, ref_model_matrix, beta_mat, size_factors
    )

    if (nrow(ref_glm_stats) > 0) {
        # Adjust the RFU counts.
        adjusted_counts = adjust_rfu_counts(
            rfu_counts,
            size_factors,
            query_model_matrix,
            ref_deseq2_dds,
            ref_glm_stats,
            beta_mat,
            glm_adjustment,
            predictor,
            fdr_cutoff,
            depth_colnames,
            depth_poly_order
        )
    } else {
        logger::log_info("No RFU counts to adjust.")
        adjusted_counts = rfu_counts[
            -seq_len(nrow(rfu_counts)), , drop = FALSE
        ]
    }

    fh = gzfile(output_tsv, "w")
    write.table(
        adjusted_counts,
        fh,
        sep = "\t",
        quote = FALSE,
        col.names = TRUE,
        row.names = TRUE
    )
    close(fh)
}


verify_all_same_length = function(lst) {
    lengths = sapply(lst, length)
    stopifnot(all(lengths == lengths[1]))
}


#' Find and trim common prefixes in a list of strings.
trim_common_prefix = function(strings) {
    shortest_length = min(nchar(strings))
    for (prefix_len in rev(seq_len(shortest_length))) {
        prefixes = unique(substr(strings, 1, prefix_len))
        if (all(prefixes == prefixes[1])) {
            return(substr(strings, prefix_len + 1, nchar(strings)))
        }
    }
}


main = function(
    rfu_counts_tsv,
    size_factors_tsv,
    metadata_tsv,
    ref_deseq2_dds,
    ref_glm_stats_tsv,
    glm_adjustment,
    predictor,
    fdr_cutoff,
    depth_colnames,
    depth_poly_order,
    output_tsv
) {
    cached_input_reader = create_cached_input_reader()

    # Compile the list of files to process and order them by the following files
    # in the nested order of the number of wildcards.
    file_list = list(
        rfu_counts_tsv = rfu_counts_tsv,
        ref_deseq2_dds = ref_deseq2_dds,
        ref_glm_stats_tsv = ref_glm_stats_tsv,
        output_tsv = output_tsv
    )
    verify_all_same_length(file_list)
    filename_df = as.data.frame(file_list) %>%
        arrange(across(all_of(names(file_list))))

    # Process each output file using its respective input files.
    for (i in seq_len(nrow(filename_df))) {
        rfu_counts_tsv = filename_df[i, "rfu_counts_tsv"]
        ref_deseq2_dds = filename_df[i, "ref_deseq2_dds"]
        ref_glm_stats_tsv = filename_df[i, "ref_glm_stats_tsv"]
        output_tsv = filename_df[i, "output_tsv"]

        trimmed_fnames = trim_common_prefix(
            c(output_tsv, rfu_counts_tsv, ref_deseq2_dds, ref_glm_stats_tsv)
        )
        msg = sprintf(
            "(%s/%s) Reading inputs for output %s: %s",
            i,
            nrow(filename_df),
            trimmed_fnames[1],
            paste(trimmed_fnames[-1], collapse = ", ")
        )
        logger::log_info(msg)

        inputs = load_inputs(
            cached_input_reader,
            rfu_counts_tsv,
            size_factors_tsv,
            metadata_tsv,
            ref_deseq2_dds,
            ref_glm_stats_tsv,
            predictor
        )
        rfu_counts = inputs$rfu_counts
        size_factors = inputs$size_factors
        metadata = inputs$metadata
        ref_deseq2_dds = inputs$ref_deseq2_dds
        ref_glm_stats = inputs$ref_glm_stats
        ref_model_matrix = inputs$ref_model_matrix
        query_model_matrix = inputs$query_model_matrix
        beta_mat = inputs$beta_mat

        process_output(
            rfu_counts,
            size_factors,
            metadata,
            ref_deseq2_dds,
            ref_glm_stats,
            ref_model_matrix,
            query_model_matrix,
            beta_mat,
            glm_adjustment,
            predictor,
            fdr_cutoff,
            depth_colnames,
            depth_poly_order,
            output_tsv
        )
    }
}


if (!interactive()) {
    # Set progressr options.
    progress_format = paste(
        ":spin :current/:total (:message) [:bar] :percent in :elapsed ETA:",
        " :eta",
        collapse = ""
    )
    progressr::handlers(
        progressr::handler_progress(
            format = progress_format, interval = 5, enable = TRUE
        )
    )
    progressr::handlers(global = TRUE)
    options(progressr.enable = TRUE)

    # Set parallelization options.
    future::plan(multicore)
    options(
        mc.cores = snakemake@threads,
        future.globals.maxSize = 21474836480  # 20 * (1024 ^ 3) bytes = 20 GB.
    )

    main(
        rfu_counts_tsv = snakemake@input[["rfu_counts_tsv"]],
        size_factors_tsv = snakemake@input[["size_factors_tsv"]],
        metadata_tsv = snakemake@input[["metadata_tsv"]],
        ref_deseq2_dds = snakemake@input[["deseq2_obj"]],
        ref_glm_stats_tsv = snakemake@input[["glm_stats_tsv"]],
        glm_adjustment = snakemake@params[["glm_adjustment"]],
        predictor = snakemake@wildcards[["predictor"]],
        fdr_cutoff = snakemake@params[["fdr_cutoff"]],
        depth_colnames = snakemake@params[["depth_colnames"]],
        depth_poly_order = snakemake@params[["depth_poly_order"]],
        output_tsv = unlist(snakemake@output)
    )
}

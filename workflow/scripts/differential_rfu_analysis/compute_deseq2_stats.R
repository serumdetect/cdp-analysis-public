#' Snakemake script for performing "differential RFU analysis" using DESeq2

suppressMessages(library(DESeq2))
suppressMessages(library(apeglm))


read_input_data = function(
    counts_mat_tsv,
    size_factors_tsv,
    metadata_tsv,
    min_subject_count,
    depth_colnames,
    disease_subset,
    predictor,
    covariates
) {
    metadata = read.table(
        metadata_tsv,
        header = T,
        sep = "\t",
        row.names = "sample_name",
        comment = "",
        check.names = F,
        quote = "\""
    )

    # Subset individuals to use by disease. Prepare the predictors data frame.
    if (!is.null(disease_subset)) {
        disease_subset = match.arg(
            disease_subset,
            c("cancer", "control", "positive", "negative")
        )
        if (disease_subset %in% c("cancer", "positive")) {
            metadata = metadata[metadata[[predictor]], ]
        } else {
            metadata = metadata[!metadata[[predictor]], ]
        }
    }

    # Remove samples with incomplete information.
    col_data_columns = c(predictor, covariates, depth_colnames)
    complete_cases = complete.cases(metadata[col_data_columns])
    message(
        sprintf(
            "Removing %s individuals with incomplete covariates.",
            sum(!complete_cases)
        )
    )
    metadata = metadata[complete_cases, , drop = F]

    # The covariates for the model are the specified covariates plus the TCR
    # repertoire depth.
    col_data = metadata[, col_data_columns, drop = F]

    # Read the RFU counts table.
    temp = read.table(
        counts_mat_tsv, header = T, row.names = 1, sep = "\t", check.names = F
    )
    counts_mat = as.matrix(temp)
    rownames(counts_mat) = rownames(temp)
    colnames(counts_mat) = colnames(temp)
    rm(temp)

    # Consolidate counts matrix columns and col_data rows.
    col_data = col_data[rownames(col_data) %in% colnames(counts_mat), , F]
    counts_mat = counts_mat[, rownames(col_data)]

    # Filter RFUs by frequency *after* the samples have been filtered by
    # completeness.
    selected_rfus = apply(
        counts_mat,
        1,
        function(rfu_counts) {
            return(sum(rfu_counts > 0) >= min_subject_count)
        }
    )
    counts_mat = counts_mat[selected_rfus, ]

    # Read the size factors.
    size_factors = read.table(size_factors_tsv, header = F, row.names = 1)
    size_factors = setNames(size_factors[[1]], rownames(size_factors))
    size_factors = size_factors[colnames(counts_mat)]

    list(counts = counts_mat, col_data = col_data, size_factors = size_factors)
}


#' Create a string representation of the polynomial term for TCR repertoire
#' depth.
create_depth_poly_string = function(depth_colname, depth_poly_order) {
    return(
        sprintf(
            "poly(%s, %s, raw = TRUE)",
            depth_colname,
            depth_poly_order
        )
    )
}


#' Create a model design that includes a polynomial term for the TCR repertoire
#' depth.
create_full_design = function(col_data, depth_colnames, depth_poly_order) {
    stopifnot(depth_colnames %in% colnames(col_data))
    tcr_depth_terms = paste(
        sapply(
            depth_colnames,
            function(depth_colname) {
                create_depth_poly_string(depth_colname, depth_poly_order)
            }
        ),
        collapse = " + "
    )
    design = update(
        reformulate(setdiff(colnames(col_data), depth_colnames)),
        formula(sprintf("~ . + %s", tcr_depth_terms))
    )
    return(design)
}


#' Helper function for creating the reduced model design.
#'
#' TCR repertoire depth is a special variable that creates a polynomial term in
#' the full design. Therefore, the reduced design should exclude all of these
#' terms.
create_reduced_design = function(
    design_full, covariate_to_remove, depth_colnames, depth_poly_order
) {
    if (!(covariate_to_remove %in% depth_colnames)) {
        design_reduced = update(
            design_full, sprintf("~ . - %s", covariate_to_remove)
        )
    } else {
        design_reduced = update(
            design_full,
            sprintf(
                "~ . - %s",
                create_depth_poly_string(covariate_to_remove, depth_poly_order)
            )
        )
    }
    return(design_reduced)
}


#' Helper function for selecting the variables to use in the DESeq2 results.
.make_deseq2_results_variables = function(
    col_data, covariate_label, depth_colnames, depth_poly_order
) {
    if (is.factor(col_data[[covariate_label]])) {
        variables = paste(
            covariate_label,
            levels(col_data[[covariate_label]])[2],
            sep = ""
        )
    } else if (is.logical(col_data[[covariate_label]])) {
        variables = paste(covariate_label, "TRUE", sep = "")
    } else if (covariate_label %in% depth_colnames) {
        # Get the names that model.matrix will make.
        depth_terms = as.formula(
            sprintf(
                "~ 0 + %s",
                create_depth_poly_string(covariate_label, depth_poly_order)
            )
        )
        variables = make.names(
            colnames(model.matrix(depth_terms, data = col_data))
        )
    } else {
        variables = covariate_label
    }
    return(variables)
}


#' Compute DESeq2 statistics.
run_deseq2 = function(
    counts_mat,
    col_data,
    size_factors,
    depth_colnames,
    depth_poly_order
) {
    message(
        sprintf(
            "Running DESeq2 on a matrix of shape (%s)",
            paste(dim(counts_mat), collapse = ", ")
        )
    )
    # Run the full DESeq2 model using all predictors.
    design_full = create_full_design(col_data, depth_colnames, depth_poly_order)
    dds = DESeqDataSetFromMatrix(
        countData = counts_mat,
        colData = col_data,
        design = design_full
    )

    message("Using pre-computed size factors...")
    sizeFactors(dds) = size_factors[colnames(dds)]

    # Estimate dispersions manually.
    dds = estimateDispersions(dds, fitType = "glmGamPoi")
    dds = estimateDispersionsMAP(dds, outlierSD = Inf, type = "glmGamPoi")

    # Compute reduced models for each covariate in order to compute LRT P value.
    # `first_dds` is the first reduced model, which is to be saved.
    results_of = list()
    first_dds = NULL
    for (covariate_label in colnames(col_data)) {
        message(
            sprintf("Computing LRT P value for predictor '%s'", covariate_label)
        )
        design_reduced = create_reduced_design(
            design_full, covariate_label, depth_colnames, depth_poly_order
        )

        # Manually compute LRT.
        DESeq2:::checkLRT(design_full, design_reduced)
        dds_reduced = nbinomLRT(
            dds,
            full = design_full,
            reduced = design_reduced,
            type = "glmGamPoi"
        )
        variables = .make_deseq2_results_variables(
            dds_reduced@colData,
            covariate_label,
            depth_colnames,
            depth_poly_order
        )
        results = results(
            dds_reduced, list(variables), cooksCutoff = FALSE
        )
        results[["fdr"]] = p.adjust(results[["pvalue"]], method = "fdr")
        results_of[[covariate_label]] = results

        # `first_dds` is the model computed with the first predictor, which is
        # assumed to be the main variable of interest.
        if (is.null(first_dds)) {
            first_dds = dds_reduced
        }
    }

    # Combine results into a data frame and order by the P value of (of the
    # first predictor).
    results = do.call(data.frame, results_of)
    first_covariate_pval_column = paste(
        names(results_of)[1], "pvalue", sep = "."
    )
    results = results[order(results[[first_covariate_pval_column]]), ]

    list(dds = first_dds, result = results)
}


main = function(
    counts_mat_tsv,
    size_factors_tsv,
    metadata_tsv,
    min_subject_count,
    depth_colnames,
    depth_poly_order,
    disease_subset,
    predictor,
    covariates,
    deseq2_output_rds,
    glm_output_tsv,
    coef_mat_tsv
) {
    message("Reading input files...")
    inputs = read_input_data(
        counts_mat_tsv,
        size_factors_tsv,
        metadata_tsv,
        min_subject_count,
        depth_colnames,
        disease_subset,
        predictor,
        covariates
    )
    counts_mat = inputs[["counts"]]
    col_data = inputs[["col_data"]]
    size_factors = inputs[["size_factors"]]

    # Make sure col_data columns are valid R variables.
    colnames(col_data) = make.names(colnames(col_data))

    message(sprintf("Running DESeq2 for %s RFUs...", nrow(counts_mat)))
    deseq2_result = run_deseq2(
        counts_mat, col_data, size_factors, depth_colnames, depth_poly_order
    )
    dds = deseq2_result[["dds"]]
    deseq2_stats = deseq2_result[["result"]]
    coef_mat = coef(dds)

    message("Saving DESeq2 results...")
    if (!is.null(deseq2_output_rds)) {
        saveRDS(dds, deseq2_output_rds)
    }

    output_gzfile = gzfile(glm_output_tsv, "w")
    write.table(
        deseq2_stats,
        output_gzfile,
        col.names = T,
        row.names = T,
        sep = "\t",
        quote = F
    )
    close(output_gzfile)

    output_gzfile = gzfile(coef_mat_tsv, "w")
    write.table(
        coef_mat,
        output_gzfile,
        col.names = T,
        row.names = T,
        sep = "\t",
        quote = F
    )
    close(output_gzfile)
}


if (!interactive()) {
    main(
        counts_mat_tsv = snakemake@input[["counts_mat"]],
        size_factors_tsv = snakemake@input[["size_factors"]],
        metadata_tsv = snakemake@input[["metadata"]],
        min_subject_count = snakemake@params[["min_subject_count"]],
        depth_colnames = snakemake@params[["depth_colnames"]],
        depth_poly_order = snakemake@params[["depth_poly_order"]],
        disease_subset = snakemake@params[["disease_subset"]],
        predictor = snakemake@wildcards[["predictor"]],
        covariates = snakemake@params[["covariates"]],
        deseq2_output_rds = snakemake@output[["deseq2_obj"]],
        glm_output_tsv = snakemake@output[["glm_stats_tsv"]],
        coef_mat_tsv = snakemake@output[["coef_mat_tsv"]]
    )
}

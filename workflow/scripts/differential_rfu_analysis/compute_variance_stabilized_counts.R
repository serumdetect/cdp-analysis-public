#' Snakemake script for computing variance-stabilized counts.

library(dplyr)
library(transformGamPoi)


main = function(
    rfu_counts_tsv,
    size_factors_tsv,
    min_subject_count,
    output_rfu_counts_tsv,
    output_mu_dispersion_tsv
) {
    # Read inputs.
    rfu_counts = read.table(
        rfu_counts_tsv,
        header = TRUE,
        row.names = 1,
        check.names = FALSE
    )
    rfu_counts = as.matrix(
        rfu_counts[
            rowSums(rfu_counts > 0) >= min_subject_count, , drop = FALSE
        ],
    )

    size_factors = read.table(
        size_factors_tsv,
        header = FALSE,
        check.names = FALSE,
        col.names = c("sample_name", "size_factor")
    ) %>%
        with(setNames(size_factor, sample_name))
    size_factors = size_factors[colnames(rfu_counts)]
    stopifnot(all(colnames(rfu_counts) == names(size_factors)))

    # Compute variance-stabilized counts.
    fit = transformGamPoi::residual_transform(
        rfu_counts,
        residual_type = "pearson",
        overdispersion = TRUE,
        size_factors = size_factors,
        offset_model = TRUE,
        overdispersion_shrinkage = TRUE,
        ridge_penalty = 0,
        return_fit = TRUE,
        subsample = 100000,
        verbose = TRUE
    )

    # Write outputs.
    output_rfu_counts_tsv_gz = gzfile(output_rfu_counts_tsv, "w")
    write.table(
        fit$Residuals,
        output_rfu_counts_tsv_gz,
        sep = "\t",
        col.names = NA,
        row.names = TRUE,
        quote = FALSE
    )
    close(output_rfu_counts_tsv_gz)

    output_mu_dispersion_tsv_gz = gzfile(output_mu_dispersion_tsv, "w")
    if (!("Intercept" %in% colnames(fit$fit$Beta))) {
        stop("Intercept not found in fit$Beta")
    }
    write.table(
        data.frame(
            beta_0 = fit$fit$Beta[, "Intercept"],
            overdispersions = fit$fit$overdispersions
        ),
        output_mu_dispersion_tsv_gz,
        sep = "\t",
        col.names = NA,
        row.names = TRUE,
        quote = FALSE
    )
    close(output_mu_dispersion_tsv_gz)
}


if (!interactive()) {
    main(
        rfu_counts_tsv = snakemake@input$rfu_counts_tsv,
        size_factors_tsv = snakemake@input$size_factors_tsv,
        min_subject_count = snakemake@params$min_subject_count,
        output_rfu_counts_tsv = snakemake@output$output_rfu_counts_tsv,
        output_mu_dispersion_tsv = snakemake@output$output_mu_dispersion_tsv
    )
}

#' Snakemake script for computing variance-stabilized counts for query
#' RFU counts based on model fitted on reference RFU counts.


library(dplyr)
library(glmGamPoi)


#' Calculate the gamma-Poisson mean based on model mean and size factors.
.calc_glm_mu = function(log_mu, size_factors) {
    mu = outer(exp(log_mu), size_factors, "*")
    return(mu)
}


#' Function for computing Pearson residuals as in glmGamPoi::residual_transform.
.calc_pearson_resid = function(rfu_counts, mu, overdispersions) {
    stopifnot(
        all(rownames(rfu_counts) == rownames(mu)),
        all(colnames(rfu_counts) == colnames(mu))
    )

    residuals = glmGamPoi:::div_zbz_dbl_mat(
        rfu_counts - mu,
        sqrt(
            mu + glmGamPoi:::multiply_vector_to_each_column(
                mu^2, overdispersions
            )
        )
    )
    rownames(residuals) = rownames(rfu_counts)
    colnames(residuals) = colnames(rfu_counts)
    return(residuals)
}


main = function(
    rfu_counts_tsv,
    size_factors_tsv,
    reference_mu_dispersion_tsv,
    output_rfu_counts_tsv
) {
    # Read inputs. Only keep RFUs that are present in the reference.
    reference_mu_dispersion = read.table(
        reference_mu_dispersion_tsv,
        header = TRUE,
        row.names = 1,
        check.names = FALSE
    )

    rfu_counts = read.table(
        rfu_counts_tsv, header = TRUE, row.names = 1, check.names = FALSE
    ) %>%
        as.matrix()

    rfu_counts = rfu_counts[
        rownames(rfu_counts) %in% rownames(reference_mu_dispersion),
        ,
        drop = FALSE
    ]
    reference_mu_dispersion = reference_mu_dispersion[
        rownames(rfu_counts), , drop = FALSE
    ]

    size_factors = read.table(
        size_factors_tsv,
        header = FALSE,
        check.names = FALSE,
        col.names = c("sample_name", "size_factor")
    ) %>%
        with(setNames(size_factor, sample_name))
    size_factors = size_factors[colnames(rfu_counts)]
    stopifnot(all(colnames(rfu_counts) == names(size_factors)))

    # Calculate query residuals.
    query_mu = .calc_glm_mu(reference_mu_dispersion[, "beta_0"], size_factors)
    rownames(query_mu) = rownames(reference_mu_dispersion)
    colnames(query_mu) = names(size_factors)
    query_resid = .calc_pearson_resid(
        rfu_counts, query_mu, reference_mu_dispersion[["overdispersions"]]
    )

    # Write outputs.
    output_rfu_counts_tsv_gz = gzfile(output_rfu_counts_tsv, "w")
    write.table(
        query_resid,
        output_rfu_counts_tsv_gz,
        sep = "\t",
        col.names = NA,
        row.names = TRUE,
        quote = FALSE
    )
    close(output_rfu_counts_tsv_gz)
}


if (!interactive()) {
    main(
        snakemake@input$rfu_counts_tsv,
        snakemake@input$size_factors_tsv,
        snakemake@input$reference_mu_dispersion_tsv,
        snakemake@output$output_rfu_counts_tsv
    )
}

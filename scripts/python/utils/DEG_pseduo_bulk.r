library(edgeR)
library(MAST)

#' Differential Gene Expression Analysis with edgeR
#'
#' This function performs a complete differential gene expression (DEG) analysis
#' using the edgeR package. It handles pseudo-bulk count data for a specific cell
#' population, fits a quasi-likelihood negative binomial model, and generates
#' diagnostic plots and a results table. It includes optional batch effect correction.
#'
#' @param expr_df A data frame or matrix of raw gene counts with genes in rows and samples in columns.
#' @param coldata_df A data frame containing metadata for each sample.
#' @param group_col The column name in `coldata_df` specifying the experimental group or condition.
#' @param DEGDesign A character string defining the contrast for the DEG comparison (e.g., "GroupA-GroupB").
#' @param out_prefix A character string used as a prefix for all output filenames (plots and tables).
#' @param fig_dir The directory path where diagnostic plots will be saved.
#' @param table_dir The directory path where the DEG results table will be saved.
#' @param cell_identity_key The column name in `coldata_df` for cell type annotations (default: "celltype").
#' @param replicate_col The column name for replicate or batch information (optional, default: NULL).
#'
#' @return A list containing the fitted model and processed data:
#' \itemize{
#'   \item \strong{fit}: The `glmQLFit` object containing the fitted model.
#'   \item \strong{design}: The design matrix used for the model fit.
#'   \item \strong{y}: The processed `DGEList` object after normalization and filtering.
#' }
#'
#' @details The function first combines the `group_col` and `cell_identity_key` columns
#'   to create a unique factor for each cell population and experimental group combination.
#'   It then filters out lowly expressed genes and normalizes the data using TMM normalization.
#'   A design matrix is built based on the `combined_group` factor, with an option to include
#'   a `replicate_col` for batch effect correction. The function then fits a
#'   quasi-likelihood negative binomial generalized linear model to the data.
#'
#'   Finally, it performs the differential expression test based on the `DEGDesign` contrast,
#'   saves the full DEG results table to a TSV file, and generates and saves three
#'   diagnostic plots: a Multi-Dimensional Scaling (MDS) plot for sample similarity,
#'   a Biological Coefficient of Variation (BCV) plot, and a Smear plot highlighting
#'   genes with an FDR less than 0.01.
#'
#' @seealso \code{\link[edgeR]{glmQLFit}}, \code{\link[edgeR]{plotMDS}}, \code{\link[edgeR]{plotBCV}},
#'   \code{\link[limma]{makeContrasts}}
#'
#' @examples
#' # Example usage (assuming expr_df and coldata_df are prepared)
#' # fit_model(
#' #   expr_df = counts_matrix,
#' #   coldata_df = sample_metadata,
#' #   group_col = "treatment",
#' #   DEGDesign = "treated-untreated",
#' #   out_prefix = "my_analysis",
#' #   fig_dir = "results/plots",
#' #   table_dir = "results/tables"
#' # )
fit_model <- function(
    expr_df, 
    coldata_df,
    group_col,
    DEGDesign,
    out_prefix,
    fig_dir,
    table_dir,
    cell_identity_key = "celltype",
    replicate_col = NULL  # replicate is optional
){
    # Check columns exist
    if (!(group_col %in% colnames(coldata_df))) stop("group_col not found in coldata_df")
    if (!(cell_identity_key %in% colnames(coldata_df))) stop("cell_identity_key not found in coldata_df")

    # Combine group + cell type for filtering
    combined_group <- factor(paste0(coldata_df[[group_col]], "_", coldata_df[[cell_identity_key]]))

    # Check that each group has at least 2 samples
    if (any(table(combined_group) < 2)) {
        warning("Some groups have fewer than 2 samples; glmQLFit may fail")
    }

    # Create DGEList using combined group
    y <- DGEList(counts = as.matrix(expr_df), group = combined_group)

    # Filter lowly expressed genes
    cat("Dimensions before filtering lowly expressed genes:\n")
    print(dim(y))
    keep <- filterByExpr(y)
    y <- y[keep, , keep.lib.sizes = FALSE]
    cat("Dimensions after filtering:\n")
    print(dim(y))

    # Normalize counts
    y <- calcNormFactors(y)

    # Build design matrix
    if (is.null(replicate_col) || !(replicate_col %in% colnames(coldata_df)) || all(is.na(coldata_df[[replicate_col]]))) {
        design <- model.matrix(~0 + combined_group)
        print("Design matrix includes only group.")
    } else {
        replicate <- factor(coldata_df[[replicate_col]])
        print(replicate)
        design <- model.matrix(~0 + combined_group + replicate)
        print(colnames(design))
        print("Design matrix includes group and replicate.")
    }

    # Estimate dispersion
    y <- estimateDisp(y, design = design)

    # Fit quasi-likelihood negative binomial model
    fit <- glmQLFit(y, design)

    png(filename = paste0(fig_dir, "/", out_prefix, "_MDS.png"))
    plotMDS(y)
    dev.off()
    png(filename = paste0(fig_dir, "/", out_prefix, "_BCV.png"))
    plotBCV(y)
    dev.off()
    print("-----1-----")
    print(DEGDesign)
    myContrast <- do.call(
    what = 'makeContrasts',
    args = list(DEGDesign, levels = design)
    )

    qlf <- glmQLFTest(fit, contrast=myContrast)
    # get all of the DE genes and calculate Benjamini-Hochberg adjusted FDR
    tt <- topTags(qlf, n = Inf)
    tt <- tt$table
    write.table(tt, file = paste0(table_dir,"/",out_prefix,"_DEG.tsv"), sep = "\t", quote = FALSE, row.names = TRUE)

    # tr <- glmTreat(fit, contrast=myContrast, lfc=1.5)
    png(filename = paste0(fig_dir, "/", out_prefix, "_Smear_0.01FDR.png"))
    plotSmear(qlf, de.tags = rownames(tt)[which(tt$FDR<0.01)])
    dev.off()

    return(list(
        fit = fit,
        design = design,
        y = y
    ))
}

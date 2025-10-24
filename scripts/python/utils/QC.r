library(decontX)
library(SingleCellExperiment)

#' ambientRNA: Decontaminate ambient RNA in single-cell data
#'
#' This function takes a count matrix and cluster assignments as data frames,
#' runs the decontX algorithm from the celda package to estimate and remove
#' ambient RNA contamination, and outputs decontaminated counts and contamination scores.
#'
#' @param count_df A data frame of raw counts (genes x cells). Row names should be gene names.
#' @param cluster_df A data frame containing cluster assignments for cells. Row names should match columns of count_df.
#' @param cluster_col The column name in cluster_df that contains cluster labels. Default is "leiden".
#' @param prefix A character string to prepend to output files. Default is "output".
#'
#' @return A SingleCellExperiment object with decontaminated counts stored in assay "decontXcounts"
#'         and contamination scores in sce$decontX_contamination.
#'         Also writes two CSV files:
#'           - <prefix>_decontaminated_counts.csv
#'           - <prefix>_contamination.csv
#'
#' @examples
#' # counts_df: genes x cells raw counts
#' # clusters_df: data frame with a column "leiden" for cluster assignments
#' sce <- ambientRNA(counts_df, clusters_df, cluster_col = "leiden", prefix = "sample1")
#'
ambientRNA <- function(count_df, cluster_df, cluster_col = "leiden", prefix = "output") {
    print("-------1------------")
    # Convert count data frame to matrix
    counts <- as.matrix(count_df)
    clusters <- as.data.frame(cluster_df)
    
    # Create SingleCellExperiment object
    sce <- SingleCellExperiment(assays = list(counts = counts))
    sce$cluster <- factor(clusters[[cluster_col]])
    
    # Run decontX to remove ambient RNA contamination
    sce <- decontX(sce)
    print("-------2------------")
    # # Construct output file names
    # counts_file <- paste0(prefix, "_decontaminated_counts.csv")
    # contamination_file <- paste0(prefix, "_contamination.csv")
    decountXMatrix = as.matrix(assay(sce, "decontXcounts")) # The complex transformation is left to r
    # # Write results to CSV
    # write.csv(assay(sce, "decontXcounts"), counts_file)
    # write.csv(sce$decontX_contamination, contamination_file)
    
    # print("DecontX finished. Results saved to:")
    # print(counts_file)
    # print(contamination_file)
    
    return(decountXMatrix)
}

ambientRNA1 <- function(count_file, cluster_file, cluster_col = "leiden", prefix = "output") {
    print("--------------1----------")
    counts_df <- read.csv(count_file, row.names = 1, check.names = FALSE)
    clusters_df <- read.csv(cluster_file, row.names = 1, check.names = FALSE)
    # Convert count data frame to matrix
    counts <- as.matrix(counts_df)
    clusters <- as.data.frame(clusters_df)
    
    # Create SingleCellExperiment object
    sce <- SingleCellExperiment(assays = list(counts = counts))
    cluster_vec <- as.character(clusters[cluster_col, ])
    names(cluster_vec) <- colnames(clusters)  
    sce$cluster <- factor(cluster_vec)
    
    # Run decontX to remove ambient RNA contamination
    sce <- decontX(sce)
    print("-------2------------")
    # Construct output file names
    counts_file <- paste0(prefix, "_decontaminated_counts.csv")
    contamination_file <- paste0(prefix, "_contamination.csv")
    
    # Write results to CSV
    write.csv(assay(sce, "decontXcounts"), counts_file)
    write.csv(sce$decontX_contamination, contamination_file)
    
    print("DecontX finished. Results saved to:")
    print(counts_file)
    print(contamination_file)
    return(sce)
}
# count_file = "/disk5/luosg/scRNAseq/count.csv"
# cluster_file = "/disk5/luosg/scRNAseq/cluster.csv"
# ambientRNA1(count_file,cluster_file)

library(reticulate)
library(zellkonverter)
library(SingleCellExperiment)
library(argparse)
# 定义函数
convert_h5ad_to_sce <- function(input_h5ad, output_rds, verbose = TRUE, python_env = NULL) {
  # 可选：手动指定 Python 环境
  if (!is.null(python_env)) {
    reticulate::use_condaenv(python_env, required = TRUE)
  }
  
  # 导入 anndata 包
  ad <- reticulate::import("anndata", convert = FALSE)
  
  if (verbose) message("Reading h5ad file: ", input_h5ad)
  adata_py <- ad$read_h5ad(input_h5ad)
  
  if (verbose) message("Converting AnnData to SingleCellExperiment...")
  sce <- zellkonverter::AnnData2SCE(adata_py, verbose = verbose)
  
  if (verbose) message("Saving SingleCellExperiment to RDS: ", output_rds)
  saveRDS(sce, file = output_rds)
  
  if (verbose) message("Done.")
  return(sce)
}
# -------------------------------
parser <- ArgumentParser(description = "Convert .h5ad to SingleCellExperiment and save as RDS")

parser$add_argument("input_h5ad", help = "Path to input h5ad file")
parser$add_argument("output_rds", help = "Path to output RDS file")
parser$add_argument("--python_env", help = "Optional Python conda environment", default = NULL)
parser$add_argument("--verbose", help = "Verbose output", action = "store_true")

args <- parser$parse_args()

# -------------------------------
# 调用函数
# -------------------------------
convert_h5ad_to_sce(
  input_h5ad = args$input_h5ad,
  output_rds = args$output_rds,
  python_env = args$python_env,
  verbose = args$verbose
)

from pathlib import Path
import sys
import os
basePath = Path(__file__).resolve()
baseDir = basePath.parent
sys.path.append(str(baseDir / "utils"))
import anndata as ad
import pandas as pd
import scanpy as sc
import logging
import rpy2.robjects as ro
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import pandas2ri
import rpy2.rinterface_lib.callbacks
from rpy2.robjects.packages import importr
import rpy2.robjects.numpy2ri as numpy2ri
from scipy import sparse
rpy2.rinterface_lib.callbacks.logger.setLevel(logging.INFO)
sc.settings.verbosity = 0
sce_pkg = importr('SingleCellExperiment')
logging.basicConfig(
    level=logging.INFO, 
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logging.info(ro.r("R.version.string"))
os.environ['R_HOME'] = "/home/luosg/miniconda3/envs/scRNAseq_rpy2/lib/R"
def ambientRNA_Remove(
    adata:ad.AnnData,
    out_prefix:str,
    r_script_path:str,
    cluster_col:str = "leiden",
):
    X = adata.X if adata.raw is None else adata.raw.X

    logging.info("chekck whether has toarray methods")
    if hasattr(X, "toarray"):
        adata.layers["raw"] = adata.X
        logging.info("convert it to dense matrix to create dense dataframe")
        X = X.toarray()
    else:
        adata.layers["raw"] = sparse.csr_matrix(adata.X)
    
    logging.info("create dataframe AnnData (Cells x Genes) -> R (Genes x Cells)")
    df_count = pd.DataFrame(
        X.T,
        index=adata.var_names if adata.raw is None else adata.raw.var_names,
        columns=adata.obs_names if adata.raw is None else adata.raw.obs_names
    )
    
    if cluster_col not in adata.obs.columns:
        logging.info(f"Cluster column '{cluster_col}' not found. Recomputing PCA, neighbors, and clustering.")
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata)
        sc.pp.pca(adata)
        sc.pp.neighbors(adata)
        sc.tl.leiden(adata)
    df_cluster = adata.obs[[cluster_col]]
    ro.r(f'source("{r_script_path}")')
    r_ambientRNA = ro.globalenv['ambientRNA']
    logging.info("Convert DataFrame from Python to R ")
    
    with localconverter(ro.default_converter + pandas2ri.converter):
        r_count = ro.conversion.py2rpy(df_count)
        r_cluster = ro.conversion.py2rpy(df_cluster)

    r_counts_matrix = r_ambientRNA(
        count_df = r_count,
        cluster_df = r_cluster,
        cluster_col = "leiden",
        prefix = out_prefix
    )
    logging.info(f"r_counts_matrix 类型： {ro.r['class'](r_counts_matrix)}")
    if r_counts_matrix is not None:
        try:
            with localconverter(ro.default_converter + numpy2ri.converter):
                py_corrected_counts = ro.conversion.rpy2py(r_counts_matrix)
            if py_corrected_counts.shape[0] == adata.n_vars and py_corrected_counts.shape[1] == adata.n_obs:
                logging.info("R (Genes x Cells) -> AnnData (Cells x Genes)")
                py_corrected_counts = py_corrected_counts.T
            if py_corrected_counts.shape == adata.shape:
                adata.X = sparse.csr_matrix(py_corrected_counts)
                adata.write_h5ad(f"{out_prefix}_deambiendRNA.h5ad")
                logging.info(f"Successfully injected corrected counts into adata.X with shape {adata.X.shape}. Saved it to {out_prefix}_deambiendRNA.h5ad")
            else:
                logging.error(f"Shape mismatch: Corrected counts ({py_corrected_counts.shape}) must match AnnData shape ({adata.shape}).")
        except Exception as e:
            logging.error(f"Error during R matrix conversion or injection: {e}")
    return adata
if __name__ == "__main__":
    infileDict = {"Muscle":"/disk5/luosg/scRNAseq/output/result/h5ad/Muscle_qc.h5ad",
              "Intestine": "/disk5/luosg/scRNAseq/output/result/h5ad/Intestine_qc.h5ad",
           "Lung":"/disk5/luosg/scRNAseq/output/result/h5ad/Lung_qc.h5ad"}
    outdir = "/disk5/luosg/scRNAseq/output/result/h5ad"
    for tissue,infile in infileDict.items():
        adata = sc.read_h5ad(infile)
        outfilePrefix = f"{outdir}/{tissue}_qc"
        logging.info(f"tissue: {tissue}, infile: {infile}, outfile: {outfilePrefix}_deambiendRNA.h5ad")
        adata = sc.read_h5ad(infile)
        ambientRNA_Remove(adata,outfilePrefix,r_script_path="/disk5/luosg/scRNAseq/workflow/scRNAseq/scripts/python/utils/QC.r")
        

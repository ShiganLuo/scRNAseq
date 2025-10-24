from ast import Raise
import scanpy as sc
import matplotlib.pyplot as plt
import datetime
import argparse
import anndata as ad
from typing import Literal
import logging
from pathlib import Path
logging.basicConfig(

  level=logging.INFO,

  format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',

  datefmt='%Y-%m-%d %H:%M:%S'

)
def setup_output_directory(flag: bool, directory: str | Path | None, item_name: str) -> Path | None:
    if not flag:
        logging.info(f"Don't save any {item_name}.")
        return None
    if directory is None:
        raise ValueError(f"{item_name}_path must be provided if {item_name}_flag=True")
    
    output_path = Path(directory)
    logging.info(f"Save {item_name} to {output_path}")

    try:
        output_path.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        logging.error(f"Could not create {item_name} directory '{output_path}': {e}")
        raise
        
    return output_path

def Run_Normalization(
        adata: ad.AnnData,
        sample: str,
        batch: str | None,
        n_neighbors: int = 50,
        n_pcs: int = 50,
        n_top_genes: int = 3000,
        resolution: float = 1,
        do_scale: bool = False,
        fig_flag: bool = False,
        fig_dir: str | None = None
) -> ad.AnnData:
    """
    Perform normalization, feature selection, dimensionality reduction, and clustering 
    on a single-cell RNA-seq AnnData object, with optional visualization outputs.

    Workflow:
    1. Store the raw count matrix in `adata.layers['count']`.
    2. Apply library-size normalization (`sc.pp.normalize_total`, target_sum=1e4) 
       and log1p transformation.
    3. Identify highly variable genes using the Seurat flavor (default: top 3000 genes).
       - Optionally plot and save the highly variable gene distribution.
    4. Run PCA on the highly variable genes.
    5. Construct the neighborhood graph using the PCA representation.
    6. Compute a UMAP embedding.
    7. Perform Leiden clustering, storing results in `adata.obs['leiden']`.
    8. Optionally plot and save UMAP visualizations colored by sample and clustering.

    Parameters
    ----------
    adata : AnnData
        AnnData object containing the single-cell expression matrix.
    sample : str
        Sample name used for figure file naming.
    n_neighbors : int, default=50
        Number of neighbors for graph construction.
    n_pcs : int, default=50
        Number of principal components to use.
    resolution : float, default=0.6
        Resolution parameter for Leiden clustering. Higher values yield more clusters.
    fig_flag : bool, default=False
        Whether to generate and save QC and result figures.
    fig_path : str | None
        Path to save figures. Required if `fig_flag=True`.

    Returns
    -------
    AnnData
        Updated AnnData object containing normalized expression values, 
        PCA/UMAP embeddings, and Leiden clustering results.
    """
    logging.info("begin normalize")
    fig_dir = setup_output_directory(fig_flag, fig_dir, "normalization figures")
    adata.layers['count'] = adata.X
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    # annotation need complete var and X; and the raw.X must be normalize
    # adata.raw = adata
    adata.raw = adata
    sc.pp.highly_variable_genes(adata, 
                            n_top_genes =  n_top_genes,
                            flavor='seurat',
                            subset=False, 
                            batch_key=batch)
    if fig_flag:
        sc.pl.highly_variable_genes(adata)
        high_var_path = fig_dir / f"{sample}_{resolution}_{n_top_genes}-high_variable_genes.png"
        plt.savefig(high_var_path,bbox_inches="tight")
        plt.close()
    
    adata = adata[:, adata.var.highly_variable].copy()

    #scale;if there exist some extrem value in X,below should be added;otherwise not
    if do_scale:
        sc.pp.scale(adata, max_value=10) # the X may be negative after this dispose

    # avoid wrong and speed;didn't affect the obs.cpy():avoid warning
    # adata = adata[:, adata.var.highly_variable].copy()
    
    try:
        sc.pp.pca(adata, mask_var="highly_variable" , n_comps = n_pcs)
    except TypeError:  # 旧版本没有 mask_var 参数
        sc.pp.pca(adata, use_highly_variable=True, n_comps = n_pcs)
    if fig_flag:
        sc.pl.pca_variance_ratio(adata, n_pcs = n_pcs, log=True)
        pca_var_path = fig_dir / f"{sample}_{resolution}_{n_top_genes}-pca_variance_ratio.png"
        plt.savefig(pca_var_path,bbox_inches="tight")
        plt.close()
    sc.pp.neighbors(adata, n_neighbors = n_neighbors, n_pcs = n_pcs ,use_rep = 'X_pca')
    sc.tl.umap(adata) # caculate umap cooordinate for plot
    sc.tl.leiden(adata, resolution=resolution, key_added='leiden',flavor="igraph", directed=False)                                                                                                                                                 
    
    if fig_flag:
        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(15, 5), constrained_layout=True)
        sc.pl.umap(adata, color=['sample'], frameon=False, ax=ax1, show=False)
        sc.pl.umap(adata, color=['leiden'], frameon=False, ax=ax2, show=False)
        umap_normalize_path = fig_dir / f"{sample}_{resolution}_{n_top_genes}-umap-normalization-two.png"
        fig.savefig(umap_normalize_path,bbox_inches="tight")
        plt.close(fig)
    logging.info("end normalize")
    return adata

def Run_batchRemove(
        adata: ad.AnnData,
        sample: str,
        batch: str | list[str],
        n_neighbors: int = 50,
        n_pcs: int = 50,
        resolution: float = 1,
        fig_flag: bool = False,
        fig_dir: str | None = None,
        method: Literal["bbknn", "harmony"] = "harmony"
) -> ad.AnnData:
    """
    Perform batch effect correction, dimensionality reduction, clustering, 
    and visualization on a single-cell AnnData object.

    Workflow:
    1. Run PCA on the expression matrix.
    2. Construct the neighborhood graph (before batch correction).
    3. Apply batch correction using one of the supported methods:
       - "bbknn": Batch Balanced KNN integration (`sc.external.pp.bbknn`).
       - "harmony": Harmony integration (`sc.external.pp.harmony_integrate`).
       - Otherwise: raise an error if unsupported method is specified.
    4. Compute a UMAP embedding of the batch-corrected space.
    5. Perform Leiden clustering, storing results in `adata.obs['leiden']`.
    6. Optionally generate and save UMAP plots:
       - Side-by-side UMAP colored by sample and clustering labels.
       - A single enhanced UMAP plot highlighting Leiden clusters.

    Parameters
    ----------
    adata : AnnData
        AnnData object containing the single-cell expression matrix.
    sample : str
        Sample name used for figure file naming.
    n_neighbors : int, default=50
        Number of neighbors for initial graph construction.
    n_pcs : int, default=50
        Number of principal components to use for PCA.
    resolution : float, default=0.6
        Resolution parameter for Leiden clustering. Higher values produce more clusters.
    fig_flag : bool, default=False
        Whether to generate and save UMAP visualizations.
    fig_dir : str | None
        Directory path for saving figures. Required if `fig_flag=True`.
    method : {"bbknn", "harmony"}, default="harmony"
        Batch correction method to use.

    Returns
    -------
    AnnData
        Updated AnnData object containing batch-corrected embeddings, 
        UMAP layout, and Leiden clustering results.
    """
    logging.info("begin batch effect correction")
    fig_dir = setup_output_directory(fig_flag, fig_dir, "batch correction figures")
    if method == 'bbknn':
        sc.external.pp.bbknn(adata, batch_key = batch)
        # no need to run sc.pp.neighbors
    elif method == 'harmony':
        sc.external.pp.harmony_integrate(adata, key = batch)
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs = n_pcs ,use_rep = 'X_pca_harmony')
    else:
        raise ValueError("This method is not supported for the time being,more batch remove methods need to be developed")
    
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution = resolution, key_added='leiden',flavor="igraph", directed=False)
    if fig_flag:
        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(15, 5), constrained_layout=True)
        sc.pl.umap(adata, color='sample', frameon=False, ax=ax1, show=False)
        sc.pl.umap(adata, color='leiden', frameon=False, ax=ax2, show=False)
        umap_two_path = fig_dir / f"{sample}_{resolution}-umap-{method}-two.png"
        fig.savefig(umap_two_path,bbox_inches="tight")    
        plt.close(fig)
        
        fig, ax = plt.subplots(figsize=(6, 8))
        sc.pl.umap(adata, color='leiden', frameon=False, ax=ax, show=False)
        umap_leiden_path = fig_dir / f"{sample}_{resolution}-umap-{method}-leiden.png"
        fig.savefig(umap_leiden_path,dpi = 300,bbox_inches="tight")
        plt.close(fig)
    logging.info("end batch effect correction")
    return adata

##############################################
def combine_group(
        ad_list: list[ad.AnnData],
        join_method: str = 'inner',
        axis_method: int = 0,
        merge_method: str = 'same'
) -> ad.AnnData:
    """
    Combine multiple AnnData objects into a single AnnData object with flexible options.

    Parameters
    ----------
    ad_list : list[ad.AnnData]
        A list of AnnData objects to be combined.
    join_method : str, optional (default='outer')
        Determines how to handle variables (genes) or observations (cells) that are not shared across all objects:
        - 'outer': keep the union of all obs/var, filling missing values with NaN
        - 'inner': keep only the intersection of obs/var
    axis_method : int, optional (default=0)
        Axis along which to concatenate:
        - 0: concatenate along observations (cells)
        - 1: concatenate along variables (genes)
    merge_method : str, optional (default='unique')
        How to handle duplicate obs/var names:
        - 'unique': automatically make duplicate names unique
        - 'same': require that duplicate names correspond to identical data

    Returns
    -------
    ad.AnnData
        A new AnnData object containing all input objects combined.
        The variable names (genes) and observation names (cells) are guaranteed to be unique.

    Notes
    -----
    - Uses `scanpy.concat` to perform the concatenation with the specified options.
    - After concatenation, `var_names_make_unique()` and `obs_names_make_unique()` are applied to ensure name uniqueness.
    """
    adata = sc.concat(ad_list, axis=axis_method, join=join_method, merge=merge_method)
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    return adata




if __name__ == '__main__':
    start = datetime.datetime.now()
    # parser = argparse.ArgumentParser(description="A script to process cell cycle data.")
    # parser.add_argument('--options', type=str, required=True, help='options to execute procedure')
    # parser.add_argument('--input', type=str, required=True, help='Path to input file')
    # parser.add_argument('--out', type=str, required=True, help='Path to out file')
    h5ad = ""
    samples = ["GBM27","GBM28","GBM29"]
    GBM_SC = []
    GBM_TE = []
    for i in samples:
        adata_SC = sc.read_h5ad(f"{h5ad}/{i}-SC-QC.h5ad")
        adata_TE = sc.read_h5ad(f"{h5ad}/{i}-TE-QC.h5ad")
        #TE may don't need QC,because it already filter after generating by scTE
        # adata_TE = sc.read_h5ad(f"{h5ad}/{i}-TE-raw.h5ad")
        GBM_SC.append(adata_SC)
        GBM_TE.append(adata_TE)
    end = datetime.datetime.now()
    logging.info("程序运行时间："+str((end-start).seconds/3600)+"h")
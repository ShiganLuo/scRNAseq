from pathlib import Path
import numpy as np
import scanpy as sc
import seaborn as sns
from scipy.stats import median_abs_deviation
import anndata as ad
import matplotlib.pyplot as plt
import scrublet as scr
from math import sqrt
import pandas as pd
import datetime
import umap
import logging
import matplotlib as mpl
logging.basicConfig(

  level=logging.INFO,

  format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',

  datefmt='%Y-%m-%d %H:%M:%S'

)
sc.settings.verbosity = 0
sc.settings.set_figure_params(
    dpi=80,
    facecolor="white",
    frameon=False,
)
###python version >= 3.9
def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier

def slice_anndata(adata, obs_mask=None, var_mask=None):
    """
    Slice an AnnData object and automatically update all layers to match the new shape.
    Supports both pandas Series and numpy arrays for obs_mask and var_mask.
    """

    # Default: select all columns if var_mask is None
    if var_mask is None:
        var_mask = slice(None) # slice == :
    else:
        var_mask = np.asarray(var_mask,dtype=bool)
    if obs_mask is None:
        obs_mask = slice(None)
    else:
        obs_mask = np.asarray(obs_mask, dtype=bool)

    # adata[["Cell_1", "Cell_10"], ["Gene_5", "Gene_1900"]]
    adata_new = adata[obs_mask, var_mask].copy()

    # Synchronize layers
    for layer_name, mat in adata.layers.items():
        new_mat = mat[obs_mask, var_mask].copy()
        adata_new.layers[layer_name] = new_mat

    return adata_new

def lowquality(
    adata: ad.AnnData,
    sample: str,
    low_count_outlier: int = 5,
    low_gene_outlier: int = 5,
    top20_outlier: int = 5,
    mt_outlier: int | None  = 3,
    mt_cutoff: int | None = 8,
    n_genes_by_counts_cutoff: tuple[int,int] | None = (200,9999999),
    total_counts_cutoff: tuple[int,int] | None = (500,9999999),
    fig_flag: bool = False,
    fig_dir: str | None = None
) -> ad.AnnData:
    """
    Filter low-quality cells in single-cell RNA-seq data based on multiple quality control (QC) metrics,
    and optionally generate diagnostic plots.

    Parameters
    ----------
    adata : AnnData
        The input single-cell AnnData object containing expression matrix and metadata.
    sample : str
        Sample name used for labeling and saving plots.
    low_count_outlier : int, default=5
        Number of median absolute deviations (MADs) for identifying outliers in
        ``log1p_total_counts`` (library size).
    low_gene_outlier : int, default=5
        Number of MADs for identifying outliers in ``log1p_n_genes_by_counts`` (gene complexity).
    top20_outlier : int, default=5
        Number of MADs for identifying outliers in ``pct_counts_in_top_20_genes`` (dominance of top genes).
    mt_outlier : int, default=3
        Number of MADs for identifying outliers in ``pct_counts_mt`` (mitochondrial fraction).
        Only used if ``use_mt_outlier=True``.
    mt_cutoff : int, default=8
        Fixed cutoff threshold for mitochondrial percentage (``pct_counts_mt > mt_cutoff``).
        Cells above this threshold are always considered outliers.
    use_mt_outlier : bool, default=True
        Whether to apply MAD-based outlier detection on ``pct_counts_mt``.
        If ``False``, mitochondrial filtering relies solely on ``mt_cutoff``.
    fig_dir : str or None, optional
        Directory path to save QC plots. Required if ``fig_flag=True``.
    fig_flag : bool, default=False
        Whether to generate and save diagnostic plots.

    Returns
    -------
    AnnData
        A filtered AnnData object with low-quality cells removed.

    Raises
    ------
    ValueError
        If ``fig_flag=True`` but ``fig_dir`` is not provided.

    Notes
    -----
    - QC metrics include mitochondrial (``MT-``), ribosomal (``RPS``/``RPL``), and
      hemoglobin (``^HB``) gene fractions, calculated using
      :func:`scanpy.pp.calculate_qc_metrics`.
    - Outliers are detected using median absolute deviation (MAD) and/or fixed thresholds.
    - Cells are flagged as low quality if:
        * They are MAD outliers for total counts, gene complexity, or top20 dominance.
        * They are mitochondrial outliers, as defined by ``mt_outlier`` (if enabled) 
          or ``mt_cutoff``.
    - If ``fig_flag=True``, histograms, violin plots, and scatter plots of QC metrics
      are saved before and after filtering.
    """
    logging.info(f"{sample}'s low quality cells filtering start")
    fig_dir = setup_output_directory(fig_flag, fig_dir, "figures")

    logging.info(f"Number of cells originally: {adata.n_obs}")
    ### metrix caculate
    # human mitcondronial gene most start with MT-, but mouse's most start with mt-
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.upper().str.startswith(("RPS", "RPL"))
    adata.var["hb"] = adata.var_names.str.upper().str.match(r"^HB[^P]")

    sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True)


    if n_genes_by_counts_cutoff is not None:
        logging.info(f"n_genes_by_counts_cutoff is open: {n_genes_by_counts_cutoff}")
        low, high = n_genes_by_counts_cutoff
        before = adata.n_obs
        adata = adata[(adata.obs["n_genes_by_counts"] > low) &
                    (adata.obs["n_genes_by_counts"] < high)].copy()
        logging.info(f"Filtered by n_genes_by_counts: {before} → {adata.n_obs}")
    if total_counts_cutoff is not None:
        logging.info(f"total_counts_cutoff is open: {total_counts_cutoff}")
        low, high = total_counts_cutoff
        before = adata.n_obs
        adata = adata[(adata.obs["total_counts"] > low) &
                    (adata.obs["total_counts"] < high)].copy()
        logging.info(f"Filtered by total_counts: {before} → {adata.n_obs}")

    logging.info(f"Number of cells after hard filtered: {adata.n_obs}")
    adata.obs["outlier"] = (
    is_outlier(adata, "log1p_total_counts", low_count_outlier)
    | is_outlier(adata, "log1p_n_genes_by_counts", low_gene_outlier)
    | is_outlier(adata, "pct_counts_in_top_20_genes", top20_outlier)
    )
    
    if mt_outlier is not None and mt_cutoff is not None:
        adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", mt_outlier) | (
        adata.obs["pct_counts_mt"] > mt_cutoff
        )
    elif mt_cutoff is not None:
        adata.obs["mt_outlier"] = adata.obs["pct_counts_mt"] > mt_cutoff
    elif mt_outlier is not None:
        adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", mt_outlier)
    else:
        pass
    

    # adata_filter = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()
    obs_mask = ((~adata.obs.outlier) & (~adata.obs.mt_outlier))
    adata_filtered = slice_anndata(adata,obs_mask)
    logging.info(f"Number of cells after filtering of low quality cells: {adata_filtered.n_obs}")

    
    if fig_flag:
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        sns.histplot(adata.obs["total_counts"], bins=100, kde=False,ax=axes[0])
        axes[0].set_title("Before filtering")
        sns.histplot(adata_filtered.obs["total_counts"], bins=100, kde=False,ax=axes[1])
        axes[1].set_title("After filtering")
        total_counts_path = fig_dir / f"{sample}-total_counts.png"
        fig.savefig(total_counts_path,dpi=300, bbox_inches='tight')
        plt.close(fig)

        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        sc.pl.violin(adata, "pct_counts_mt",show=False,ax=axes[0])
        axes[0].set_title("Before filtering")
        sc.pl.violin(adata_filtered, "pct_counts_mt",show=False,ax=axes[1])
        axes[1].set_title("After filtering")
        pct_counts_mt_path = fig_dir / f"{sample}-pct_counts_mt.png"
        plt.savefig(pct_counts_mt_path,dpi=300, bbox_inches='tight')
        plt.close(fig)

        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt",show=False,ax=axes[0])
        axes[0].set_title("Before filtering")
        sc.pl.scatter(adata_filtered, "total_counts", "n_genes_by_counts", color="pct_counts_mt",show=False,ax=axes[1])
        axes[1].set_title("After filtering")
        for ax in fig.axes[::-1]:
            if ax not in list(axes):
                fig.delaxes(ax)
        try:
            coll = axes[1].collections[0]            # PathCollection
            cmap = coll.get_cmap()
            norm = coll.norm
        except Exception:
            cmap = None
            # 取数据范围作为 norm 回退方案
            vmin = adata.obs["pct_counts_mt"].min()
            vmax = adata.obs["pct_counts_mt"].max()
            norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
        sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=axes, orientation="vertical", fraction=0.03, pad=0.04)
        cbar.set_label("pct_counts_mt (%)")
        tnp_path = fig_dir / f"{sample}-scatter.png"
        plt.savefig(tnp_path, dpi=300, bbox_inches='tight')
        plt.close(fig)
    logging.info(f"{sample}'s low quality cells filtering completed")
    return adata_filtered 


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

def Doublet_scrub(
        adata: ad.AnnData,
        sample: str,
        scrub_estimated_rate: float = 0.06,
        sim_doublet_ratio: int = 2,
        scrub_min_counts: int = 2,
        scrub_min_cells: int = 3,
        scrub_min_gene_variability_pctl: int = 85,
        scrub_n_prin_comps: int = 30,
        fig_flag: bool = False,
        fig_dir: str | None = None,
        table_flag: bool = False,
        table_dir: str | None = None
) -> ad.AnnData:
    """
    Detect and filter potential doublets in single-cell RNA-seq data using Scrublet,
    automatically adjusting the expected doublet rate based on the initial estimation.
    
    This function first runs Scrublet with a default expected doublet rate (0.06) 
    to obtain an initial estimate of the doublet fraction. It then reruns Scrublet 
    using the estimated doublet rate to improve detection sensitivity, especially 
    for samples with higher-than-expected doublet rates.
    
    Parameters
    ----------
    adata : AnnData
        Input AnnData object containing single-cell gene expression counts.
    sample : str
        Sample name, used for labeling plots or tables.
    sim_doublet_ratio : int, default=2
        Ratio of simulated doublets to observed cells.
    scrub_min_counts : int, default=2
        Minimum number of counts required for a gene to be considered.
    scrub_min_cells : int, default=3
        Minimum number of cells a gene must be expressed in to be considered.
    scrub_min_gene_variability_pctl : int, default=85
        Percentile of gene variability used to filter highly variable genes.
    scrub_n_prin_comps : int, default=30
        Number of principal components used for doublet simulation and detection.
    fig_flag : bool, default=False
        Whether to generate diagnostic plots (histogram and TSNE embedding).
    fig_dir : str or None, optional
        Directory to save diagnostic plots. Required if fig_flag=True.
    table_flag : bool, default=False
        Whether to save doublet scores and predictions as a table.
    table_dir : str or None, optional
        Directory to save output table. Required if table_flag=True.
    
    Returns
    -------
    adata_filtered : AnnData
        Filtered AnnData object with predicted doublets removed.
        Two new columns are added to adata.obs:
        - 'doublet_scores': doublet likelihood scores for each cell.
        - 'predicted_doublets': boolean labels indicating predicted doublets.
    
    Notes
    -----
    - The function uses a two-step Scrublet approach:
        1. Initial estimation using default expected rate.
        2. Recalculation using the estimated doublet rate for more accurate detection.
    - Diagnostic plots help visualize doublet score distributions and predicted doublets.
    - This approach improves filtering accuracy for samples with high doublet rates.
    """
    logging.info(f"Begin filtering doublets for sample {sample}")
    fig_dir = setup_output_directory(fig_flag, fig_dir, "figures")
    table_dir = setup_output_directory(table_flag, table_dir, "tables")

    counts_matrix = adata.to_df()
    n_cells = adata.shape[0]
    logging.info(f"Number of cells originally: {adata.n_obs}")
    # initial Scrublet run with default expected rate (0.06)
    logging.info(f"initial Scrublet run with the expected rate {scrub_estimated_rate}")
    scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=scrub_estimated_rate,
                              n_neighbors=round(0.5 * sqrt(n_cells)),
                              sim_doublet_ratio=sim_doublet_ratio)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(
        min_counts=scrub_min_counts,
        min_cells=scrub_min_cells,
        min_gene_variability_pctl=scrub_min_gene_variability_pctl,
        n_prin_comps=scrub_n_prin_comps
    )

    if fig_flag:
        logging.info(f"begin plot histogram and umap for initial Scrublet run")
        scrub.plot_histogram() # origin scrublet use numpy 1.x()
        histogram_path = fig_dir / f"{sample}-initial-histogram.png"
        plt.savefig(histogram_path, dpi=300, bbox_inches='tight')
        plt.close()

        X_umap = umap.UMAP(n_neighbors=15, min_dist=0.5, random_state=0).fit_transform(scrub.manifold_obs_)
        scrub.set_embedding('UMAP', X_umap)
        scrub.plot_embedding('UMAP', order_points=True)
        umap_path = fig_dir / f"{sample}-initial-umap.png"
        plt.savefig(umap_path,dpi=300, bbox_inches='tight')
        plt.close()

    # Save table if requested
    out_df = pd.DataFrame({
        'barcodes': counts_matrix.index,
        'doublet_scores': doublet_scores,
        'predicted_doublets': predicted_doublets
    })
    if table_flag:
        logging.info(f"save table to {table_dir}/{sample}-doublet.txt")
        doublet_table_path = table_dir / f"{sample}-doublet.txt"
        out_df.to_csv(doublet_table_path, index=False, header=True)

    # Add to AnnData and filter doublets
    adata.obs['doublet_scores'] = doublet_scores
    adata.obs['predicted_doublets'] = predicted_doublets
    obs_mask = ~adata.obs['predicted_doublets']
    # adata_filtered = adata[~adata.obs['predicted_doublets']]
    adata_filtered = slice_anndata(adata,obs_mask)
    logging.info(f"Number of cells after filtering of low quality cells: {adata_filtered.n_obs}")
    logging.info("Doublet filtering completed.")
    return adata_filtered

def Doublet_scrub_scanpy_official(
    adata: ad.AnnData,
    sample: str,
    scrub_estimated_rate: float = 0.06,
    fig_flag: bool = False,
    fig_dir: str | None = None,
    table_flag: bool = False,
    table_dir: str | None = None,
    **kwargs
) -> ad.AnnData:
    """
    Detect and filter doublets using sc.pp.scrublet for maximum Scanpy integration.
    
    Parameters
    ----------
    adata : AnnData
        Input AnnData object containing single-cell gene expression counts.
    sample : str
        Sample name, used for labeling plots or tables.
    scrub_estimated_rate : float, default=0.06
        Expected doublet rate for Scrublet estimation.
    fig_flag : bool, default=False
        Whether to generate diagnostic plots (histogram and UMAP embedding).
    fig_dir : str or None, optional
        Directory to save diagnostic plots. Required if fig_flag=True.
    table_flag : bool, default=False
        Whether to save doublet scores and predictions as a table.
    table_dir : str or None, optional
        Directory to save output table. Required if table_flag=True.
    **kwargs : 
        Additional keyword arguments passed to sc.pp.scrublet.
        
    Returns
    -------
    adata_filtered : AnnData
        Filtered AnnData object with predicted doublets removed.
    """
    logging.info(f"Begin filtering doublets for sample {sample} using sc.pp.scrublet.")
    
    fig_dir_path = setup_output_directory(fig_flag, fig_dir, "figures")
    table_dir_path = setup_output_directory(table_flag, table_dir, "tables")
    
    sc.pp.scrublet(
        adata,
        expected_doublet_rate=scrub_estimated_rate,
        **kwargs 
    )

    if fig_flag and fig_dir_path:
        logging.info("Generating Scrublet diagnostic plots.")
        
        sc.pl.scrublet_score_distribution(
            adata,
            show=False
        )
        histogram_path = fig_dir_path / f"{sample}-doublet-histogram.png"
        plt.suptitle(f"{sample} - Doublet Score Distribution", y=1.02)
        plt.savefig(histogram_path, dpi=300, bbox_inches='tight')
        plt.close()        
    if table_flag and table_dir_path:
        logging.info(f"Saving doublet scores table to {table_dir_path}/{sample}-doublet.csv")
        out_df = adata.obs[['doublet_scores', 'predicted_doublets']].copy()
        out_df.index.name = 'barcodes'
        doublet_table_path = table_dir_path / f"{sample}-doublet.csv"
        out_df.to_csv(doublet_table_path) 

    n_total = adata.shape[0]
    n_doublets = np.sum(adata.obs['predicted_doublets'])
    
    adata_filtered = adata[~adata.obs['predicted_doublets'], :].copy() 
    
    logging.info(f"Doublet filtering completed. Removed {n_doublets} doublets ({n_doublets/n_total:.2%}).")
    
    return adata_filtered

def Doublet_Solo():
    logging.info("This is a placeholder for the Doublet_Solo function.")


if __name__ == '__main__':
    start = datetime.datetime.now()
    # parser = argparse.ArgumentParser(description="A script to process cell cycle data.")
    # parser.add_argument('--options', type=str, required=True, help='options to execute procedure')
    # parser.add_argument('--input', type=str, required=True, help='Path to input file')
    # parser.add_argument('--out', type=str, required=True, help='Path to out file')
    h5ad = "/home/lsg/Data/glioblastoma/output/new/h5ad"
    samples = ["GBM27","GBM28","GBM29"]
    for i in samples:
        adata_SC = sc.read_h5ad(f"{h5ad}/{i}-SC-raw.h5ad")
        adata_TE = sc.read_h5ad(f"{h5ad}/{i}-TE-raw.h5ad")
        adata_SC = lowquality(adata_SC,f"{i}-SC")
        adata_TE = lowquality(adata_TE,f"{i}-TE")
        adata = lowquality(adata,f"{i}-SC",fig_flag = True,fig_dir = "/disk5/luosg/scRNAseq/output/result/QC")
        adata = Doublet_scrub(adata,f"{i}-TE",fig_flag = True,fig_dir = "/disk5/luosg/scRNAseq/output/result/doublet/plot",
                                table_flag=True,table_dir="/disk5/luosg/scRNAseq/output/result/doublet/table")
        adata_SC.write(f'{h5ad}/{f"{i}-SC"}-QC.h5ad')
        adata_TE.write(f'{h5ad}/{f"{i}-TE"}-QC.h5ad')
    end = datetime.datetime.now()
    print("程序运行时间："+str((end-start).seconds/3600)+"h")
    
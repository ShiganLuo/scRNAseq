import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pertpy as pt
import scanpy as sc
import seaborn as sns
import logging
import anndata as ad
import os
from pathlib import Path
os.environ["QT_QPA_PLATFORM"] = "offscreen"
from schist.inference import fit_model
logging.basicConfig(
	level=logging.INFO,
	format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
	datefmt='%Y-%m-%d %H:%M:%S'
)

from typing import Optional, List

def run_scCODA(
    adata: ad.AnnData,
    fig_prefix: str,
    cell_type_identifier: str = "cell_label",
    sample_identifier: str = "batch",
    covariate_obs: Optional[List[str]] = None,
    feature_name: str = "condition",
    reference_cell_type: str = "automatic",
    formulas: List[str] = ["condition"],
    modality_key: str = "coda",
    rng_key: int = 1234
):
    """
    Run a complete scCODA analysis pipeline on a single-cell dataset.
    
    This function performs:
    1. Preparation of a scCODA-compatible data object from AnnData.
    2. Visualization of cell type abundance using boxplots and stacked barplots.
    3. Preparation of the dataset with covariates and reference cell type.
    4. Bayesian inference using the NUTS sampler to identify compositional changes.
    5. Plotting of inferred effects as barplots.

    Parameters
    ----------
    adata : AnnData
        Input single-cell dataset containing cell type annotations and covariates.
    fig_prefix : str
        Prefix for saved figure filenames.
    cell_type_identifier : str, default "cell_label"
        Column in `adata.obs` that contains cell type labels.
    sample_identifier : str, default "batch"
        Column in `adata.obs` that identifies the sample or batch.
    covariate_obs : list of str, optional, default ["condition"]
        Covariate columns in `adata.obs` to include in the model.
    feature_name : str, default "condition"
        The feature to visualize and model for compositional changes.
    reference_cell_type : str, default "automatic"
        Reference cell type for scCODA analysis. Can be "automatic" or a specific cell type.
    formula : str, default "condition"
        Formula describing covariates for the Bayesian model (R-style),such as "C(condition, Treatment('CKO'))".
    modality_key : str, default "coda"
        Key for storing results in the AnnData layers.
    rng_key : int, default 1234
        Random seed for reproducibility.

    Returns
    -------
    sccoda_model : pt.tl.Sccoda
        The fitted scCODA model object.
    sccoda_data : pt.tl.CompositionalData
        The processed data object ready for further analysis.
    """

    n_cell_types =  len(adata.obs[cell_type_identifier].unique())
    n_condition = len(adata.obs[feature_name].unique())
    logging.info(f"cell types: {n_cell_types}")
    logging.info(f"condition types: {n_condition}")
    boxplot_figsize = (max(8, n_cell_types * 1), 5)
    barplot_figsize = (max(n_condition * 1,4),max(2, n_cell_types / 4))
    logging.info(f"boxplot_figsize is {boxplot_figsize}")
    logging.info(f"barplot_figsize is {barplot_figsize}")
        # Set default covariates if not provided
    if covariate_obs is None:
        covariate_obs = [feature_name]

    ### Initialize scCODA model
    sccoda_model = pt.tl.Sccoda()

    ### Load AnnData into scCODA data structure
    sccoda_data = sccoda_model.load(
        adata,
        type="cell_level",
        generate_sample_level=True,
        cell_type_identifier=cell_type_identifier,
        sample_identifier=sample_identifier,
        covariate_obs=covariate_obs,
    )

    ### Boxplot visualization of cell type composition
    sccoda_model.plot_boxplots(
        sccoda_data,
        modality_key=modality_key,
        feature_name=feature_name,
        figsize=boxplot_figsize,
        add_dots=True,
        args_swarmplot={},  # Can add custom seaborn parameters
    )
    output_path = Path(f"{fig_prefix}_description_boxplot.png")
    output_path.parent.mkdir(parents=True, exist_ok=True) 
    plt.tight_layout()
    plt.savefig(str(output_path), dpi=300,bbox_inches="tight")
    plt.close()

    ### Stacked barplot visualization of composition
    fig = sccoda_model.plot_stacked_barplot(
        sccoda_data,
        modality_key=modality_key,
        feature_name=feature_name,
        figsize=barplot_figsize,
        return_fig=True
    )
    # matplotlib Figure
    for ax in fig.axes:  # 遍历所有子图
        handles, labels = ax.get_legend_handles_labels()
        if handles:
            ax.legend(
                handles=handles,
                labels=labels,
                title=None,
                fontsize=8,          # 图例字体大小
                title_fontsize=10,   # 图例标题字体
                loc='center left',   # 图例位置
                bbox_to_anchor=(1, 0.5)
            )
    output_path = Path(f"{fig_prefix}_stacked_barplot.png")
    output_path.parent.mkdir(parents=True, exist_ok=True) 
    fig.savefig(str(output_path), dpi=300, bbox_inches="tight")
    plt.close()

    #########################################
    cell_type_order = sorted(adata.obs[cell_type_identifier].unique().tolist())
    for formula in formulas:
        logging.info(f"dispose {formula}")
        if "_vs_" in formula:
            ctrl = formula.split("_vs_")[1]
            expr = formula.split("_vs_")[0]
            if not set([ctrl, expr]).issubset(set(adata.obs[feature_name])):
                logging.info(f"Skipping {formula}: not all groups present in adata")
                continue
            formula_r = f"C({feature_name}, Treatment('{ctrl}'))"
            adata_filtered = adata[adata.obs[feature_name].isin([ctrl, expr])].copy()
            if adata_filtered.n_obs == 0:
                logging.info(f"Skipping {expr} vs {ctrl}: no cells found.")
                continue
        
            # cell_types_in_groups = adata_filtered.obs.groupby(feature_name)[cell_type_identifier].unique()
            # common_cell_types = set(cell_types_in_groups[ctrl]).intersection(cell_types_in_groups[expr])
            # adata_filtered = adata_filtered[adata_filtered.obs[cell_type_identifier].isin(common_cell_types)].copy()
            # if adata_filtered.n_obs == 0:
            #     logging.info(f"Skipping {expr} vs {ctrl}: no cells remain after filtering uncommon cell types.")
            #     continue

            sccoda_data = sccoda_model.load(
                adata_filtered,
                type="cell_level",
                generate_sample_level=True,
                cell_type_identifier=cell_type_identifier,
                sample_identifier=sample_identifier,
                covariate_obs=covariate_obs,
            )
            ### Prepare data for Bayesian inference
            sccoda_data = sccoda_model.prepare(
                sccoda_data,
                modality_key=modality_key,
                formula=formula_r,
                reference_cell_type=reference_cell_type,
            )

            ### Run NUTS sampler for Bayesian inference
            sccoda_model.run_nuts(
                sccoda_data,
                modality_key=modality_key,
                rng_key=rng_key
            )
            effect_df = sccoda_model.credible_effects(sccoda_data, modality_key="coda")
            credible_effects = effect_df[effect_df]
            if credible_effects.empty:
                logging.info(f"No credible effects found for {formula}, skipping plot.")
                continue

            # Plot effects barplot and get the FacetGrid
            g = sccoda_model.plot_effects_barplot(sccoda_data, return_fig=True)
            for ax in g.axes.flatten():  # g.axes is ndarray
                if ax is not None:
                    present_cell_types = [tick.get_text() for tick in ax.get_xticklabels()]
                    new_labels = [ct if ct in present_cell_types else "" for ct in cell_type_order]
                    ax.set_xticks(range(len(cell_type_order)))
                    ax.set_xticklabels(new_labels, rotation=60, ha='right')
            output_path = Path(f"{fig_prefix}_{formula}_effects_barplot.png")
            output_path.parent.mkdir(parents=True, exist_ok=True) 
            g.fig.savefig(str(output_path), dpi=300, bbox_inches="tight")
            plt.close(g.fig)
        else:
            # "C(experiment, Treatment('CKO'))"
            sccoda_data = sccoda_model.load(
                adata,
                type="cell_level",
                generate_sample_level=True,
                cell_type_identifier=cell_type_identifier,
                sample_identifier=sample_identifier,
                covariate_obs=covariate_obs,
            )
            ### Prepare data for Bayesian inference
            sccoda_data = sccoda_model.prepare(
                sccoda_data,
                modality_key=modality_key,
                formula=formula,
                reference_cell_type=reference_cell_type,
            )

            ### Run NUTS sampler for Bayesian inference
            sccoda_model.run_nuts(
                sccoda_data,
                modality_key=modality_key,
                rng_key=rng_key
            )
            
            # Plot effects barplot and get the FacetGrid
            g = sccoda_model.plot_effects_barplot(sccoda_data, return_fig=True)
            for ax in g.axes.flatten():   # g.axes is ndarray
                if ax is not None:
                    plt.setp(ax.get_xticklabels(), rotation=60, ha='right')
            output_path = Path(f"{fig_prefix}_{formula}_effects_barplot.png")
            output_path.parent.mkdir(parents=True, exist_ok=True) 
            g.fig.savefig(str(output_path), dpi=300, bbox_inches="tight")
            plt.close(g.fig)        
    logging.info(f"scCODA analysis completed. Figures saved with prefix: {fig_prefix}")

    return sccoda_model, sccoda_data


    






if __name__ == "__main__":
    adata = sc.read_h5ad("/disk5/luosg/scRNAseq/output/result1/Intestine_annotate.h5ad")
    run_scCODA(adata,"/disk5/luosg/scRNAseq/output/result1/Intestine/Compositional/Intestine",
               cell_type_identifier="celltype",
               sample_identifier="sample",
               covariate_obs=["experiment"],
               feature_name="experiment",
               formulas=["CKO_vs_WT","E2_vs_CKO","TRA_vs_CKO"])
    adata = sc.read_h5ad("/disk5/luosg/scRNAseq/output/result1/Lung_annotate.h5ad")

    run_scCODA(adata,"/disk5/luosg/scRNAseq/output/result1/Lung/Compositional/Lung",
               cell_type_identifier="celltype",
               sample_identifier="sample",
               covariate_obs=["experiment"],
               feature_name="experiment",
               formulas=["CKO_vs_WT","E2_vs_CKO","TRA_vs_CKO"])
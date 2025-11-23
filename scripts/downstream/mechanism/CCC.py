from pdb import run
from numpy import outer
import scanpy as sc
import liana as li
import anndata as ad
from liana.method import cellphonedb
from liana.method import rank_aggregate
from typing import Callable
import warnings
from pathlib import Path
warnings.filterwarnings("ignore")

def run_liana_one_condition(
        adata:ad.AnnData,
        outfile:str,
        method:str="cellphonedb",
        condition_key:str="condition",
        target_condition:str="stim",
        cell_type_key:str="celltype",
    ) -> None:
    """
    Run LIANA on your data to infer cell-cell communication.
    Parameters
    ----------
    - adata : ad.AnnData
        Annotated data matrix.
    - method : str, optional
        The method to use for cell-cell communication inference. Default is "cellphonedb".    
    Returns
        -------
        None
    """
    adata_stim = adata[adata.obs[condition_key] == target_condition].copy()
    cellphonedb(
        adata_stim, groupby=cell_type_key, use_raw=True, return_all_lrs=True, verbose=False
    )
    adata_stim.uns["liana_res"].drop_duplicates(
        ["ligand_complex", "receptor_complex"]
    ).to_csv(outfile,index=False,sep="\t")
    return adata_stim

def liana_rank_aggregate(
        adata:ad.AnnData,
        outfile:str,
        condition_key:str="condition",
        target_condition:str="stim",
        cell_type_key:str="celltype",
    ) -> None:
    """
    Aggregate LIANA results using rank aggregation.
    Parameters
    ----------
    - adata : ad.AnnData
        Annotated data matrix.
    - cell_type_key : str, optional
        The key in adata.obs that contains cell type information. Default is "celltype".    
    Returns
        -------
        None
    """
    adata_stim = adata[adata.obs[condition_key] == target_condition].copy()
    rank_aggregate(
        adata_stim, groupby=cell_type_key, return_all_lrs=True, use_raw=True, verbose=False
    )
    adata_stim.uns["liana_res"].drop_duplicates(
        ["ligand_complex", "receptor_complex"]
    ).to_csv(outfile,index=False,sep="\t")
    return adata_stim

def filter_fun_single(row):
    return row["cellphone_pvals"] < 0.01

def filter_fun_aggregate(row):
    return row["specificity_rank"] < 0.05

def dotplot_example(
        adata_stim:ad.AnnData,
        source_cell_types:list,
        target_cell_types:list,
        outImage:str,
        filter_fun:Callable[[any], bool] = filter_fun_single,
        color:str="lr_means",
        size:str="cellphone_pvals",
        
) -> None:

    fig = li.pl.dotplot(
        adata=adata_stim,
        colour=color,
        size=size,
        inverse_size=True,  # we inverse sign since we want small p-values to have large sizes
        # We choose only the cell types which we wish to plot
        source_labels=source_cell_types,
        target_labels=target_cell_types,
        filter_fun = filter_fun,
        # as this type of methods tends to result in large numbers
        # of predictions, we can also further order according to
        # expression magnitude
        orderby="lr_means",
        orderby_ascending=False,  # we want to prioritize those with highest expression
        top_n=20,  # and we want to keep only the top 20 interactions
        figure_size=(9, 5),
        size_range=(1, 6),
    )
    fig.save(outImage, width=8, height=6, dpi=300)


def run_steady_state_CCC(
        adata:ad.AnnData,
        outdir:str,
        source_cell_types:list,
        target_cell_types:list,
):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    outdir = str(outdir)
    for experiment in adata.obs['experiment'].unique():
        outfile = f"{outdir}/liana_dotplot_{experiment}.tsv"
        adata_stim = run_liana_one_condition(
            adata,
            outfile=outfile,
            method="cellphonedb",
            condition_key="experiment",
            target_condition=experiment,
            cell_type_key="celltype",
        )
        # Filter paired cell types based on their presence in the current experiment
        liana_res = adata_stim.uns['liana_res']
        filtered_pairs = liana_res[
            (liana_res['source'].isin(source_cell_types)) &
            (liana_res['target'].isin(target_cell_types))
        ][['source', 'target']].drop_duplicates()
        filtered_source = filtered_pairs['source'].unique().tolist()
        filtered_target = filtered_pairs['target'].unique().tolist()
        print("Filtered source:", filtered_source)
        print("Filtered target:", filtered_target)
        try:
            outeImage = f"{outdir}/liana_dotplot_{experiment}.png"
            dotplot_example(
                adata_stim,
                filtered_source,
                filtered_target,
                outeImage,
            )
        except Exception as e:
            print(f"Failed to create dotplot for {experiment}: {e}")
            continue

        outfile_aggregate = f"{outdir}/liana_dotplot_{experiment}_aggregate.tsv"
        adata_stim_aggregate = liana_rank_aggregate(
            adata,
            outfile=outfile_aggregate,
            condition_key="experiment",
            target_condition=experiment,
            cell_type_key = "celltype"
        )
        liana_res = adata_stim.uns['liana_res']
        filtered_pairs = liana_res[
            (liana_res['source'].isin(source_cell_types)) &
            (liana_res['target'].isin(target_cell_types))
        ][['source', 'target']].drop_duplicates()
        filtered_source = filtered_pairs['source'].unique().tolist()
        filtered_target = filtered_pairs['target'].unique().tolist()
        print("Filtered source:", filtered_source)
        print("Filtered target:", filtered_target)
        try:
            outeImage = f"{outdir}/liana_dotplot_{experiment}_aggregate.png"
            dotplot_example(
                adata_stim_aggregate,
                filtered_source,
                filtered_target,
                outeImage,
                filter_fun_aggregate,
                color = "magnitude_rank",
                size = "specificity_rank"
            )
        except Exception as e:
            print(f"Failed to create dotplot for {experiment} aggregate: {e}")
            continue

def upperGeneNames(
        adata:ad.AnnData
):
    """
    Example to convert gene names to uppercase in AnnData and its raw attribute.
    Also demonstrates how to run LIANA and generate dotplots for specific cell types.
    """
    adata.var_names = adata.var_names.str.upper()
    # adata raw is read only with gene names in lowercase, need to convert to uppercase
    raw_X = adata.raw.X.copy()
    raw_var = adata.raw.var.copy()
    raw_var.index = raw_var.index.str.upper()
    adata.raw = ad.AnnData(
        X=raw_X,
        var=raw_var
    )
    return adata
if __name__ == "__main__":
    # adata_Intestine = sc.read_h5ad("/disk5/luosg/scRNAseq/output/combine/Intestine/Intestine_annotate.h5ad")
    # adata_Intestine = upperGeneNames(adata_Intestine)
    # run_steady_state_CCC(
    #     adata_Intestine,
    #     outdir="/disk5/luosg/scRNAseq/output/combine/Intestine/mechanisms/CCC",
    #     source_cell_types = ["B cell", "T cytotoxic cell","Macrophage"],
    #     target_cell_types = ["Paneth", "Goblet", "Enterocyte cell","Enteroendocrine cell","Intestinal stem cell","Tuft"]
    # )

    adata_Lung = sc.read_h5ad("/disk5/luosg/scRNAseq/output/combine/Lung/Lung_annotate.h5ad")
    adata_Lung = upperGeneNames(adata_Lung)
    run_steady_state_CCC(
        adata_Lung,
        outdir="/disk5/luosg/scRNAseq/output/combine/Lung/mechanisms/CCC/T_cytotoxic_cell_Il7r",
        source_cell_types = ["T cytotoxic cell(Il7r)"],
        target_cell_types = ["T cytotoxic cell","T cytotoxic cell(Il7r)","Macrophage","Neutrophil"]
    )
    # source_cell_types = ["Macrophage","Neutrophil","T cytotoxic cell(Il7r)","B cell","T cytotoxic cell","Fibroblast","Macrophage(Mki67)"],
    # target_cell_types = ["endothelial cells","smooth muscle cells","Fibroblast",
    #                          "Airway goblet cells","Pulmonary alveolar type II cells",
    #                          "Lymphatic endothelial cell","Pulmonary alveolar type I cells"]



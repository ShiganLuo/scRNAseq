import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import logging
import anndata as ad
import scanpy as sc
from functools import reduce
from typing import Optional, TypedDict

logging.basicConfig(

    level=logging.INFO,

    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',

    datefmt='%Y-%m-%d %H:%M:%S'
)

def TEfamily(
    DEG_file: str,
    rmsk_file: str,
    out: str
):
    """
    Analyze and visualize transposable element (TE) family regulation 
    based on differential expression results.

    This function performs the following steps:
    1. Filter TEs in the DEG file using the RepeatMasker (rmsk) annotation file.
    2. Select significant TEs based on thresholds: |logFC| > 0.58 and -log10(PValue) > 1.3.
    3. Classify significant TEs into 'Up' or 'Down' regulated groups depending on logFC.
    4. Count the number of significant TEs per TE family and regulation direction.
    5. Save the results as a tab-delimited file (<out>_TEfamily.tsv).
    6. Generate and save a bar plot (<out>_TEfamily.png) comparing upregulated vs 
       downregulated TE families, with counts annotated.

    Parameters
    ----------
    DEG_file : str
        Path to the differential expression results file containing logFC and PValue.
    rmsk_file : str
        Path to the RepeatMasker annotation file used for filtering TE entries.
    out : str
        Prefix for output files (.tsv and .png).

    Raises
    ------
    ValueError
        If no significant TE entries are found after filtering.

    Outputs
    -------
    - A TSV file containing counts of upregulated and downregulated TE families.
    - A bar plot PNG file visualizing TE family regulation.
    """ 
    logging.info("filter TE according rmsk file from DEG file")
    new_df = filterTE(DEG_file,rmsk_file)
    # 筛选显著
    logging.info("filter meaningful TE")
    significant_df = new_df[
        (abs(new_df["logFC"]) > 0.58) & 
        (-np.log10(new_df["PValue"]) > 1.3)
    ].copy()
    
    if significant_df.empty:
        raise ValueError("The filtered DataFrame of significant TEs is empty.")
    
    # 直接分组计数：同时保留方向信息
    significant_df["regulation"] = np.where(significant_df["logFC"] > 0, "Up", "Down")
    significant_df.to_csv(f"{out}_TEsubfamily.tsv",sep="\t", header=True,index=False)
    logging.info(f"meaningful TE subfamily  has been saved to {out}_TEsubfamily.tsv")
    result = (
        significant_df.groupby(["family", "regulation"])
        .size()
        .reset_index(name="count")
    )
    
    logging.info(f"save statistic to {out}_TEfamily.tsv")
    result_pivot = result.pivot(index="family", columns="regulation", values="count").fillna(0).astype(int)
    result_pivot.to_csv(f"{out}_TEfamily.tsv", sep="\t", header=True)
    logging.info(f"meaningful TE family  satistic has been saved to {out}_TEfamily.tsv")
    sns.set_style("white")  # 去掉网格线
    fig, ax = plt.subplots(figsize=(12, 6))

    # 分开取调色板：深色给Up，浅色给Down
    up_families = result[result["regulation"] == "Up"]["family"].unique()
    down_families = result[result["regulation"] == "Down"]["family"].unique()
    up_colors = sns.color_palette("dark", len(up_families))
    down_colors = sns.color_palette("pastel", len(down_families))

    # 上调在左侧
    for i, (fam, cnt) in enumerate(result[result["regulation"] == "Up"][["family", "count"]].values):
        ax.bar(i, cnt, color=up_colors[i])
        ax.text(i, cnt + 0.1, str(cnt), ha="center", va="bottom", fontsize=9)

    # 下调在右侧
    offset = len(up_families) + 2  # 空2列作为分隔
    for i, (fam, cnt) in enumerate(result[result["regulation"] == "Down"][["family", "count"]].values):
        ax.bar(offset + i, cnt, color=down_colors[i])
        ax.text(offset + i, cnt + 0.1, str(cnt), ha="center", va="bottom", fontsize=9)

    # 设置 x 轴刻度
    xticks = list(range(len(up_families))) + list(range(offset, offset + len(down_families)))
    xlabels = list(up_families) + list(down_families)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels, rotation=45, ha="right")

    ax.set_ylabel("Count", fontsize=12)
    ax.set_title("Top TE Families: Upregulated vs Downregulated", fontsize=16)

    # 重新计算 Y 轴上限，横线位置稍高于最大值

    y_max = max(result["count"])
    # 在上方标注 Up 和 Down
    ax.hlines(y=y_max*1.1, xmin=0, xmax=len(up_families)-1, color="black", linewidth=1.5)
    ax.text((len(up_families)-1)/2, y_max*1.13, "Upregulated", ha="center", va="bottom", fontsize=12)

    ax.hlines(y=y_max*1.1, xmin=offset, xmax=offset+len(down_families)-1, color="black", linewidth=1.5)
    ax.text(offset + (len(down_families)-1)/2, y_max*1.13, "Downregulated", ha="center", va="bottom", fontsize=12)

    y_max = max(result["count"]) * 1.25
    ax.set_ylim(0, y_max)
    plt.tight_layout()
    plt.savefig(f"{out}_TEfamily.png", dpi=300)
    plt.close()
    logging.info(f"save plot to {out}_TEfamily.png")

def filterTE(    
    DEG_file: str,
    rmsk_file: str
) -> pd.DataFrame:
    df_DEG = pd.read_csv(DEG_file, sep="\t")
    df_DEG.reset_index(names="subfamily", inplace=True)
    df_rmsk = pd.read_csv(rmsk_file, sep="\t", header=None)

    df_annotation = df_rmsk[[10, 11]].drop_duplicates(keep="first")
    df_annotation.columns = ["subfamily", "family"]
    
    new_df = pd.merge(df_DEG, df_annotation, on="subfamily", how="inner")
    return new_df

class GroupDict(TypedDict):
    CKO_vs_WT: str
    TRA_vs_CKO: Optional[str]
    E2_vs_CKO: Optional[str]

def heatmap_DEG(
    group_dict: GroupDict,
    adata: ad.AnnData,
    cell: str,
    out_image: str,
    rmsk_file: str = None
):
    """
    Generate a heatmap of significant genes or transposable elements (TEs) 
    across samples for a specific cell type.

    This function identifies overlapping significant features (genes/TEs) from 
    multiple differential expression result files, filters the AnnData object 
    for a given cell type, and visualizes the expression of the selected features 
    in a heatmap.

    Workflow:
    1. For each DEG result file in `DEG_list`:
       - If `rmsk_file` is None, process both genes and TEs.
       - If `rmsk_file` is provided, filter and keep only TE entries.
       - Select significant features with thresholds: |logFC| > 0.58 
         (fold change > 1.5) and -log10(PValue) > 1.3 (p < 0.05).
       - Retain only subfamily, logFC, and PValue columns.
    2. Identify overlapping significant features across all DEG files.
    3. Subset the AnnData object (`adata`) to include only the given cell type.
    4. Plot a heatmap of the overlapping features grouped by sample, 
       scaled per feature (standard_scale="var").
    5. Save the heatmap to `out_image`.

    Parameters
    ----------
    DEG_list : list of str
        List of file paths to differential expression (DEG) result tables.
    adata : ad.AnnData
        AnnData object containing expression data.
    cell : str
        Target cell type to subset from `adata`.
    out_image : str
        Path to save the output heatmap image (.png).
    rmsk_file : str, optional
        Path to RepeatMasker annotation file; if provided, only TE features are used.

    Outputs
    -------
    - Heatmap PNG file showing the expression of overlapping significant 
      genes/TEs across samples in the specified cell type.
    """
    logging.info("overlap genes of different DEG result file")
    df_list = []
    for group,infile in group_dict.items():
        if rmsk_file is None:
            logging.info("rmsk file is not provided, will plot gene and TE")
            df = pd.read_csv(infile,sep="\t")
            df.reset_index(names = "subfamily",inplace=True)
        else:
            logging.info("rmsk file is provided, will plot TE")
            df = filterTE(infile,rmsk_file)
        df = df[
            (abs(df["logFC"]) > 0.58) & 
            (-np.log10(df["PValue"]) > 1.3)
        ].copy()
        logging.info("filter meaningful feature |foldchange| > 1.5, p < 0.05")
        df = df[["subfamily","logFC","PValue"]]
        logging.info("rename the df's colname for filter")
        df.rename(columns={col: f"{group}_{col}" for col in df.columns if col != "subfamily"},inplace=True)
        df_list.append(df)
    df = reduce(lambda left, right: pd.merge(left,right,on="subfamily",how="inner"),df_list)
    logging.info("filter df, mode:\n1. CKO expression is higher than WT\n2. TRA expression is lower than CKO\n3. E2 expression is lower than CKO")
    
    keys_to_check = ["CKO_vs_WT","TRA_vs_CKO","E2_vs_CKO"]
    if all(k in group_dict.keys() for k in keys_to_check):
        df = df[(df["CKO_vs_WT_logFC"] > 0) & (df["TRA_vs_CKO_logFC"] < 0) & (df["E2_vs_CKO_logFC"] < 0)]
    elif "TRA_vs_CKO" in group_dict.keys():
        logging.info("E2_vs_CKO_logFC is not in you group_dict")
        df = df[(df["CKO_vs_WT_logFC"] > 0) & (df["TRA_vs_CKO_logFC"] < 0)]
    else:
        logging.info("TRA_vs_CKO_logFC is not in you group_dict")
        df = df[(df["CKO_vs_WT_logFC"] > 0) & (df["E2_vs_CKO_logFC"] < 0)]
    
    logging.info(f"filter adata, get the {cell} anndata obejct")
    adata = adata[adata.obs["celltype"] == cell] 
    genes_to_plot = df["subfamily"].to_list()
    logging.info(f"Find {len(genes_to_plot)} genes that fit the pattern.")
    groups = adata.obs["sample"].unique()
    plt.figure()
    sc.pl.heatmap(
        adata, 
        var_names=genes_to_plot, 
        groupby='sample',
        cmap='viridis',
        standard_scale="var",
        swap_axes=True,
        show_gene_labels=True,
        figsize=(len(groups) * 2, len(genes_to_plot) * 1.2)
    )
    fig = plt.gcf()
    fig.suptitle(cell,fontsize=16, y=0.90, x=0.5)
    # fig = plt.gcf()
    # plt.subplots_adjust(bottom=0.28, top=0.95, left=0.15, right=0.95)
    plt.savefig(out_image,dpi = 300,bbox_inches="tight")
    plt.close(fig)
    logging.info(f"heatmap plot saved to {out_image}")
    return fig


def heatmap_sc(
        adata:ad.AnnData,
        cell:str,
        out_image:str,
        rmsk_file:str
):
    df_rmsk = pd.read_csv(rmsk_file, sep="\t", header=None)
    genes_to_keep = df_rmsk[10].drop_duplicates(keep="first")

    logging.info(f"filter {cell} and TE in rmsk file")
    adata_subset = adata[
        (adata.obs["celltype"] == cell), 
        adata.var_names.isin(genes_to_keep)
    ].copy()

    if adata_subset.shape[1] == 0:
        logging.error("Warning: No genes found after filtering with df_annotation.")
        raise ValueError("Warning: No genes found after filtering with df_annotation.")
    else:
        logging.info(f"Filtered adata subset contains {adata_subset.shape[1]} genes.")

    if hasattr(adata.X, "toarray"):
        logging.info(f"adata.X is sparse matrix,convert it to a dense matrix for calculating the average expression")
        adata_subset.X = adata_subset.X.toarray()

    logging.info("Use the pandas groupby method to calculate the average expression of each group")
    
    pseudo_count = 1
    threshold = 0.58
    sample_to_check = ['WT-chang-10XSC3', 'CKO-chang-10XSC3', 'TRA-chang-10XSC3', 'E2-chang-10XSC3']

    unique_samples_set = set(adata_subset.obs["sample"].unique())
    flag = [sample in unique_samples_set for sample in sample_to_check]
    
    present_samples = [sample for sample, is_present in zip(sample_to_check, flag) if is_present]
    logging.info(f"Sample presence samples: {present_samples}")

    if not ('CKO-chang-10XSC3' in present_samples and 'WT-chang-10XSC3' in present_samples):
        logging.error("Required baseline samples 'CKO-chang-10XSC3' and 'WT-chang-10XSC3' are missing.")
        raise ValueError("Missing required samples.")
    group_averages = pd.DataFrame({
        group: np.array(adata_subset[adata_subset.obs['sample'] == group].X.mean(axis=0)).flatten()
        for group in present_samples
    }, index=adata_subset.var_names)

    logging.info("filter mode: np.log2(A + pseudo_count) - np.log2(B + pseudo_count)")
    filter_condition = np.log2(group_averages['CKO-chang-10XSC3'] + pseudo_count) - np.log2(group_averages['WT-chang-10XSC3'] + pseudo_count) > threshold
    for sample in present_samples:
        if sample != 'CKO-chang-10XSC3' and sample != 'WT-chang-10XSC3':
            logFC = np.log2(group_averages['CKO-chang-10XSC3'] + pseudo_count) - np.log2(group_averages[sample] + pseudo_count)
            filter_condition &= (logFC > threshold)
    filtered_genes = group_averages[filter_condition]
    
    
    # 1. CKO expression is higher than WT
    # 2. TRA expression is lower than CKO
    # 3. E2 expression is lower than CKO
    genes_to_plot = filtered_genes.index.tolist()
    if len(genes_to_plot) == 0:
        logging.error("the targeted gene is null, end the procedure")
        raise ValueError("the targeted gene is null, end the procedure")
    logging.info(f"Find {len(genes_to_plot)} genes that fit the pattern.")
    sc.pl.heatmap(
        adata_subset, 
        var_names=genes_to_plot, 
        groupby='sample',
        cmap='viridis',
        standard_scale="var",
        swap_axes=True,
        show_gene_labels=True
    )
    fig = plt.gcf()
    # fig.suptitle(cell,y=0.95, x=0.5)
    plt.savefig(out_image,dpi=300,bbox_inches="tight")
    plt.close(fig)
    logging.info(f"heatmap has been saved to {out_image}")

def heatmap_bulk(
        adata:ad.AnnData,
        cell:str,
        tissue:str,
        out_image:str,
        rmsk_file:str
):
    df_rmsk = pd.read_csv(rmsk_file, sep="\t", header=None)
    genes_to_keep = df_rmsk[10].drop_duplicates(keep="first")

    logging.info(f"filter {cell} and TE in rmsk file")
    adata_subset = adata[
        (adata.obs["celltype"] == cell), 
        adata.var_names.isin(genes_to_keep)
    ].copy()

    if adata_subset.shape[1] == 0:
        logging.error("Warning: No genes found after filtering with df_annotation.")
        raise ValueError("Warning: No genes found after filtering with df_annotation.")
    else:
        logging.info(f"Filtered adata subset contains {adata_subset.shape[1]} genes.")

    if hasattr(adata.X, "toarray"):
        logging.info(f"adata.X is sparse matrix,convert it to a dense matrix for calculating the average expression")
        adata_subset.X = adata_subset.X.toarray()

    logging.info("Use the pandas groupby method to calculate the average expression of each group")
    
    pseudo_count = 1
    threshold = 0.58
    sample_to_check = [f'WT_{tissue}_10XSC3', f'CKO_{tissue}_10XSC3', f'TRA_{tissue}_10XSC3', f'E2_{tissue}_10XSC3']

    unique_samples_set = set(adata_subset.obs["sample"].unique())
    flag = [sample in unique_samples_set for sample in sample_to_check]
    
    present_samples = [sample for sample, is_present in zip(sample_to_check, flag) if is_present]
    logging.info(f"Sample presence samples: {present_samples}")

    if not (f'CKO_{tissue}_10XSC3' in present_samples and f'WT_{tissue}_10XSC3' in present_samples):
        logging.error(f"Required baseline samples 'CKO_{tissue}_10XSC3' and 'WT_{tissue}_10XSC3' are missing.")
        raise ValueError("Missing required samples.")
    group_averages = pd.DataFrame({
        group: np.array(adata_subset[adata_subset.obs['sample'] == group].X.mean(axis=0)).flatten()
        for group in present_samples
    }, index=adata_subset.var_names)

    logging.info("filter mode: np.log2(A + pseudo_count) - np.log2(B + pseudo_count)")
    filter_condition = np.log2(group_averages[f'CKO_{tissue}_10XSC3'] + pseudo_count) - np.log2(group_averages[f'WT_{tissue}_10XSC3'] + pseudo_count) > threshold
    for sample in present_samples:
        if sample != f'CKO_{tissue}_10XSC3' and sample != f'WT_{tissue}_10XSC3':
            logFC = np.log2(group_averages[f'CKO_{tissue}_10XSC3'] + pseudo_count) - np.log2(group_averages[sample] + pseudo_count)
            filter_condition &= (logFC > threshold)
    filtered_genes = group_averages[filter_condition]
    
    # 1. CKO expression is higher than WT
    # 2. TRA expression is lower than CKO
    # 3. E2 expression is lower than CKO
    genes_to_plot = filtered_genes.index.tolist()
    if len(genes_to_plot) == 0:
        logging.error("the targeted gene is null, end the procedure")
        raise ValueError("the targeted gene is null, end the procedure")
    logging.info(f"Find {len(genes_to_plot)} genes that fit the pattern.")
    sc.pl.heatmap(
        adata_subset, 
        var_names=genes_to_plot, 
        groupby='sample',
        cmap='viridis',
        standard_scale="var",
        swap_axes=True,
        show_gene_labels=True
    )
    fig = plt.gcf()
    # fig.suptitle(cell,y=0.95, x=0.5)
    plt.savefig(out_image,dpi=300,bbox_inches="tight")
    plt.close(fig)
    logging.info(f"heatmap has been saved to {out_image}")

if __name__ == '__main__':
    DEG_file = "/disk5/luosg/scRNAseq/output/result/DEG/Intestine/table/combined_groupCKO_chang_10XSC3_B_cell-combined_groupWT_chang_10XSC3_B_cell_DEG.tsv"
    rmsk_file = "/ChIP_seq_2/Data/index/Mus_musculus/UCSC/mm39/rmsk_mm39.txt.gz"
    TEfamily(DEG_file,rmsk_file,"a")
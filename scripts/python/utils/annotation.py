import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import datetime
# import omicverse as ov
import celltypist
import argparse
from celltypist import models
import anndata as ad
import seaborn as sns
import logging
logging.basicConfig(
	level=logging.INFO,
	format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
	datefmt='%Y-%m-%d %H:%M:%S'
)
def Show_Markers(
        adata: ad.AnnData,
        fig_dir: str,
        table_dir: str,
        out: str,
        fig_flag: bool = False,
) -> ad.AnnData:
    """
    Perform marker gene identification and visualization for single-cell clusters.

    功能:
    1. 检查并设置 `adata.raw`（保证差异分析基于 raw 数据）。
    2. 对各 cluster (leiden) 进行差异基因分析 (Wilcoxon)。
    3. 保存所有 marker 基因表格，以及各 cluster top5 marker。
    4. 根据需要绘制相关图表 (聚类相关性、marker violin、matrixplot、dotplot、heatmap等)。
    5. 对数据进行标准化(scale)，并在 `adata.layers['scaled']` 中保存，同时基于标准化数据再次运行差异分析。
    
    参数:
        adata (ad.AnnData): AnnData 对象，含有 Leiden 聚类结果。
        fig_dir (str): 输出图像保存路径。
        table_dir (str): 输出表格保存路径。
        out (str): 文件名前缀。
        fig_flag (bool, optional): 是否绘制图像，默认 False。

    返回:
        ad.AnnData: 更新后的 AnnData 对象。
    """

    # --- Step 1: 确保 raw 存在 ---
    if adata.raw is None:
        logging.warning("Warning: adata.raw is None, will set raw = adata.copy()")
        adata.raw = adata.copy()
    
    # --- Step 2: 可选绘制 Leiden 相关性图 ---
    if fig_flag:
        sc.pl.correlation_matrix(adata, 'leiden', figsize=(5, 3.5), show=False)
        plt.savefig(f"{fig_dir}/{out}-correlation.png", dpi=300)
        plt.close()

    # --- Step 3: 差异分析 (基于 raw) ---
    sc.tl.rank_genes_groups(
        adata,
        groupby='leiden',
        method='wilcoxon',
        use_raw=True
    )
    markers = sc.get.rank_genes_groups_df(
        adata, group=None, pval_cutoff=0.05, log2fc_min=0.25
    )
    markers.to_csv(f"{table_dir}/{out}-all_markers.csv", index=False)

    # 保存 top5 marker
    top5 = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(5)
    top5.to_csv(f"{table_dir}/{out}-top5_markers.csv", index=False)

    # --- Step 4: 可选绘制 Top5 violin 图 ---
    if fig_flag:
        for col in top5.columns:  # 遍历 cluster
            fig = plt.figure(figsize=(6, 8), dpi=300)
            sc.pl.rank_genes_groups_violin(
                adata,
                groups=str(col),  # 当前 cluster
                n_genes=5,
                show=False
            )
            plt.tight_layout()
            plt.axis("off")
            fig.savefig(f"{fig_dir}/top5-markers-{col}.png")
            plt.close(fig)

    # --- Step 5: 数据标准化并存储 ---
    adata.layers['scaled'] = sc.pp.scale(adata, copy=True).X

    # 再次运行差异分析 (基于 scaled 数据)
    sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon', use_raw=False)

    # --- Step 6: 可选绘制多种 marker 图 ---
    if fig_flag:
        plots = {
            "matrix": sc.pl.rank_genes_groups_matrixplot,
            "violin": sc.pl.rank_genes_groups_stacked_violin,
            "dot": sc.pl.rank_genes_groups_dotplot,
            "heatmap": sc.pl.rank_genes_groups_heatmap
        }

        # 各图参数配置
        plot_params = {
            "matrix": dict(n_genes=3, use_raw=False, vmin=-3, vmax=3, cmap='bwr', layer='scaled'),
            "violin": dict(n_genes=3, cmap='bwr'),
            "dot": dict(n_genes=3, values_to_plot='logfoldchanges', min_logfoldchange=3,
                        vmax=7, vmin=-7, cmap='bwr'),
            "heatmap": dict(n_genes=10, use_raw=False, swap_axes=True, 
                            show_gene_labels=False, vmin=-3, vmax=3, cmap='bwr')
        }

        for name, func in plots.items():
            func(adata, show=False, **plot_params[name])
            plt.savefig(f"{fig_dir}/{out}-genes-{name}.png", dpi=300)
            plt.close()

    return adata

# def auto_omi(
#         adata: ad.AnnData,
#         out: str,
#         fig_dir: str,
#         tissue: str = "All"
# ) -> ad.AnnData:
#     print("begin auto annotate")
#     if adata.raw is None:
#         print("Warning: adata.raw is None")
#         try:
#             scsa=ov.single.pySCSA(adata = adata,
#                                 foldchange = 1.5,
#                                 pvalue = 0.01,
#                                 celltype = 'normal',
#                                 target = 'cellmarker',
#                                 tissue = tissue,
#                                 model_path ='/disk5/luosg/scRNAseq/data/pySCSA_2024_v1_plus.db'                    
#                 )
#         except Exception as e:   # Exception 是所有内置错误的基类
#             print(f"{tissue} is not in pySCSA_2024_v1_plus.db")
#             scsa=ov.single.pySCSA(adata = adata,
#                             foldchange = 1.5,
#                             pvalue = 0.01,
#                             celltype = 'normal',
#                             target = 'cellmarker',
#                             tissue = 'All',
#                             model_path ='/disk5/luosg/scRNAseq/data/pySCSA_2024_v1_plus.db'                    
#             )
#     else:
#         adata_full = adata.raw.to_adata()
#         try:
#             scsa=ov.single.pySCSA(adata = adata_full,
#                                 foldchange = 1.5,
#                                 pvalue = 0.01,
#                                 celltype = 'normal',
#                                 target = 'cellmarker',
#                                 tissue = tissue,
#                                 model_path ='/disk5/luosg/scRNAseq/data/pySCSA_2024_v1_plus.db'                    
#                 )
#         except Exception as e:   # Exception 是所有内置错误的基类
#             print(f"{tissue} is not in pySCSA_2024_v1_plus.db")
#             scsa=ov.single.pySCSA(adata = adata_full,
#                             foldchange = 1.5,
#                             pvalue = 0.01,
#                             celltype = 'normal',
#                             target = 'cellmarker',
#                             tissue = 'All',
#                             model_path ='/disk5/luosg/scRNAseq/data/pySCSA_2024_v1_plus.db'                    
#             )    

    
#     scsa.cell_anno(clustertype='leiden',
#                 cluster='all',rank_rep=True)
#     # scsa.cell_anno_print()
#     scsa.cell_auto_anno(adata,clustertype='leiden',
#                         key='scsa_celltype_cellmarker')
#     ov.utils.embedding(adata,
#                     basis='X_umap',
#                     color=[ "leiden","scsa_celltype_cellmarker"],
#                     title=['Cell type'],
#                     palette=ov.palette()[1:],
#                     show=False,frameon='small',wspace=0.35)
#     plt.tight_layout()
#     plt.legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize=8)
#     plt.savefig(f"{fig_dir}/{out}_autoAnnotation.png",dpi = 300)
#     plt.close()
#     print("end auto annotate")
#     return adata

def auto_cet(
        adata:ad.AnnData,
        fig_dir:str,
        out:str
) -> ad.AnnData:
    adata_celltypist = adata.copy()  # make a copy of our adata
    adata_celltypist.X = adata.layers["count"]  # set adata.X to raw counts
    sc.pp.normalize_per_cell(
        adata_celltypist, counts_per_cell_after=10**4
    )  # normalize to 10,000 counts per cell
    sc.pp.log1p(adata_celltypist)  # log-transform
    # make .X dense instead of sparse, for compatibility with celltypist:
    adata_celltypist.X = adata_celltypist.X.toarray()
    models.download_models(
    force_update=True, model=["Immune_All_Low.pkl", "Immune_All_High.pkl"]
    )
    model_low = models.Model.load(model="Immune_All_Low.pkl")
    model_high = models.Model.load(model="Immune_All_High.pkl")
    predictions_high = celltypist.annotate(
    adata_celltypist, model=model_high, majority_voting=True
    )
    predictions_high_adata = predictions_high.to_adata()
    adata.obs["celltypist_cell_label_coarse"] = predictions_high_adata.obs.loc[
    adata.obs.index, "majority_voting"
    ]
    adata.obs["celltypist_conf_score_coarse"] = predictions_high_adata.obs.loc[
        adata.obs.index, "conf_score"
    ]
    predictions_low = celltypist.annotate(
    adata_celltypist, model=model_low, majority_voting=True
    )
    predictions_low_adata = predictions_low.to_adata()
    adata.obs["celltypist_cell_label_fine"] = predictions_low_adata.obs.loc[
    adata.obs.index, "majority_voting"
    ]
    adata.obs["celltypist_conf_score_fine"] = predictions_low_adata.obs.loc[
        adata.obs.index, "conf_score"
    ]
    sc.pl.umap(
    adata,
    color=["celltypist_cell_label_coarse", "celltypist_conf_score_coarse"],
    frameon=False,
    sort_order=False,
    wspace=1,
    )
    plt.savefig(f"{fig_dir}/{out}-1.jpg")
    sc.pl.umap(
    adata,
    color=["celltypist_cell_label_fine", "celltypist_conf_score_fine"],
    frameon=False,
    sort_order=False,
    wspace=1,
    )
    plt.savefig(f"{fig_dir}/{out}-2.jpg")
    sc.pl.dendrogram(adata, groupby="celltypist_cell_label_fine")
    plt.savefig(f"{fig_dir}/{out}-3.jpg")

def handful_annotate(
        adata: ad.AnnData,
        marker_genes: dict[str, list[str]],
        fig_dir: str,
        out:str
) -> ad.AnnData:
    """
    Annotate single-cell RNA-seq clusters using a predefined set of marker genes,
    visualize marker expression, and assign major cell type identities.

    Workflow:
    1. Check which input marker genes are present in the dataset and keep only those.
    2. Generate a dendrogram and dotplot showing expression patterns of marker genes 
       across Leiden clusters.
    3. For each cell type, plot UMAPs highlighting expression of its marker genes 
       (normalized per gene with vmax capped at 99th percentile).
    4. Map Leiden clusters to predefined major cell type annotations and store in 
       `adata.obs['major_celltype']`.
    5. Generate UMAP plots colored by clusters and annotated cell types.
    6. Save figures and annotated AnnData object (.h5ad).

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix of shape n_obs × n_vars.
    marker_genes : dict[str, list[str]]
        Dictionary mapping cell type names to a list of marker genes.
    fig_dir : str
        Directory to save generated figures and output files.
    out : str
        Prefix for naming saved figures and output files.

    Returns
    -------
    AnnData
        The input AnnData object with an additional column 
        `obs['major_celltype']` containing the assigned annotations.
    """

    # 检查输入数据
    logging.info("Variable (genes):")
    logging.info(adata.var.head())
    logging.info("Observations (cells):")
    logging.info(adata.obs.head())

    
    if adata.raw is None:
        print("Warning: adata.raw is None")
        # marker_genes_in_data = {
        #     ct: [g for g in markers if g in adata.var.index]
        #     for ct, markers in marker_genes.items()
        # }
        # ignore upper or lower
        adata_genes_map = {g.upper(): g for g in adata.raw.var.index}  # 或 adata.var.index
        marker_genes_in_data = {}
        for ct, markers in marker_genes.items():
            filtered = []
            for g in markers:
                g_upper = g.upper()
                if g_upper in adata_genes_map:
                    filtered.append(adata_genes_map[g_upper])  # 用原名保留给绘图
            marker_genes_in_data[ct] = filtered
    else: 
        # marker_genes_in_data = {
        #     ct: [g for g in markers if g in adata.raw.var.index]
        #     for ct, markers in marker_genes.items()
        # }
        adata_genes_map = {g.upper(): g for g in adata.raw.var.index}  # 或 adata.var.index
        marker_genes_in_data = {}
        for ct, markers in marker_genes.items():
            filtered = []
            for g in markers:
                g_upper = g.upper()
                if g_upper in adata_genes_map:
                    filtered.append(adata_genes_map[g_upper])  # 用原名保留给绘图
            marker_genes_in_data[ct] = filtered
    logging.info("Filtered marker genes in data:", marker_genes_in_data)   
    

    # 计算层次聚类树并绘制 dotplot
    sc.tl.dendrogram(adata, groupby="leiden")
    sc.pl.dotplot(
        adata,
        groupby="leiden",
        var_names=marker_genes_in_data,
        dendrogram=True,
        standard_scale="var",  # 每个基因标准化到 [0, 1]
        show=False,
        use_raw=True
    )
    plt.tight_layout()
    plt.savefig(f"{fig_dir}/{out}-dot.png", dpi=300)

    # 绘制每个细胞类型的 marker 基因在 UMAP 上的分布
    for ct, genes in marker_genes_in_data.items():
        logging.info(f"Cell type: {ct}, genes: {genes}")
        if genes:
            sc.pl.umap(
                adata,
                color=["leiden"] + genes,
                vmin=0,
                legend_loc="on data",
                vmax="p99",  # 上限设为 99 分位数，避免极端值影响
                sort_order=False,  # 不强制高表达覆盖低表达
                frameon=False,
                cmap="Reds",
                show=False,
                use_raw=True
            )
            plt.tight_layout()
            plt.savefig(f"{fig_dir}/{out}-{ct}.png", dpi=300)


    return adata

def show_annotation(
    adata: ad.AnnData,
    cluster2annotation:dict[str,str],
    fig_dir:str,
    out:str
):
    print("begin map cluster")
    adata.obs['celltype'] = adata.obs['leiden'].map(cluster2annotation).astype('category')
    # prepare for color map
    leiden_to_celltype = adata.obs.groupby('leiden')['celltype'].first().to_dict()
    celltypes = adata.obs['celltype'].unique().tolist()
    palette_celltype = sns.color_palette("tab20", len(celltypes))
    celltype_to_color = {ct: c for ct, c in zip(celltypes, palette_celltype)}
    leiden_to_color = {l: celltype_to_color[ct] for l, ct in leiden_to_celltype.items()}

    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(15, 5), constrained_layout=True)
    sc.pl.umap(adata, color='leiden',palette=leiden_to_color, frameon=False, legend_loc="on data",size=8, ax=ax1, show=False)
    sc.pl.umap(adata, color='celltype',palette=celltype_to_color, frameon=False,size=10, ax=ax2, show=False)
    fig.savefig(f"{fig_dir}/{out}-handfulann.png",dpi=300,bbox_inches="tight")    
    plt.close(fig)
    logging.info("end map cluster")
    return adata


def dot_marker(
        adata:ad.AnnData,
        marker_genes:dict[str,list]
):
    adata_genes_map = {g.upper(): g for g in adata.raw.var.index}  # 或 adata.var.index
    marker_genes_in_data = {}
    for ct, markers in marker_genes.items():
        filtered = []
        for g in markers:
            g_upper = g.upper()
            if g_upper in adata_genes_map:
                filtered.append(adata_genes_map[g_upper])  # 用原名保留给绘图
        marker_genes_in_data[ct] = filtered
    logging.info("Filtered marker genes in data:", marker_genes_in_data) 
    # 画 dotplot / heatmap 看表达谱
    sc.pl.dotplot(adata, var_names=marker_genes_in_data, groupby='leiden', standard_scale='var', show=False,use_raw=True)
    plt.savefig("Muscle.png",dpi=300)
    return adata


if __name__ == '__main__':
    start = datetime.datetime.now()
    # parser = argparse.ArgumentParser(description="A script to process cell cycle data.")
    # parser.add_argument('--options', type=str, required=True, help='options to execute procedure')
    # parser.add_argument('--input', type=str, required=True, help='Path to input file')
    # parser.add_argument('--out', type=str, required=True, help='Path to out file')
    h5ad = ""
    SC = sc.read_h5ad(f"{h5ad}/SC-bbknn.h5ad")
    # TE = sc.read_h5ad(f"{h5ad}/TE-bbknn.h5ad")
    # GBM = sc.read_h5ad(f"{h5ad}/GBM-0.3-bbknn.h5ad")
    # ad_dict = {"SC":SC,"TE":TE,"GBM":GBM}
    ad_dict = {"SC":SC}
    for out,adata in ad_dict.items():
        # auto_omi(adata,f"{out}-auto")
        handful_annotate(adata,out)
        # Show_Markers(adata,out)
        # auto_cet(adata,out)
    end = datetime.datetime.now()
    logging.info("程序运行时间："+str((end-start).seconds/3600)+"h")
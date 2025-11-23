import anndata as ad
import scanpy as sc
import numpy as np
import logging
import sys
import pandas as pd
import decoupler
import matplotlib.pyplot as plt
import seaborn as sns
import textwrap

logging.basicConfig(
	level=logging.INFO,
	format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
	stream=sys.stdout,  # 指定输出到 stdout 而不是 stderr
	datefmt='%Y-%m-%d %H:%M:%S'
)
def GSEA_rankGene(
        adata:ad.AnnData,
        control:str,
        gsea_outfile:str,
        reactome_file:str = "/disk5/luosg/scRNAseq/data/reactiome.csv",
        reference:str = "rest",
        groupby:str="celltype",
        method:str = "t-test",
        use_raw:bool = False,
        decoupler_times = 10000000

):
    logging.info(f"parements: {control} vs {reference}, groupby: {groupby},method: {method}, use_raw: {use_raw}")
    sc.tl.rank_genes_groups(adata, 
                            groupby = groupby, 
                            groups = control, 
                            method=method,
                            reference = reference,
                            key_added=method,
                            use_raw=use_raw)
    celltype_condition = control
    # extract scores
    t_stats = (
        # Get dataframe of DE results for condition vs. rest
        sc.get.rank_genes_groups_df(adata, celltype_condition, key="t-test")
        # Subset to highly variable genes
        .set_index("names")
        # .loc[adata.var["highly_variable"]]
        # Sort by absolute score
        .sort_values("scores", key=np.abs, ascending=False)[
            # Format for decoupler
            ["scores"]
        ]
        .rename_axis([celltype_condition], axis=1)
    )
    # sample:stim cell population:FCGR3A+ Monocytes vs rest all cell
    logging.info("Convert gene names from lowercase to uppercase for gsea")
    t_stats.index = t_stats.index.str.upper()

    logging.info(f"the deg shape:{t_stats.shape}")
    reactome = pd.read_csv(reactome_file,sep=",")
    geneset_size = reactome.groupby("geneset").size()
    logging.info(f"filter pathway from {reactome_file}, standard: gene number is in (15,500)")
    gsea_genesets = geneset_size.index[(geneset_size > 15) & (geneset_size < 500)]
    logging.info(f"rename reactome dataframe: genesymbol:target,geneset:source for decoupler")
    reactome.rename(columns={"genesymbol":"target","geneset":"source"},inplace=True)

    scores, pvals = decoupler.mt.gsea(t_stats.T,reactome[reactome["source"].isin(gsea_genesets)],times=decoupler_times)
    gsea_result = pd.concat({"score": scores.T, "pval": pvals.T}, axis=1).droplevel(level=1,axis=1).sort_values("pval")
    gsea_result.reset_index(names = "source",inplace=True)
    gsea_result.to_csv(gsea_outfile,sep="\t",index=False)
    return gsea_result


def GSEAPlot(
        gsea_result: pd.DataFrame,
        gsea_obvious_outfile: str,
        outplot: str,
        pval_cutoff: float = 0.05
):
    # 过滤
    df = gsea_result[gsea_result["pval"] < pval_cutoff].copy()
    df.to_csv(gsea_obvious_outfile, sep="\t", index=False)

    # 处理 p-value
    df["pval_plot"] = -np.log10(df["pval"].replace(0, 1e-30))

    # 排序
    df = df.rename(columns={"score": "ES"}).sort_values("ES", ascending=True)

    # ✅ 自动 wrap 超长 pathway 名
    def wrap_label(x, width=40):
        return "\n".join(textwrap.wrap(x, width=width))

    df["source_wrapped"] = df["source"].apply(lambda x: wrap_label(x, width=35))

    # ✅ 根据 pathway 数量和名称长度动态调整画布大小
    n_pathways = df.shape[0]
    max_label_len = df["source"].str.len().max()

    height = max(4, n_pathways * 0.5)             # 行数影响高度
    width = 10 + max(0, (max_label_len - 35) / 3) # 过长名称自动加宽

    plt.figure(figsize=(width, height))

    # 绘图
    sns.scatterplot(
        data=df,
        x="ES",
        y="source_wrapped",
        size="pval_plot",
        hue="ES",
        sizes=(40, 300),
        palette="coolwarm",
        legend="brief"
    )

    # 参考线
    plt.axvline(0, linestyle='--', color='gray')

    # 图例放右侧中间
    plt.legend(
        title=None,
        bbox_to_anchor=(1.02, 0.5),
        loc="center left",
        frameon=False,
        borderaxespad=0
    )

    plt.xlabel("Enrichment Score (ES)", fontsize=12)
    plt.ylabel("Pathway", fontsize=12)

    plt.title("GSEA Enrichment Analysis", fontsize=14)

    plt.tight_layout()
    plt.savefig(outplot, dpi=300, bbox_inches="tight")
    plt.close()


def runGSEA(
        adata:ad.AnnData,
        control:str,
        gsea_outfile:str,
        gsea_obviusfile:str,
        gsea_plot:str
):
    gsea_result = GSEA_rankGene(adata,control=control,gsea_outfile=gsea_outfile)
    GSEAPlot(gsea_result,gsea_obvious_outfile=gsea_obviusfile,outplot=gsea_plot)


if __name__ == "__main__":

    runGSEA()
    



import glob
import os
import pandas as pd
import decoupler
from pathlib import Path
import scanpy as sc
import sys
basePath = Path(__file__).resolve()
baseDir = basePath.parent
sys.path.append(str(baseDir / "utils"))
from workflow.scRNAseq.scripts.downstream.condition.GSEA import GSEAPlot
adata_dict = {
    "Intestine": sc.read_h5ad("/disk5/luosg/scRNAseq/output/combine/Intestine/Intestine_annotate.h5ad"),
    "Lung": sc.read_h5ad("/disk5/luosg/scRNAseq/output/combine/Lung/Lung_annotate.h5ad")
}
TE = "/disk5/luosg/Reference/UCSC/mouse/mm39/rmsk_mm39.txt.gz"
df_TE = pd.read_csv(TE,sep="\t")
TE_name = df_TE["repName"].unique()
del df_TE
TE_name
fileStr = 'combined_group*'
outdir = Path("/disk5/luosg/scRNAseq/output/combine")
reactome = pd.read_csv("/disk5/luosg/scRNAseq/data/reactiome.csv",sep=",")
geneset_size = reactome.groupby("geneset").size()
gsea_genesets = geneset_size.index[(geneset_size > 15) & (geneset_size < 500)]
reactome.rename(columns={"genesymbol":"target","geneset":"source"},inplace=True)

def prepare_file():
    for group in adata_dict.keys():
        search_path = outdir / f"{group}/DEG/table"
        fig_dir = outdir / f"{group}/DEG/volcano/Gene"
        files = glob.glob(os.path.join(str(search_path), fileStr))
        for file in files:
            file_name = file.split('.')[0].split('/')[-1]
            df = pd.read_csv(file,sep="\t")
            df.reset_index(names="gene",inplace=True)
            df = df[~df["gene"].isin(TE_name)]
            print(df.head())
            # df = df[(df["PValue"] < 0.05) & (df["FDR"] < 0.05) & (abs(df["logFC"]) > 0.58)]
            DEG_gene_dir = outdir / f"{group}/GSEA/gene"
            DEG_gene_dir.mkdir(parents=True,exist_ok=True)
            DEG_gene_outfile = outdir / f"{group}/GSEA/gene/{file_name}_gene.tsv"
            df.to_csv(DEG_gene_outfile,sep="\t",index=False)
def runGSEA():
    """
    single cell用高可变基因
    bulk得用全部差异基因
    """
    for group in adata_dict.keys():
        search_path = outdir / f"{group}/GSEA/gene"
        files = glob.glob(os.path.join(str(search_path), fileStr))
        for file in files:
            file_name = file.split('.')[0].split('/')[-1]
            df = pd.read_csv(file,sep="\t",index_col=0)
            df = df[["logFC"]]
            df.index = df.index.str.upper()
            df.columns = ["score"]
            figdir = outdir / f"{group}/GSEA/plot"
            figdir.mkdir(parents=True,exist_ok=True)
            outJpeg = outdir / f"{group}/GSEA/plot/{file_name}.png"
            gsea_outdir = outdir / f"{group}/GSEA/table"
            gsea_outdir.mkdir(parents=True,exist_ok=True)
            gsea_outfile = gsea_outdir / f"{file_name}_gsea.tsv"
            gsea_obvious_outfile = gsea_outdir / f"{file_name}_gsea_obvious.tsv"
            try:
                scores, pvals = decoupler.mt.gsea(df.T,reactome[reactome["source"].isin(gsea_genesets)])
            except Exception as e:
                print(f"{file_name}\n{e}")
                continue
            gsea_result = pd.concat({"score": scores.T, "pval": pvals.T}, axis=1).droplevel(level=1,axis=1).sort_values("pval")
            gsea_result.reset_index(names = "source",inplace=True)
            gsea_result.to_csv(gsea_outfile,sep="\t",index=False)
            GSEAPlot(gsea_result=gsea_result,gsea_obvious_outfile=gsea_obvious_outfile,outplot=outJpeg)

if __name__ == "__main__":
    # prepare_file()
    runGSEA()
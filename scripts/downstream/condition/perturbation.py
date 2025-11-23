import pertpy as pt
import scanpy as sc
import logging
import sys
import anndata as ad
import matplotlib.pyplot as plt
from pathlib import Path
logging.basicConfig(
	level=logging.INFO,
	format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
	stream=sys.stdout,  # 指定输出到 stdout 而不是 stderr
	datefmt='%Y-%m-%d %H:%M:%S'
)

def runAngur(
        adata:ad.AnnData,
        fig_dir:str,
        condition_label:str,
        treatment_label:str,
        cell_type_key:str="celltype",
        condition_key:str="condition"
):
    logging.info("Angur analysis begin")
    if adata.raw is not None:
        adata_raw = adata.raw.to_adata()
        adata_raw.obs = adata.obs.copy()
        adata_raw.obsm = adata.obsm.copy()
        adata_raw.uns = adata.uns.copy()
        adata = adata_raw
        del adata_raw
    fig_dir = Path(fig_dir)
    fig_dir.mkdir(parents=True,exist_ok=True)
    logging.info("train data")
    ag_rfc = pt.tl.Augur("random_forest_classifier")
    loaded_data = ag_rfc.load(adata, label_col=condition_key, cell_type_col=cell_type_key,condition_label = condition_label, treatment_label = treatment_label)
    
    v_adata, v_results = ag_rfc.predict(
    loaded_data, subsample_size=20, n_threads=4, select_variance_features=True, span=1)
    h_adata, h_results = ag_rfc.predict(loaded_data, subsample_size=20, select_variance_features=False, n_threads=4)
    scatter = ag_rfc.plot_scatterplot(v_results, h_results)
    plt.savefig(fig_dir / f"{treatment_label}_vs_{condition_label}-gene-all_vs_high.png",dpi=300,bbox_inches="tight")
    plt.close()
    v_results["summary_metrics"]

    logging.info("plot result")
    ag_rfc.plot_lollipop(v_results)
    plt.savefig(fig_dir / f"{treatment_label}_vs_{condition_label}-lollipop_gene-all.png",dpi=300,bbox_inches="tight")
    plt.close()
    important_features = ag_rfc.plot_important_features(v_results)
    plt.savefig(fig_dir / f"{treatment_label}_vs_{condition_label}-important_features-all.png",dpi=300,bbox_inches="tight")
    plt.close()

    logging.info("Angur analysis end")

if __name__ == "__main__":
    adata_Intestine = sc.read_h5ad("/disk5/luosg/scRNAseq/output/combine/Intestine/Intestine_annotate.h5ad")
    adata_Lung = sc.read_h5ad("/disk5/luosg/scRNAseq/output/combine/Lung/Lung_annotate.h5ad")
    runAngur(adata_Intestine,"/disk5/luosg/scRNAseq/output/combine/Intestine/pertubation",condition_label="WT",treatment_label="CKO",condition_key="experiment")
    runAngur(adata_Intestine,"/disk5/luosg/scRNAseq/output/combine/Intestine/pertubation",condition_label="CKO",treatment_label="TRA",condition_key="experiment")
    runAngur(adata_Intestine,"/disk5/luosg/scRNAseq/output/combine/Intestine/pertubation",condition_label="CKO",treatment_label="E2",condition_key="experiment")
    runAngur(adata_Lung,"/disk5/luosg/scRNAseq/output/combine/Lung/pertubation",condition_label="WT",treatment_label="CKO",condition_key="experiment")
    runAngur(adata_Lung,"/disk5/luosg/scRNAseq/output/combine/Lung/pertubation",condition_label="CKO",treatment_label="TRA",condition_key="experiment")
    runAngur(adata_Lung,"/disk5/luosg/scRNAseq/output/combine/Lung/pertubation",condition_label="CKO",treatment_label="E2",condition_key="experiment")




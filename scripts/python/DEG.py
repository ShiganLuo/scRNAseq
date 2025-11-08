from pathlib import Path
import sys
import os
from typing import Literal
os.environ["R_HOME"] = "/home/luosg/miniconda3/envs/scRNAseq_rpy2_1/lib/R"
os.environ["PATH"] = "/home/luosg/miniconda3/envs/scRNAseq_rpy2_1/bin:" + os.environ.get("PATH", "")
os.environ["LD_LIBRARY_PATH"] = "/home/luosg/miniconda3/envs/scRNAseq_rpy2_1/lib:" + os.environ.get("LD_LIBRARY_PATH", "")
basePath = Path(__file__).resolve()
baseDir = basePath.parent
sys.path.append(str(baseDir / "utils"))
import anndata as ad
import numpy as np
import matplotlib.pyplot as plt
import logging
import pandas as pd
import scanpy as sc
import glob
from utils.TE import TEfamily
from rpy2.robjects.conversion import localconverter
import rpy2.robjects.packages as rpackages
import rpy2.robjects as ro
import rpy2.rinterface_lib.callbacks
from rpy2.robjects import pandas2ri
sc.settings.verbosity = 0
rpy2.rinterface_lib.callbacks.logger.setLevel(logging.ERROR)


logging.basicConfig(
    level=logging.INFO, 
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
# logging.debug("This is a debug message.") 
# logging.info("This is an info message.")
# logging.warning("This is a warning message.")
# logging.error("This is an error message.")
basePath = Path().resolve()
baseDir = basePath.parent
sys.path.append(str(baseDir / "utils"))
from utils.DEG_pseudo_bulk import Pre_AddAnnotations2RawCountsAnndata,Pre_PseduoBulk,Pre_FeatureRename

def makePseudoBulk(
    adataRaw:ad.AnnData,
    adataAnn:ad.AnnData,
    fig_dir: str,
    out:str,
    fig_flag:bool=False,
    obs_to_keep:list = ['sample',"celltype"]
):
    """Generates a pseudo-bulk AnnData object from single-cell data.

    This function processes single-cell data by integrating annotations, creating 
    pseudo-bulk samples, and performing necessary standardization steps for 
    downstream differential gene expression analysis.

    Parameters
    ----------
    adataRaw
        An AnnData object containing raw single-cell count data.
    adataAnn
        An AnnData object with cell type annotations in its `obs` attribute.
    fig_dir
        The directory to save the generated PCA plot.
    out
        The filename (without extension) for the saved PCA plot.
    fig_flag
        A boolean flag to determine whether to generate and save a PCA plot.
        Defaults to `False`.
    obs_to_keep
        A list of observation columns to be included in the final pseudo-bulk object.
        Defaults to `['sample', 'celltype']`.

    Returns
    -------
    ad.AnnData
        The pseudo-bulk AnnData object containing aggregated counts, normalized data,
        and relevant observation metadata, ready for analysis.
    
    Notes
    -----
    - The function first combines annotations from `adataAnn` with the raw counts
      in `adataRaw`.
    - It renames characters in specified observation columns to ensure compatibility
      with tools like edgeR.
    - Pseudo-bulk aggregation is performed for each unique cell type.
    - The resulting pseudo-bulk data is normalized, log-transformed, and scaled
      for PCA.
    - Library size and log library size are added to `adata.obs`.
    - If `fig_flag` is `True`, a PCA plot is generated and saved for quality control.
    """
    logging.info(f"Begin make PseduoBulk adata!")
    ######## migrate celltyep annotation to adataRaw
    logging.info(f'begin Migrate celltype annotation from adataAnn to adataRaw and change name in {obs_to_keep} for edgeR')
    adata = Pre_AddAnnotations2RawCountsAnndata(adataRaw,adataAnn)
    # edgeR can't identify " " and "+"
    adata = Pre_FeatureRename(adata,obs_to_keep)
    logging.info(f"Migrate complete! New adata {adata}")
    logging.info(f"The max count of adata : {np.max(adata.X)}")
    ######## PseudoBulk
    logging.info(f"Make PseudoBulk adata formally")
    adata.obs["sample"] = adata.obs["sample"].astype("category")
    adata.obs["celltype"] = adata.obs["celltype"].astype("category")
    celltype = adata.obs["celltype"].cat.categories[0]
    adata_list = []
    for i, celltype in enumerate(adata.obs["celltype"].cat.categories[0:]):
        print(
            f'Processing {celltype} ({i+2} out of {len(adata.obs["celltype"].cat.categories)})...'
        )
        adata_cell_type = Pre_PseduoBulk(adata, celltype, obs_to_keep=obs_to_keep,replicates_per_patient=2)
        adata_list.append(adata_cell_type)
    adata_pb = ad.concat(adata_list,index_unique=None)
    logging.info(f"PseudoBulk adata make complete Preliminarily")
    ####### check new pseduoBulk adata
    logging.info(f"Begin check PseduoBulk adata")
    adata_pb.layers['counts'] = adata_pb.X.copy()
    print(adata_pb)
    print(np.max(adata_pb.X))
    sc.pp.normalize_total(adata_pb, target_sum=1e4)
    sc.pp.log1p(adata_pb)
    sc.pp.pca(adata_pb)
    print(np.max(adata_pb.X))
    print(np.max(adata_pb.layers['counts']))
    ## add log_lib_size
    # If counts is an ndarray or sparse matrix
    counts = adata_pb.layers["counts"]
    if hasattr(counts, "toarray"):
        counts = counts.toarray()

    lib_size = np.sum(counts, axis=1).flatten()
    adata_pb.obs["lib_size"] = lib_size

    # log_lib_size: first convert to numpy array, then log
    log_lib_size = np.log(np.array(adata_pb.obs["lib_size"], dtype=float))
    adata_pb.obs["log_lib_size"] = log_lib_size
    print(adata_pb.obs["log_lib_size"])
    if fig_flag:
        logging.info(f"Plot PCA plot of column in PseduoBulk adata's obs")
        sc.pl.pca(adata_pb, color=adata_pb.obs, ncols=1, size=300)
        plt.savefig(f"{fig_dir}/{out}_obsPCA.png")
    
    logging.info(f"After normalizing and log1p; The max count of adata:{np.max(adata_pb.X)}")
    adata_pb.X = adata_pb.layers['counts'].copy()
    logging.info(f"After copying layers 'counts'; The max count of adata:{np.max(adata_pb.X)}")
    logging.info(f"PseduoBulk adata check complete !")
    logging.info(f"PseduoBulk adata make complete!")
    return adata_pb

def runEdgeRForCell(
    adata: ad.AnnData,
    cell: str,
    fig_dir: str,
    table_dir: str,
    out: str,
    DEGDesign:str,
    cell_type_key: str = "celltype",
    r_script_path: str = "/disk5/luosg/scRNAseq/workflow/scRNAseq/scripts/python/utils/DEG_pseduo_bulk.r"
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Performs pseudo-bulk differential gene expression (DEG) analysis for a specific cell type
    using a pre-defined R script and the rpy2 library.

    This function isolates a single cell type from an AnnData object, aggregates the data
    to a pseudo-bulk level, and then calls a master R function (`fit_model`) to handle the
    complete DEG analysis workflow, including filtering, normalization, model fitting,
    statistical testing, and saving plots and results.

    Args:
        adata (ad.AnnData): The input AnnData object containing single-cell data.
        cell (str): The specific cell type to be analyzed.
        fig_dir (str): The directory path to save diagnostic plots generated by the R script.
        table_dir (str): The directory path to save the differential expression results table.
        out (str): A base name for output files (plots and tables).
        DEGDesign (str): A character string defining the contrast for the DEG comparison in R.
                         e.g., "GroupA-GroupB".
        cell_type_key (str, optional): The key in `adata.obs` that contains the cell type
                                       annotations. Defaults to "celltype".
        r_script_path (str, optional): The full path to the R script containing the `fit_model`
                                       function. Defaults to a predefined path.

    Returns:
        tuple[pd.DataFrame, pd.DataFrame]: A tuple containing the prepared expression matrix
                                           (genes x cells) and metadata DataFrame for the
                                           subsetted cell type.

    Raises:
        KeyError: If the specified `cell_type_key` is not found in `adata.obs`.
        ValueError: If the specified `cell` type is not found in the `cell_type_key` column.
        rpy2.rinterface_lib.embedded.RRuntimeError: If the R script fails to load or if
                                                    R functions encounter a runtime error.
    """
    if cell_type_key not in adata.obs.columns:
        logging.error(f"The specified key '{cell_type_key}' is not found in adata.obs.")
        raise KeyError(f"Required key '{cell_type_key}' not found.")

    if cell not in adata.obs[cell_type_key].unique():
        logging.error(f"The specified cell type '{cell}' is not present in the data for key '{cell_type_key}'.")
        raise ValueError(f"Cell type '{cell}' not found in `adata.obs['{cell_type_key}']`.")

    logging.info(f"Subsetting data for cell type '{cell}' using key '{cell_type_key}'.")
    adata_mono = adata[adata.obs[cell_type_key] == cell].copy()

    X = adata_mono.X
    if hasattr(X, "toarray"):
        logging.info("Converting sparse matrix to a dense array.")
        X = X.toarray()
    
    expr_df = pd.DataFrame(
        X.T,  # edgeR and similar tools expect genes x cells
        index=adata_mono.var_names,
        columns=adata_mono.obs_names,
        dtype=np.float32  # Use a memory-efficient numerical type
    )
    logging.info(f"Expression matrix shape: {expr_df.shape}")

    coldata_df = adata_mono.obs.copy()

    try:
        coldata_df['replicate'] = coldata_df.index.to_series().str.extract(r'_([0-9]+)$')[0].astype(int)
    except KeyError:
        logging.warning("Could not extract 'replicate' from cell index. Check cell barcode format.")
        pass

    logging.info(f"Extracted metadata DataFrame with shape: {coldata_df.shape}")
    
    ro.r(f'source("{r_script_path}")')
    logging.info("Convert DataFrame from Python to R ")
    with localconverter(ro.default_converter + pandas2ri.converter):
        r_expr_df = ro.conversion.py2rpy(expr_df) # can't transpose
        r_coldata_df = ro.conversion.py2rpy(coldata_df)
    
    logging.info("Call R custom function fit_model to conduct differential gene analysis")
    logging.info(ro.r('Sys.getenv("LD_LIBRARY_PATH")'))
    fit_model_r_func = ro.globalenv['fit_model']
    outs = fit_model_r_func(
        r_expr_df,
        r_coldata_df,
        group_col="sample",
        DEGDesign=DEGDesign,
        out_prefix=out,
        fig_dir=fig_dir,
        table_dir=table_dir,
        replicate_col="replicate"
    )
    
    logging.info("Differential gene analysis complete")
    return expr_df, coldata_df

def DEGDesignGenerate(
        adata:ad.AnnData
):
    adata.obs['sample_celltype'] = adata.obs['sample'].astype(str) + '_' + adata.obs['celltype'].astype(str)
    grouped_by_celltype = adata.obs.groupby('celltype',observed=False)

    result_dict = {}
    for celltype, group in grouped_by_celltype:
        result_dict[celltype] = group['sample_celltype'].tolist()
        
    deduplicated_dict = {}
    for key, value_list in result_dict.items():
        unique_list = list(set(value_list))
        deduplicated_dict[key] = unique_list

    # print(deduplicated_dict)

    result_dict = {}

    # 遍历输入字典的每一项
    for cell_type, cell_list in deduplicated_dict.items():
        
        # 初始化一个空列表，用于存储当前 cell_type 的所有组合
        combinations = []
        
        # 筛选出每种类型的字符串
        cko_strings = [s for s in cell_list if s.startswith('CKO')]
        wt_strings = [s for s in cell_list if s.startswith('WT')]
        tra_strings = [s for s in cell_list if s.startswith('TRA')]
        e2_strings = [s for s in cell_list if s.startswith('E2')]

        # # 生成 CKO-WT 组合
        # for cko_s in cko_strings:
        #     for wt_s in wt_strings:
        #         combinations.append(f"combined_group{cko_s}-combined_group{wt_s}")
                
        # # 生成 TRA-CKO 组合
        # for tra_s in tra_strings:
        #     for cko_s in cko_strings:
        #         combinations.append(f"combined_group{tra_s}-combined_group{cko_s}")

        # # 生成 E2-CKO 组合
        # for e2_s in e2_strings:
        #     for cko_s in cko_strings:
        #         combinations.append(f"combined_group{e2_s}-combined_group{cko_s}")

        # 生成TRA-WT组合
        for tra_s in tra_strings:
            for wt_s in wt_strings:
                combinations.append(f"combined_group{tra_s}-combined_group{wt_s}")
        
        #生成E2-WT组合
        for e2_s in e2_strings:
            for wt_s in wt_strings:
                combinations.append(f"combined_group{e2_s}-combined_group{wt_s}")
        # 将生成的组合列表作为值，以 cell_type 为键存入结果字典
        result_dict[cell_type] = combinations
                
    print(result_dict)
    return result_dict

def main(
        adata:ad.AnnData,
        fig_dir:str,
        table_dir:str
):
    logging.info(ro.r('.libPaths()'))
    logging.info(ro.r('system.file(package = "edgeR")'))
    adata_design = DEGDesignGenerate(adata)
    for cell,designs in adata_design.items():
        for design in designs:
            logging.info(f"{cell}:{design}")
            expr_df,coldata_df = runEdgeRForCell(adata,
                                        cell = cell,
                                        fig_dir=fig_dir,
                                        table_dir=table_dir,
                                        out=design,
                                        DEGDesign=design)

def runVolcano(
        infile:str,
        out:str,
        mode:Literal["TE","Gene","None"] = "TE",
        TE:str = "/disk5/luosg/Reference/UCSC/mouse/mm39/rmsk_mm39.txt.gz",
        r_script_path: str = "/disk5/luosg/scRNAseq/workflow/scRNAseq/scripts/python/utils/volcano.r"
):
    logging.info("Begin read DEG result and rmsk file")
    df = pd.read_csv(infile,sep="\t")
    df["gene_name"] = df.index
    logging.info(f"mode: {mode}")
    if mode == "TE":
        TEfamily(DEG_file=infile,rmsk_file=TE,out=out)
        df_TE = pd.read_csv(TE,sep="\t",header=None)
        keep_vec = df_TE[10]
        logging.info("Filter out none TE gene_name according to rmsk file")
        transposon_index_to_keep = df.index.intersection(keep_vec.unique()) #remove duplicated value
        df_filtered = df.loc[transposon_index_to_keep]
    elif mode == "Gene":
        df_TE = pd.read_csv(TE,sep="\t",header=None)
        except_vec = df_TE[10]
        logging.info("Filter out TE gene_name according to rmsk file")
        transposon_index_to_keep = df.index[~df.index.isin(except_vec.unique())]
        df_filtered = df.loc[transposon_index_to_keep]
    else:
        df_filtered = df
    print(df_filtered.head())
    ro.r(f'source("{r_script_path}")')
    logging.info("Convert DataFrame from Python to R ")
    with localconverter(ro.default_converter + pandas2ri.converter):
        r_df = ro.conversion.py2rpy(df_filtered)
    
    logging.info("Call R custom function volcano to plot")
    volcano_r_func = ro.globalenv['volcano']
    volcano_r_func(r_df,outjpeg = f"{out}.jpeg")
    


if __name__ == '__main__':
    indir = "/disk5/luosg/scRNAseq/output/h5ad_QC"
    outdir = "/disk5/luosg/scRNAseq/output/combine"
    samplesDict = { 'Intestine': ['CKO-chang-10XSC3', 'E2-chang-10XSC3', 'TRA-chang-10XSC3', 'WT-chang-10XSC3'],
        'Lung': ['CKO-fei-10XSC3', 'E2-fei-10XSC3', 'TRA-fei-10XSC3', 'WT-fei-10XSC3']}
    # 'Muscle': ['CKO-jirou-10XSC3', 'E2-jirou-10XSC3', 'TRA-jirou-10XSC3', 'WT-jirou-10XSC3']
    ## to bulk adata
    # for group in samplesDict.keys():
    #     adataRaw = ad.read_h5ad(f"{outdir}/{group}/{group}_qc.h5ad")
    #     adataAnn = ad.read_h5ad(f"{outdir}/{group}/{group}_annotate.h5ad")
    #     obs_to_keep = ["celltype", "sample"]
    #     fig_dir = Path(f"{outdir}/{group}/DEG")
    #     fig_dir.mkdir(parents=True, exist_ok=True)
    #     adata = makePseudoBulk(adataRaw,adataAnn,fig_dir=str(fig_dir),out=group,fig_flag=True)
    #     adata.write_h5ad(f"{outdir}/{group}/{group}_bulk.h5ad")

    #### DEG


    # adata = ad.read_h5ad("/disk5/luosg/scRNAseq/output/result/DEG/Intestine/h5ad/Intestine_bulk.h5ad")
    # expr_df,coldata_df = runEdgeRForCell(adata,
    #                                      cell = "B_cell",
    #                                      fig_dir="/disk5/luosg/scRNAseq/output/result/DEG/Intestine/plot",
    #                                      table_dir="/disk5/luosg/scRNAseq/output/result/DEG/Intestine/table",
    #                                      out="Intestine_B",
    #                                      DEGDesign="combined_groupCKO_chang_10XSC3_B_cell-combined_groupWT_chang_10XSC3_B_cell")
    
    # DEGDesignGenerate(Intestine)
    # adataDcit = {}
    # Intestine = ad.read_h5ad("/disk5/luosg/scRNAseq/output/result/DEG/Intestine/h5ad/Intestine_bulk.h5ad")
    # Lung = ad.read_h5ad("/disk5/luosg/scRNAseq/output/result/DEG/Lung/h5ad/Lung_bulk.h5ad")
    # # Muscle = ad.read_h5ad("/disk5/luosg/scRNAseq/output/result/DEG/Muscle/h5ad/Muscle_bulk.h5ad")
    # adataDcit['Intestine'] = Intestine
    # adataDcit['Lung'] = Lung
    # # adataDcit['Muscle'] = Muscle
    # outdir = Path("/disk5/luosg/scRNAseq/output/combine")
    # for tissue,adata in adataDcit.items():
    #     fig_dir = outdir / f"{tissue}/DEG/plot"
    #     fig_dir.mkdir(parents=True, exist_ok=True)
    #     table_dir = outdir /f"{tissue}/DEG/table"
    #     table_dir.mkdir(parents=True, exist_ok=True)
    #     main(adata,
    #         fig_dir=str(fig_dir),
    #         table_dir=str(table_dir))

    # runVolcano("/disk5/luosg/scRNAseq/output/result/DEG/Intestine/table/combined_groupCKO_chang_10XSC3_B_cell-combined_groupWT_chang_10XSC3_B_cell_DEG.tsv"
    #            ,"a.jepg")
    fileStr = 'combined_group*'
    outdir = Path("/disk5/luosg/scRNAseq/output/combine")
    for group in samplesDict.keys():
        search_path = outdir / f"{group}/DEG/table"
        fig_dir = outdir / f"{group}/DEG/volcano/Gene"
        fig_dir.mkdir(parents=True,exist_ok=True)
        files = glob.glob(os.path.join(str(search_path), fileStr))
        for file in files:
            file_name = file.split('.')[0].split('/')[-1]
            outJpeg = f"{fig_dir}/{file_name}"
            runVolcano(file,outJpeg,mode="Gene")
        
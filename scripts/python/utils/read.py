import numpy as np
import scanpy as sc
import anndata as ad
import pandas as pd
import scipy as sp
import argparse
import datetime
import logging
import os
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

def sparsifyRead(filename:str,sample:str):
    data = pd.read_csv(filename, index_col=0, header=0)
    data.index = data.index.astype(str)
    genes = data.columns
    cells = data.index
    data = sp.sparse.csr_matrix(data.to_numpy())
    data.astype('float32')

    '''
    oh = open('gene_names.{0}.tsv'.format(os.path.split(filename)[1]), 'w')
    for g in genes:
        oh.write('%s\n' % g)
    oh.close()
    '''

    print('Loaded {0}'.format(filename))
    andata = ad.AnnData(data, obs={'obs_names': cells}, var={'var_names': genes})
    del data
    andata.obs['n_genes'] = (andata.X > 0).sum(axis=1).A1  # 计算每个细胞的基因数量
    andata.obs['n_counts'] = andata.X.sum(axis=1).A1 #axis=1沿行方向
    andata.obs['sample'] = np.full(andata.n_obs, sample)
    # ad.var_names_make_unique()
    # ad.obs_names_make_unique()
    return andata
# the output of cellranger
def cellrangerRead(
        filePath: str,
        sample: str,
        method: str = 'h5'
) -> ad.AnnData:
    """
    Reads 10x Genomics data (Cell Ranger output) using either the sparse Matrix Market 
    format (.mtx) or the more efficient HDF5 format (.h5), and annotates the result 
    with sample information.

    The function validates the input method and ensures the filePath is the correct 
    type (file or directory) for the chosen reading method.

    Parameters
    ----------
    filePath
        `str`. The path to the Cell Ranger output.
        If ``method='mtx'``, this must be the **directory** containing matrix.mtx, barcodes.tsv, etc.
        If ``method='h5'``, this must be the **path to the single .h5 file**.
    sample
        `str`. The name of the sample to assign to all cells in the resulting AnnData object.
    method
        `str`, optional (default: 'h5'). The file format to read. Must be 'mtx' or 'h5'. 
        'h5' is generally much faster.

    Returns
    -------
    ad.AnnData
        An AnnData object containing the raw count matrix in ``.X``, where gene symbols 
        are used as variable names, and an observation column ``adata.obs['sample']`` 
        is added with the specified sample name.

    Raises
    ------
    ValueError
        If an unsupported reading method is provided (i.e., not 'mtx' or 'h5').
    FileNotFoundError
        If the ``filePath`` path does not exist, or if the path exists but is the wrong 
        type (e.g., trying to read an .h5 file as an 'mtx' directory).

    Notes
    -----
    When using ``method='mtx'``, ``cache=True`` is used for faster subsequent reading.
    """
    if method not in ('mtx', 'h5'):
        raise ValueError(f"Invalid method '{method}'. Must be 'mtx' or 'h5'.")

    if not os.path.exists(filePath):
        raise FileNotFoundError(f"The path specified for Cell Ranger output does not exist: '{filePath}'")

    if method == 'h5':
        if not os.path.isfile(filePath):
            raise FileNotFoundError(f"For method 'h5', the path must point to a file, but '{filePath}' is a directory.")
    elif method == 'mtx':
        if not os.path.isdir(filePath):
            raise FileNotFoundError(f"For method 'mtx', the path must point to a directory, but '{filePath}' is a file.")


    # Read data based on the chosen method
    if method == 'h5':
        # The sc.read_10x_h5 function is for the single file format and is usually faster
        # For h5, var_names are typically already set correctly or handled internally.
        adata = sc.read_10x_h5(filePath)
    elif method == 'mtx':
        # The sc.read_10x_mtx function is for the directory format
        adata = sc.read_10x_mtx(
            filePath,
            var_names="gene_symbols",
            cache=True)

    # Annotate with sample information
    adata.obs['sample'] = np.full(adata.n_obs, sample)
    return adata

# the output of scTE
def scTERead(filename:str,sample:str):
    logging.info(f"begin read the output of scTE:\n file path:{filename}\tsampe name:{sample}")
    data = pd.read_csv(filename, index_col=0, header=0)
    data.index = data.index.astype(str)
    genes = data.columns
    cells = data.index
    data = sp.sparse.csr_matrix(data.to_numpy())
    data.astype('float32')
    adata = ad.AnnData(data, obs={'obs_names': cells}, var={'var_names': genes})
    del data
    adata.obs['sample'] = np.full(adata.n_obs, sample)
    # ad.obs['n_TE'] = (ad.X > 0).sum(axis=1).A1  # 计算每个细胞的转座子数量
    # ad.obs['n_counts'] = ad.X.sum(axis=1).A1 #axis=1沿行方向
    # ad.var_names_make_unique()
    # ad.obs_names_make_unique()
    print("read complete")
    return adata

if __name__ == '__main__':
    start = datetime.datetime.now()
    # parser = argparse.ArgumentParser(description="A script to process cell cycle data.")
    # parser.add_argument('--options', type=str, required=True, help='options to execute procedure')
    # parser.add_argument('--input', type=str, required=True, help='Path to input file')
    # parser.add_argument('--out', type=str, required=True, help='Path to out file')
    h5ad = "/home/lsg/Data/glioblastoma/output/new/h5ad"
    samples = ["GBM27","GBM28","GBM29"]
    for i in samples:
        SC_file = f"/home/lsg/Data/glioblastoma/output/cellranger/{i}/outs/filtered_feature_bc_matrix"
        TE_file = f"/home/lsg/Data/glioblastoma/output/scTE/cellranger_yes/{i}.csv"
        print(f"read: {SC_file}")
        print(f"read: {TE_file}")
        ad_SC = cellrangerRead(SC_file,i)
        ad_TE = scTERead(TE_file,i)
        ad_SC.write_h5ad(f"{h5ad}/{i}-SC-raw.h5ad")
        ad_TE.write_h5ad(f"{h5ad}/{i}-TE-raw.h5ad")
    end = datetime.datetime.now()
    print("程序运行时间："+str((end-start).seconds/3600)+"h")
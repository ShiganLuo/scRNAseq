import anndata as ad
import pandas as pd
import random
import scanpy as sc
import numpy as np
import logging
import matplotlib.pyplot as plt
logging.basicConfig(
    level=logging.INFO, 
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
# logging.debug("This is a debug message.") 
# logging.info("This is an info message.")
# logging.warning("This is a warning message.")
# logging.error("This is an error message.")

import numpy as np
import scipy.sparse as sp
import logging

def is_raw_counts(adata, layer=None):
    """
    Determine whether adata.X or a specified layer looks like raw counts matrix.
    Works for dense and sparse matrices.

    Parameters
    ----------
    adata : AnnData
    layer : str or None

    Returns
    -------
    bool
    """
    if layer is None:
        X = adata.X
        name = "adata.X"
    else:
        if layer not in adata.layers:
            raise ValueError(f"Layer '{layer}' does not exist in adata.layers")
        X = adata.layers[layer]
        name = f"adata.layers['{layer}']"

    # ============ 稀疏矩阵处理 ============
    if sp.issparse(X):
        # 负值判断
        if X.min() < 0:   # csr/csc 稀疏矩阵支持 min()
            logging.info(f"{name} contains negative values; not raw counts.")
            return False

        # 判断是否整数：稀疏矩阵只需要判断 data 数组即可
        data = X.data
        if not np.all(np.floor(data) == data):
            logging.info(f"{name} contains decimals in sparse data; not raw counts.")
            return False

        logging.info(f"{name} is likely raw counts (sparse non-negative ints).")
        return True

    # ============ dense matrix 处理 ============
    else:
        if np.min(X) < 0:
            logging.info(f"{name} contains negative values; not raw counts.")
            return False

        if not np.all(np.floor(X) == X):
            logging.info(f"{name} contains decimals; not raw counts.")
            return False

        logging.info(f"{name} is likely raw counts (dense non-negative ints).")
        return True



def Pre_AddAnnotations2RawCountsAnndata(
    adataRaw: ad.AnnData,
    adataAnn: ad.AnnData,
    annotation_keys: list[str] = ["celltype"],
    raw_layers: str = "counts",
    min_overlap_ratio: float = 0.8  # Minimum fraction of overlapping cells
) -> ad.AnnData:
    """
    Copy multiple annotation columns from one AnnData object to another raw counts AnnData object.
    Only keeps cells that have annotations in adataAnn.
    """
    if not is_raw_counts(adataRaw):
        logging.error("adataRaw.X is not raw counts data.")
        if raw_layers in adataRaw.layers:
            if is_raw_counts(adataRaw, layer = raw_layers):
                logging.info(f"Using adataRaw.layers['{raw_layers}'] as raw counts data.")
                adataRaw.X = adataRaw.layers[raw_layers]
            else:
                raise ValueError("adataRaw.X is not raw counts data, and no suitable 'counts' layer found.")
        else:
            raise ValueError(f"adataRaw.X is not raw counts data, and no {raw_layers} layer found.")
    else:
        logging.info("adataRaw.X is confirmed as raw counts data.")

    # 检查注释列是否存在
    missing_cols = [key for key in annotation_keys if key not in adataAnn.obs.columns]
    if missing_cols:
        logging.error(f"Columns {missing_cols} not found in adataAnn.obs")
        raise KeyError(f"Columns {missing_cols} not found in adataAnn.obs")

    # 找到公共细胞
    common_cells = adataRaw.obs_names.intersection(adataAnn.obs_names)
    overlap_ratio = len(common_cells) / len(adataRaw.obs_names)

    if overlap_ratio < min_overlap_ratio:
        logging.error(
            f"Overlap between adataRaw and adataAnn is too small: "
            f"{overlap_ratio:.2%} < {min_overlap_ratio:.2%}"
        )
        raise ValueError(
            f"Overlap between adataRaw and adataAnn is too small: "
            f"{overlap_ratio:.2%} < {min_overlap_ratio:.2%}"
        )

    # 只保留有注释的细胞
    adataRaw = adataRaw[common_cells].copy()

    # 重新对齐 adataAnn
    adataAnn_sub = adataAnn[common_cells].copy()
    adataAnn_sub = adataAnn_sub[adataRaw.obs_names]

    # 添加注释
    for key in annotation_keys:
        adataRaw.obs[key] = adataAnn_sub.obs[key]

    logging.info(
        f"Annotations {annotation_keys} added for {len(common_cells)} cells "
        f"(after aligning indices, overlap ratio: {overlap_ratio:.2%})."
    )

    return adataRaw

def makePseudoBulk(
    adataRaw:ad.AnnData,
    adataAnn:ad.AnnData,
    fig_dir: str,
    out:str,
    fig_flag:bool=False,
    obs_to_keep:list = ['sample',"celltype"],
    raw_layers: str = "counts"
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
    adata = Pre_AddAnnotations2RawCountsAnndata(adataRaw,adataAnn,raw_layers=raw_layers)
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

    sc.pp.normalize_total(adata_pb, target_sum=1e4)
    sc.pp.log1p(adata_pb)
    sc.pp.pca(adata_pb)
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
    if fig_flag:
        logging.info(f"Plot PCA plot of column in PseduoBulk adata's obs")
        sc.pl.pca(adata_pb, color=adata_pb.obs, ncols=1, size=300)
        plt.savefig(f"{fig_dir}/{out}_obsPCA.png", dpi=300, bbox_inches='tight')
    
    logging.info(f"After normalizing and log1p; The max count of adata:{np.max(adata_pb.X)}")
    adata_pb.X = adata_pb.layers['counts'].copy()
    logging.info(f"After copying layers 'counts'; The max count of adata:{np.max(adata_pb.X)}")
    logging.info(f"PseduoBulk adata check complete !")
    logging.info(f"PseduoBulk adata make complete!")
    return adata_pb

def Pre_FeatureRename(
        adata:ad.AnnData,
        features:list
) -> ad.AnnData :
    """Rename features in an AnnData object by replacing specific characters.

    This function iterates through a list of features in the `adata.obs` DataFrame
    and replaces spaces (' '), plus signs ('+'), periods ('.'), and hyphens ('-')
    with underscores ('_') or removes them. This is useful for standardizing
    feature names to a valid format for R downstream analysis or visualization.

    Parameters
    ----------
    adata
        Annotated data matrix.
    features
        A list of strings representing the names of the features (columns)
        in `adata.obs` to be processed.

    Returns
    -------
    ad.AnnData
        The AnnData object with the specified feature names modified.
    
    Examples
    --------
    >>> import anndata as ad
    >>> import pandas as pd
    >>> obs_data = pd.DataFrame({
    ...     'sample name': ['A + B', 'C-D', 'E.F'],
    ...     'gene_expression': [1, 2, 3]
    ... })
    >>> adata = ad.AnnData(obs=obs_data)
    >>> features_to_clean = ['sample name']
    >>> new_adata = Pre_FeatureRename(adata, features_to_clean)
    >>> print(new_adata.obs)
      sample_name  gene_expression
    0       A_B            1
    1       C_D            2
    2       E_F            3
    """
    for feature in features:
        if feature in adata.obs.columns:
            adata.obs[feature] = adata.obs[feature].astype(str).str.replace(" ", "_")
            adata.obs[feature] = adata.obs[feature].astype(str).str.replace("+", "")
            adata.obs[feature] = adata.obs[feature].astype(str).str.replace(".", "_")
            adata.obs[feature] = adata.obs[feature].astype(str).str.replace("-", "_")
        else:
            logging.warning(f"Feature '{feature}' not found in adata.obs. Skipping.")
    
    return adata


def Pre_PseduoBulk(
    adata: ad.AnnData,
    cell_identity: str,
    obs_to_keep: list | None = None,  # additional metadata columns to keep, e.g., gender, age
    donor_key: str = "sample",
    cell_identity_key: str = "celltype",
    NUM_OF_CELL_PER_DONOR: int = 30,
    replicates_per_patient: int = 1
) -> ad.AnnData:
    """
    Aggregate single-cell expression data into pseudo-bulk samples per donor for a given cell type,
    optionally filtering donors with too few cells and keeping selected metadata.

    This function performs the following steps:
    1. Subset the AnnData object to the specified cell identity.
    2. Identify donors with sufficient number of cells (more than NUM_OF_CELL_PER_DONOR).
       Donors with fewer cells are dropped.
    3. For each donor, optionally create replicates by randomly splitting cells.
    4. Aggregate gene expression by summing counts for each gene per donor (or replicate).
       Additional metadata columns can be preserved by taking the first value.
    5. Return a new AnnData object with aggregated pseudo-bulk samples and selected metadata.

    Parameters
    ----------
    adata : AnnData
        The input AnnData object containing single-cell expression data.
    cell_identity : str
        The cell type to subset and aggregate.
    obs_to_keep : list of str, optional
        List of additional metadata columns from `.obs` to preserve in the aggregated data.
    donor_key : str, default "sample"
        Column name in `.obs` identifying donors.
    cell_identity_key : str, default "celltype"
        Column name in `.obs` containing cell type annotations.
    NUM_OF_CELL_PER_DONOR : int, default 30
        Minimum number of cells required per donor; donors with fewer cells are dropped.
    replicates_per_patient : int, default 1
        Number of pseudo-replicates to create per donor by random splitting.

    Returns
    -------
    AnnData
        A new AnnData object containing aggregated pseudo-bulk expression data for the specified cell type,
        with selected metadata columns preserved.
    # pd.DataFrame
    # A DataFrame of aggregated pseudo-bulk expression data, where each row corresponds
    # to a donor replicate and includes summed gene expression plus preserved metadata.
    """

    if obs_to_keep is None:
        obs_to_keep = []

    # Subset the AnnData to the given cell identity
    adata_cell_pop = adata[adata.obs[cell_identity_key] == cell_identity].copy()
    logging.info(f"Subsetted {adata_cell_pop.n_obs} cells of type '{cell_identity}'.")

    # Determine donors to drop based on NUM_OF_CELL_PER_DONOR
    size_by_donor = adata_cell_pop.obs.groupby([donor_key],observed=False).size()
    donors_to_drop = [
        donor
        for donor in size_by_donor.index
        if size_by_donor[donor] <= NUM_OF_CELL_PER_DONOR
    ]
    if len(donors_to_drop) > 0:
        logging.info(f"Dropping donors with insufficient cells: {donors_to_drop}")


    # Prepare empty DataFrame to store aggregated data
    df = pd.DataFrame(columns=[*adata_cell_pop.var_names, *obs_to_keep])

    adata_cell_pop.obs[donor_key] = adata_cell_pop.obs[donor_key].astype("category")

    # Loop over donors and aggregate
    for i, donor in enumerate(donors := adata_cell_pop.obs[donor_key].cat.categories):
        logging.info(f"\tProcessing donor {i+1} / {len(donors)}: {donor}")
        if donor not in donors_to_drop:
            adata_donor = adata_cell_pop[adata_cell_pop.obs[donor_key] == donor]

            # Create replicates by randomly splitting cells
            indices = list(adata_donor.obs_names)
            random.shuffle(indices)
            indices = np.array_split(np.array(indices), replicates_per_patient)

            for j, rep_idx in enumerate(indices):
                adata_replicate = adata_donor[rep_idx]

                # Define aggregation: sum gene expression, keep first value for additional metadata
                agg_dict = {gene: "sum" for gene in adata_replicate.var_names}
                for obs in obs_to_keep:
                    agg_dict[obs] = "first"

                # Create DataFrame for this replicate
                df_donor = pd.DataFrame(adata_replicate.X.toarray(),
                                        index=adata_replicate.obs_names,
                                        columns=adata_replicate.var_names)
                df_donor = df_donor.join(adata_replicate.obs[obs_to_keep])

                # Aggregate by donor
                df_donor = df_donor.groupby(donor_key,observed=False).agg(agg_dict)
                df_donor[donor_key] = donor
                df.loc[f"{donor}_{cell_identity}_{j}"] = df_donor.loc[donor]
    
    logging.info("Aggregation complete.")

    # Create new AnnData object from aggregated DataFrame
    adata_cell_pop = sc.AnnData(
        X=df[adata_cell_pop.var_names].values.astype(np.float32),
        obs=df.drop(columns=adata_cell_pop.var_names),
        var=adata_cell_pop.var.copy()
    )
    logging.info(f"Created aggregated AnnData with {adata_cell_pop.n_obs} pseudo-bulk samples.")
    return adata_cell_pop



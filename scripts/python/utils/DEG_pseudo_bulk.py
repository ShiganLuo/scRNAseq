import anndata as ad
import pandas as pd
import random
import scanpy as sc
import numpy as np
import logging
logging.basicConfig(
    level=logging.INFO, 
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
# logging.debug("This is a debug message.") 
# logging.info("This is an info message.")
# logging.warning("This is a warning message.")
# logging.error("This is an error message.")

def Pre_AddAnnotations2RawCountsAnndata(
    adataRaw: ad.AnnData,
    adataAnn: ad.AnnData,
    annotation_keys: list[str] = ["celltype"],
    min_overlap_ratio: float = 0.8  # Minimum fraction of overlapping cells
) -> ad.AnnData:
    """
    Copy multiple annotation columns from one AnnData object to another raw counts AnnData object.

    This function adds specified annotation columns (e.g., 'celltype', 'leiden') from 
    `adataAnn` to `adataRaw.obs`. It first checks if the cell indices (obs_names) of 
    the two AnnData objects are identical. If they match exactly, the annotations are 
    copied directly. Otherwise, the function aligns annotations using the intersection 
    of cell indices. Cells not present in `adataAnn` will have missing values (pd.NA) 
    in the new columns.

    If the overlap between `adataRaw` and `adataAnn` is smaller than `min_overlap_ratio`, 
    a ValueError is raised to prevent partial or incorrect annotation.

    Parameters
    ----------
    adataRaw : AnnData
        The AnnData object containing raw counts data.
    adataAnn : AnnData
        The AnnData object containing annotation information.
    annotation_keys : list of str, default ["celltype"]
        List of column names in `adataAnn.obs` to copy into `adataRaw.obs`.
    min_overlap_ratio : float, default 0.8
        Minimum fraction of cells that must overlap between `adataRaw` and `adataAnn`.
        If the overlap is smaller, a ValueError is raised.

    Returns
    -------
    AnnData
        The original `adataRaw` with new annotation columns added in `.obs`.

    Example
    -------
    adataRaw = PreAddAnnotations2RawCountsAnndata(
        adataRaw, adataAnnotated, annotation_keys=["celltype", "leiden"]
    )
    """
    # Check if all annotation columns exist
    missing_cols = [key for key in annotation_keys if key not in adataAnn.obs.columns]
    if missing_cols:
        logging.error(f"Columns {missing_cols} not found in adataAnn.obs")
        raise KeyError(f"Columns {missing_cols} not found in adataAnn.obs")

    # Fast path: indices match exactly
    if (adataRaw.obs_names == adataAnn.obs_names).all():
        for key in annotation_keys:
            adataRaw.obs[key] = adataAnn.obs[key]
        logging.info(f"Annotations {annotation_keys} added directly (indices match).")
    else:
        # Align by intersection
        common_cells = adataRaw.obs_names.intersection(adataAnn.obs_names)
        overlap_ratio = len(common_cells) / len(adataRaw.obs_names)
        if overlap_ratio < min_overlap_ratio:
            logging.info(f"Overlap between adataRaw and adataAnn is too small: {overlap_ratio:.2%} < {min_overlap_ratio:.2%}")
            raise ValueError(
                f"Overlap between adataRaw and adataAnn is too small: "
                f"{overlap_ratio:.2%} < {min_overlap_ratio:.2%}"
            )

        # Initialize new columns with missing values
        for key in annotation_keys:
            adataRaw.obs[key] = pd.NA

        # Assign values for overlapping cells
        for key in annotation_keys:
            adataRaw.obs.loc[common_cells, key] = adataAnn.obs.loc[common_cells, key]

        logging.info(f"Annotations {annotation_keys} added for {len(common_cells)} cells (after aligning indices).")

    return adataRaw

def Pre_AddAnnotations2RawCountsAnndata_FullGene(
    adataRaw: ad.AnnData,
    adataAnn: ad.AnnData,
    annotation_keys: list[str] = ["celltype"],
    min_overlap_ratio: float = 0.8
) -> ad.AnnData:
    """
    将 adataAnn 的注释信息注入 adataRaw，仅按细胞对齐；
    保留 adataRaw 的全部基因，删除不匹配细胞。
    对 varm、layers 形状不匹配项自动跳过。

    Parameters
    ----------
    adataRaw : AnnData
        原始计数数据（将被更新）
    adataAnn : AnnData
        含注释和分析结果的 AnnData
    annotation_keys : list of str
        要从 adataAnn.obs 复制的列
    min_overlap_ratio : float
        最小允许细胞重叠比例，低于该值将报错

    Returns
    -------
    AnnData
        含完整基因、并注入注释后的 AnnData
    """

    # 检查注释列是否存在
    missing = [k for k in annotation_keys if k not in adataAnn.obs.columns]
    if missing:
        raise KeyError(f"以下注释列在 adataAnn.obs 中不存在: {missing}")

    # 对齐细胞（仅取交集）
    common_cells = adataRaw.obs_names.intersection(adataAnn.obs_names)
    overlap_ratio = len(common_cells) / len(adataRaw.obs_names)
    if overlap_ratio < min_overlap_ratio:
        raise ValueError(f"细胞重叠比例过低: {overlap_ratio:.2%} < {min_overlap_ratio:.2%}")

    # 只保留公共细胞（保留全部基因）
    adataRaw = adataRaw[common_cells, :].copy()
    adataAnn = adataAnn[common_cells, :].copy()

    # === 复制 obs 注释 ===
    for key in annotation_keys:
        adataRaw.obs[key] = adataAnn.obs[key]

    # === 尝试复制 obsm / varm / uns / layers ===
    def safe_copy_attr(src, dest, attr_name, axis_check=None):
        """安全复制 obsm/varm/layers，形状不匹配直接跳过"""
        src_attr = getattr(src, attr_name)
        dest_attr = getattr(dest, attr_name)
        for key, val in src_attr.items():
            try:
                if axis_check == "obs" and val.shape[0] != dest.n_obs:
                    logging.warning(f"跳过 {attr_name}['{key}'] (obs 维度不匹配)")
                    continue
                if axis_check == "var" and val.shape[0] != dest.n_vars:
                    logging.warning(f"跳过 {attr_name}['{key}'] (var 维度不匹配)")
                    continue
                dest_attr[key] = val
            except Exception as e:
                logging.warning(f"复制 {attr_name}['{key}'] 时出错: {e}")

    safe_copy_attr(adataAnn, adataRaw, "obsm", axis_check="obs")
    safe_copy_attr(adataAnn, adataRaw, "varm", axis_check="var")
    safe_copy_attr(adataAnn, adataRaw, "layers", axis_check=None)

    # === 复制 uns（直接覆盖同名项）===
    for key, val in adataAnn.uns.items():
        adataRaw.uns[key] = val

    logging.info(f"成功注入注释 {annotation_keys}，保留 {adataRaw.n_vars} 个基因，{adataRaw.n_obs} 个细胞。")
    return adataRaw


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



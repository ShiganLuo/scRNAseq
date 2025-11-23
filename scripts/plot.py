import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import anndata as ad
import logging


def pieChart(
        adata: ad.AnnData,
        TARGET_CLUSTER:int,
        outfile:str,
        CLUSTER_KEY = 'leiden',
        SAMPLE_KEY = 'sample',
):
    
    adata_cluster_14 = adata[adata.obs[CLUSTER_KEY].astype(int) == TARGET_CLUSTER].copy()
    sample_counts = adata_cluster_14.obs[SAMPLE_KEY].value_counts()

    sample_proportions = sample_counts / sample_counts.sum()
    logging.info(sample_proportions)

    plt.figure(figsize=(8, 8))
    wedges, texts, autotexts = plt.pie(
        sample_proportions.values,
        labels=sample_proportions.index,
        autopct='%1.1f%%',  # 显示百分比
        startangle=90,      # 从顶部开始绘制
        wedgeprops={'edgecolor': 'black'}
    )

    plt.title(f'Cluster {TARGET_CLUSTER} Sample Composition Proportions', fontsize=16)
    plt.axis('equal')

    plt.savefig(outfile)
def pieChartForWhole(
                     adata_dict:dict[str,ad.AnnData],
                     outfile:str,
                     SAMPLE_KEY="sample",
):
    all_proportions = []
    all_sample_ids = []
    first_condition = list(adata_dict.keys())[0]
    reference_samples = set(adata_dict[first_condition].obs[SAMPLE_KEY].unique())
    for condition, adata in adata_dict.items():
        current_samples = set(adata.obs[SAMPLE_KEY].unique())
        if current_samples != reference_samples:
            print(f"Warning: Samples in condition '{condition}' do not match the reference samples ({first_condition}).")

        all_sample_ids.extend(list(current_samples))
    unique_samples = sorted(list(set(all_sample_ids)))
    for condition, adata in adata_dict.items():
        sample_counts = adata.obs[SAMPLE_KEY].value_counts()
        sample_proportions = sample_counts / sample_counts.sum()
        df_temp = sample_proportions.reset_index()
        df_temp.columns = [SAMPLE_KEY, 'proportion']
        df_temp['condition'] = condition
        
        all_proportions.append(df_temp)

    df_combined = pd.concat(all_proportions)
    df_combined = df_combined.astype({
    'sample': 'object',      # 将 category 转换为 object (string)
    'proportion': 'float64',
    'condition': 'object'
    })
    print(df_combined)
    print(df_combined.dtypes)
    pivot_df = df_combined.pivot_table(
        index='condition', 
        columns='sample', 
        values='proportion', 
        fill_value=0,
        aggfunc='sum'
    )
    pivot_df = pivot_df.reindex(['before qc', 'after qc'])
    print(pivot_df)
    fig, ax = plt.subplots(figsize=(10, 6))

    pivot_df.plot(
        kind='bar',
        stacked=True,
        ax=ax,
        colormap='tab20',
        edgecolor='black'
    )

    ax.set_title('Sample Contribution Proportion within Each Condition', fontsize=16)
    ax.set_ylabel('Proportion within Condition (%)', fontsize=12)
    ax.set_xlabel('Condition', fontsize=12)

    ax.set_ylim(0, 1) 

    plt.xticks(rotation=0)

    ax.legend(title='Sample ID', bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.tight_layout()
    fig.savefig(outfile)

if __name__ == "__main__":
    # adata = sc.read_h5ad("/disk5/luosg/scRNAseq/output/combine/Lung/Lung_annotate.h5ad")
    # pieChart(adata,14,"/disk5/luosg/scRNAseq/output/combine/Lung/result/Lung_14_pie.png")
    adata_dict = {
        "before qc": sc.read_h5ad("/disk5/luosg/scRNAseq/output/combine/Muscle/Muscle_deambiendRNA.h5ad"),
        "after qc": sc.read_h5ad("/disk5/luosg/scRNAseq/output/combine/Muscle/Muscle_qc.h5ad")
    } 
    pieChartForWhole(adata_dict,"/disk5/luosg/scRNAseq/output/combine/Muscle/result/qc_ba_sampe_proportion.png")
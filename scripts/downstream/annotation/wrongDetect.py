import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import silhouette_samples
from sklearn.neighbors import NearestNeighbors
import anndata as ad
def detectWrongCluster(
        adata: ad.AnnData,
        score_dist:int = 0.01,
        score_purity:int = 0.09,
        score_sil:int = 0.9


):
    # ==== 使用 PCA 或 Harmony 空间 ====
    # if 'X_pca_harmony' in adata.obsm.keys():
    #     X = adata.obsm['X_pca_harmony'][:, :30]
    # else:
    #     X = adata.obsm['X_pca'][:, :30]
    X = adata.obsm['X_umap']
    clusters = adata.obs['celltype'].astype(str)

    # ==== 1️⃣ cluster 中心距离（z-score）====
    centroids = {}
    for c in clusters.unique():
        centroids[c] = X[clusters == c].mean(axis=0)

    dist = np.array([np.linalg.norm(X[i] - centroids[clusters[i]]) for i in range(X.shape[0])]) # 欧式距离
    zscore = np.zeros_like(dist)
    for c in clusters.unique():
        idx = np.where(clusters == c)[0]
        z = (dist[idx] - dist[idx].mean()) / (dist[idx].std(ddof=1) + 1e-9) # z-score
        zscore[idx] = z

    adata.obs['center_dist'] = dist
    adata.obs['center_dist_z'] = zscore

    # ==== 2️⃣ kNN purity（邻居同簇比例）====
    k = 30
    nbrs = NearestNeighbors(n_neighbors=k+1).fit(X)
    distances, indices = nbrs.kneighbors(X)
    purity = []
    for i in range(indices.shape[0]):
        neigh = indices[i, 1:]
        same = np.sum(clusters.values[neigh] == clusters.values[i])
        purity.append(same / k)
    adata.obs[f'knn_purity_k{k}'] = purity

    # ==== 3️⃣ silhouette coefficient ====
    sil = silhouette_samples(X, clusters.values)
    adata.obs['silhouette'] = sil

    # ==== 4️⃣ 综合指标 miscluster_score ====
    def minmax(x):
        x = np.array(x, dtype=float)
        return (x - np.nanmin(x)) / (np.nanmax(x) - np.nanmin(x) + 1e-9)

    adata.obs['score_dist'] = minmax(np.clip(adata.obs['center_dist_z'], -5, 5))
    adata.obs['score_sil'] = 1 - minmax(np.clip(adata.obs['silhouette'], -1, 1))
    adata.obs['score_purity'] = 1 - minmax(adata.obs[f'knn_purity_k{k}'])

    adata.obs['miscluster_score'] = (
        score_dist * adata.obs['score_dist'] + 
        score_purity * adata.obs['score_purity'] +
        score_sil * adata.obs['score_sil']
    )

    # 阈值：0.6 通常较稳健，可调整
    adata.obs['is_suspect'] = adata.obs['miscluster_score'] > 0.49

    # ==== 5️⃣ 输出结果 ====
    suspects = adata.obs[adata.obs['is_suspect']]
    print(f"可疑细胞数量: {suspects.shape[0]} / {adata.n_obs}")
    suspects[['sample','experiment','leiden','center_dist_z',f'knn_purity_k{k}','silhouette','miscluster_score']].to_csv('suspect_cells.csv')

    # ==== 6️⃣ 可视化 ====
    fig = sc.pl.umap(
        adata,
        color=['leiden', 'miscluster_score', 'is_suspect', 'pct_counts_mt'],
        size=20,
        cmap='viridis',
        return_fig=True
    )
    return fig

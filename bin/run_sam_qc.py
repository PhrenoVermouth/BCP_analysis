#!/usr/bin/env python

import argparse
from typing import Dict, List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sps
from samalg import SAM
from scipy import stats

def vectorized_bs_mean_diff(
    X_target, X_next, alpha: float = 0.1, n_bs: int = 500, rng: np.random.Generator | None = None
):
    """
    Update: Optimized Bootstrap using Matrix Multiplication (Weighted Sums).
    Avoids memory copying of the large expression matrix.
    """

    # 1. Preprocessing
    target = _to_dense_matrix(X_target)
    next_cluster = _to_dense_matrix(X_next)

    n_genes = target.shape[1] if target.size else next_cluster.shape[1]
    if target.size == 0 or next_cluster.size == 0 or n_genes == 0:
        zeros = np.zeros(n_genes)
        return zeros, zeros, zeros

    rng = rng or np.random.default_rng()
    n_cells_t = target.shape[0]
    n_cells_n = next_cluster.shape[0]

    # 2. Generate weight matrix (Weights Matrix)
    # Instead of "moving" the data, generate a count matrix of shape (n_bs, n_cells).
    
    # Target weight
    w_t = np.zeros((n_bs, n_cells_t), dtype=np.float32)
    for i in range(n_bs):
        # current_datetime = datetime.now()
        # print('Finish n_bs ' + str(i))
        # print(current_datetime)
        indices = rng.integers(0, n_cells_t, size=n_cells_t)
        # Count the occurrences of each index.
        w_t[i, :] = np.bincount(indices, minlength=n_cells_t)

    # Next Cluster weight
    w_n = np.zeros((n_bs, n_cells_n), dtype=np.float32)
    for i in range(n_bs):
        indices = rng.integers(0, n_cells_n, size=n_cells_n)
        w_n[i, :] = np.bincount(indices, minlength=n_cells_n)

    # 3. Core: Matrix Multiplication
    # Bootstrap mean = (Weight matrix x Original expression matrix) / Total number of cells 
    # dimensions: (n_bs, n_cells) @ (n_cells, n_genes) -> (n_bs, n_genes)
    means_target_bs = (w_t @ target) / n_cells_t
    means_next_bs = (w_n @ next_cluster) / n_cells_n

    # 4. Calculate differences and confidence intervals.
    diffs = means_target_bs - means_next_bs
    real_diff = target.mean(axis=0) - next_cluster.mean(axis=0)
    
    lower = np.percentile(diffs, (alpha / 2) * 100, axis=0)
    upper = np.percentile(diffs, (1 - alpha / 2) * 100, axis=0)
    
    return real_diff, lower, upper

# def vectorized_bs_mean_diff(
#     X_target, X_next, alpha: float = 0.1, n_bs: int = 500, rng: np.random.Generator | None = None
# ):
#     """Bootstrap mean differences for all genes at once using vectorized resampling."""
#     target = _to_dense_matrix(X_target)
#     next_cluster = _to_dense_matrix(X_next)

#     n_genes = target.shape[1] if target.size else next_cluster.shape[1]
#     if target.size == 0 or next_cluster.size == 0 or n_genes == 0:
#         zeros = np.zeros(n_genes)
#         return zeros, zeros, zeros

#     rng = rng or np.random.default_rng()
#     n_cells_t, n_cells_n = target.shape[0], next_cluster.shape[0]
#     means_target_bs = np.empty((n_genes, n_bs))
#     means_next_bs = np.empty_like(means_target_bs)

#     for i in range(n_bs):
#         idx_t = rng.integers(0, n_cells_t, size=n_cells_t)
#         idx_n = rng.integers(0, n_cells_n, size=n_cells_n)
#         means_target_bs[:, i] = target[idx_t].mean(axis=0)
#         means_next_bs[:, i] = next_cluster[idx_n].mean(axis=0)

#     diffs = means_target_bs - means_next_bs
#     real_diff = target.mean(axis=0) - next_cluster.mean(axis=0)
#     lower = np.percentile(diffs, (alpha / 2) * 100, axis=1)
#     upper = np.percentile(diffs, (1 - alpha / 2) * 100, axis=1)
#     return real_diff, lower, upper

def _to_dense_matrix(array_like) -> np.ndarray:
    if sps.issparse(array_like):
        return array_like.toarray()
    return np.asarray(array_like)


def _sort_cluster_labels(labels: List[str]) -> List[str]:
    def _sort_key(label: str):
        try:
            return int(label)
        except ValueError:
            return label

    return sorted(labels, key=_sort_key)


def _ensure_marker_dict(markers, clusters: List[str]) -> Dict[str, List[str]]:
    marker_dict: Dict[str, List[str]] = {}
    if isinstance(markers, dict):
        for cluster in clusters:
            values = markers.get(cluster)
            if values is None:
                values = markers.get(str(cluster), [])
            marker_dict[str(cluster)] = list(values)
    else:
        for idx, cluster in enumerate(clusters):
            marker_dict[str(cluster)] = list(markers[idx]) if len(markers) > idx else []
    return marker_dict

def run_qc2(
    sample_id,
    input_h5ad,
    whitelist,
    n_hvg,
    n_neighbors,
    n_pcs,
    marker_genes_file=None,
):
    """Run single-cell data QC after SoupX correction using precomputed whitelist."""
    adata = sc.read_h5ad(input_h5ad)
    adata.var_names_make_unique()
    sc.pp.normalize_total(adata, target_sum=1e6)
    sc.pp.log1p(adata)

    keep_cells = pd.read_csv(whitelist, header=None)[0].tolist()
    adata = adata[adata.obs_names.isin(keep_cells), :]

    sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg)
    sc.pp.pca(adata, n_comps=n_pcs)
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)

    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    adata.obs.rename(columns={'total_counts': 'n_counts', 'n_genes_by_counts': 'n_genes'}, inplace=True)

    adata.raw = adata.copy()
    sam = SAM(adata)
    sam.run(
        projection='umap',
        preprocessing='StandardScaler',
        k=n_neighbors,
        npcs=n_pcs,
        n_genes=n_hvg,
    )
    sam.leiden_clustering(res=0.5)

    cluster_labels = sam.adata.obs['leiden_clusters'].astype(str)
    clusters_sorted = _sort_cluster_labels(cluster_labels.unique().tolist())

    cluster_stats = (
        sam.adata.obs.assign(leiden_clusters=cluster_labels)
        .groupby('leiden_clusters')[['n_counts', 'n_genes']]
        .median()
        .loc[clusters_sorted]
    )

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw={'width_ratios': [3, 1]})
    
    sc.pl.umap(sam.adata,color=['leiden_clusters'],legend_loc='on data',s = 40, ax=axes[0],show=False, title=f'{sample_id} - Leiden Clustering (UMAP) - QC2')

    y_pos = np.arange(len(clusters_sorted))
    bar_height = 0.2
    axes[1].barh(
        y_pos - bar_height / 2,
        cluster_stats['n_counts'],
        height=bar_height,
        label='n_counts'
    )
    axes[1].barh(
        y_pos + bar_height / 2,
        cluster_stats['n_genes'],
        height=bar_height,
        label='n_genes'
    )
    axes[1].set_yticks(y_pos)
    axes[1].set_yticklabels(clusters_sorted)
    axes[1].invert_yaxis()
    axes[1].set_xlabel('Median value')
    axes[1].set_title('Median n_counts and n_genes by Leiden cluster')
    axes[1].legend()

    plt.tight_layout()
    plt.savefig(
        f'4.{sample_id}_umap_leiden_QC2_mqc.png', dpi=300, bbox_inches='tight'
    )
    plt.close(fig)


    ranked_genes = list(sam.adata.uns['ranked_genes'])
    strategy_one_genes = ranked_genes[:500]
    strategy_one_markers = {
        cluster: strategy_one_genes for cluster in clusters_sorted
    }

    X = sam.adata.X
    labels = cluster_labels.values
    n_genes_total = sam.adata.n_vars
    cluster_means = np.vstack(
        [np.asarray(X[labels == c].mean(axis=0)).ravel() for c in clusters_sorted]
    )
    second_best_cluster_idx = (
        np.argsort(cluster_means, axis=0)[-2, :]
        if len(clusters_sorted) > 1
        else np.zeros(n_genes_total, dtype=int)
    )

    rng = np.random.default_rng()

    strategy_two_results: Dict[str, pd.DataFrame] = {}
    for idx, cluster in enumerate(clusters_sorted):
        y_mask = labels == cluster
        a = np.zeros(n_genes_total)
        b = np.zeros(n_genes_total)
        c_ci = np.zeros(n_genes_total)
        p_values = np.ones(n_genes_total)

        target_data = _to_dense_matrix(X[y_mask, :])
        rival_cluster_data: Dict[int, np.ndarray] = {}
        bs_results: Dict[int, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}

        for rival_idx in np.unique(second_best_cluster_idx):
            rival_mask = labels == clusters_sorted[rival_idx]
            rival_data = _to_dense_matrix(X[rival_mask, :])
            rival_cluster_data[rival_idx] = rival_data
            bs_results[rival_idx] = vectorized_bs_mean_diff(
                target_data, rival_data, alpha=0.1, n_bs=500, rng=rng
            )

        for j in range(n_genes_total):
            rival_idx = second_best_cluster_idx[j]
            rival_data = rival_cluster_data[rival_idx]

            real_diff, lower, upper = bs_results[rival_idx]
            a[j] = real_diff[j]
            b[j] = lower[j]
            c_ci[j] = upper[j]

            x_target = target_data[:, j]
            x_next = rival_data[:, j]
            if x_target.size > 0 and x_next.size > 0:
                _, p_values[j] = stats.mannwhitneyu(
                    x_target, x_next, alternative="two-sided"
                )

        anchor = pd.DataFrame(
            {
                'genes': sam.adata.var_names,
                'a': a,
                'b': b,
                'c': c_ci,
                'p': p_values,
            }
        ).set_index('genes')
        strategy_two_results[cluster] = anchor

    strategy_two_markers: Dict[str, List[str]] = {}
    for cluster, anchor in strategy_two_results.items():
        ranked = anchor.sort_values(by=['p', 'b'], ascending=[True, False])
        strategy_two_results[cluster] = ranked
        strategy_two_markers[cluster] = ranked.index.tolist()

    intersection_markers: Dict[str, List[str]] = {}
    strategy_one_gene_set = set(strategy_one_genes)
    for cluster in clusters_sorted:
        secondary = strategy_two_markers[cluster]
        intersection = [g for g in secondary if g in strategy_one_gene_set]
        intersection_markers[cluster] = intersection

    intersection_rows = []
    for cluster in clusters_sorted:
        anchor = strategy_two_results[cluster]
        cluster_intersection = intersection_markers[cluster]
        for rank, gene in enumerate(cluster_intersection, start=1):
            intersection_rows.append(
                {
                    'cluster': cluster,
                    'strategy': 'intersection',
                    'rank': rank,
                    'gene': gene,
                    'p_value': anchor.loc[gene, 'p'],
                    'b_value': anchor.loc[gene, 'b'],
                }
            )

    if intersection_rows:
        markers_df = pd.DataFrame(intersection_rows)
        markers_df.to_csv(f'{sample_id}_markers.csv', index=False)

    sam.adata.uns['marker_genes'] = {
        'strategy_one': strategy_one_markers,
        'strategy_two': {
            k: v_df.sort_values(by='b', ascending=False)
            .reset_index()
            .to_dict(orient='list')
            for k, v_df in strategy_two_results.items()
        },
        'strategy_two_top': strategy_two_markers,
        'intersection': intersection_markers,
    }

    genes_for_dotplot: List[str] = []
    for cluster in clusters_sorted:
        candidate_genes = intersection_markers[cluster][:3]
        for gene in candidate_genes:
            if gene not in genes_for_dotplot:
                genes_for_dotplot.append(gene)

    if genes_for_dotplot:
        sc.pl.dotplot(
            sam.adata,
            genes_for_dotplot,
            'leiden_clusters',
            standard_scale='var',
            dendrogram=True,
            show=False,
        )
        plt.savefig(
            f'5.{sample_id}_marker_genes_dotplot_mqc.png',
            dpi=300,
            bbox_inches='tight'
        )
        plt.close()

    raw_mat = pd.DataFrame(
        sam.adata.raw.X.todense(),
        index=sam.adata.obs_names,
        columns=sam.adata.raw.var_names
    )
    raw_mat.to_csv(f'{sample_id}_raw_counts.csv')
    sam.adata.obs['leiden_clusters'].to_csv(f'{sample_id}_leiden_clusters.csv', header=True)
    sam.adata.write_h5ad(f'{sample_id}_filtered_QC2.h5ad')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run SAM QC on SoupX output")
    parser.add_argument('--sample_id', required=True, help="Sample identifier")
    parser.add_argument('--input_h5ad', required=True, help="SoupX corrected h5ad")
    parser.add_argument('--whitelist', required=True, help="Path to whitelist of non-doublet cells")
    parser.add_argument('--n_hvg', type=int, default=3000, help="Number of highly variable genes")
    parser.add_argument('--n_neighbors', type=int, default=15, help="Number of nearest neighbors in PCA space")
    parser.add_argument('--n_pcs', type=int, default=50, help="Number of principal components")
    parser.add_argument('--marker_genes', type=str, default=None, help="Optional path to marker genes list (one gene per line)")
    args = parser.parse_args()

    run_qc2(
        args.sample_id,
        args.input_h5ad,
        args.whitelist,
        args.n_hvg,
        args.n_neighbors,
        args.n_pcs,
        args.marker_genes,
    )

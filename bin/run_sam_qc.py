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


def bs_mean_diff(x: np.ndarray, y: np.ndarray, alpha: float = 0.1, n_bs: int = 500):
    """Bootstrap the mean difference between two vectors."""
    if x.size == 0 or y.size == 0:
        return 0.0, 0.0, 0.0
    rng = np.random.default_rng()
    diffs = np.empty(n_bs)
    for i in range(n_bs):
        diffs[i] = rng.choice(x, size=x.size, replace=True).mean() - rng.choice(
            y, size=y.size, replace=True
        ).mean()

    mean_diff = x.mean() - y.mean()
    lower = np.percentile(diffs, (alpha / 2) * 100)
    upper = np.percentile(diffs, (1 - alpha / 2) * 100)
    return mean_diff, lower, upper


def _to_dense_vector(array_like) -> np.ndarray:
    if sps.issparse(array_like):
        return array_like.toarray().ravel()
    return np.asarray(array_like).ravel()


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

    sc.pl.umap(
        sam.adata,
        color='leiden_clusters',
        show=False,
        title=f'{sample_id} - Leiden Clustering (UMAP) - QC2'
    )
    plt.savefig(f'4.{sample_id}_umap_leiden_QC2.png', dpi=300, bbox_inches='tight')
    plt.close()

    cluster_labels = sam.adata.obs['leiden_clusters'].astype(str)
    clusters_sorted = _sort_cluster_labels(cluster_labels.unique().tolist())

    strategy_one_markers_raw, strategy_one_scores = sam.identify_marker_genes_rf(
        labels='leiden_clusters', clusters=clusters_sorted, n_genes=10
    )
    strategy_one_markers = _ensure_marker_dict(strategy_one_markers_raw, clusters_sorted)
    for cluster in clusters_sorted:
        strategy_one_markers[cluster] = strategy_one_markers[cluster][:10]

    X = sam.adata.X
    labels = cluster_labels.values
    n_genes_total = sam.adata.n_vars
    cluster_means = np.vstack(
        [np.asarray(X[labels == c].mean(axis=0)).ravel() for c in clusters_sorted]
    )

    strategy_two_results: Dict[str, pd.DataFrame] = {}
    for idx, cluster in enumerate(clusters_sorted):
        y_mask = labels == cluster
        a = np.zeros(n_genes_total)
        b = np.zeros(n_genes_total)
        c_ci = np.zeros(n_genes_total)
        p_values = np.ones(n_genes_total)
        for j in range(n_genes_total):
            gene_means = cluster_means[:, j]
            order = np.argsort(gene_means)
            second_idx = order[-2] if len(clusters_sorted) > 1 else idx
            y_next_mask = labels == clusters_sorted[second_idx]

            x_target = _to_dense_vector(X[y_mask, j])
            x_next = _to_dense_vector(X[y_next_mask, j])

            a[j], b[j], c_ci[j] = bs_mean_diff(x_target, x_next, alpha=0.1, n_bs=500)
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
        top_genes = anchor.sort_values(by='b', ascending=False).head(10).index.tolist()
        strategy_two_markers[cluster] = top_genes

    intersection_markers: Dict[str, List[str]] = {}
    for cluster in clusters_sorted:
        primary = strategy_one_markers[cluster]
        secondary = strategy_two_markers[cluster]
        intersection = [g for g in secondary if g in primary][:10]
        intersection_markers[cluster] = intersection

    marker_rows = []
    for cluster in clusters_sorted:
        for rank, gene in enumerate(strategy_one_markers[cluster], start=1):
            marker_rows.append(
                {
                    'cluster': cluster,
                    'strategy': 'strategy_one',
                    'rank': rank,
                    'gene': gene,
                }
            )
        for rank, gene in enumerate(strategy_two_markers[cluster], start=1):
            marker_rows.append(
                {
                    'cluster': cluster,
                    'strategy': 'strategy_two',
                    'rank': rank,
                    'gene': gene,
                }
            )
        for rank, gene in enumerate(intersection_markers[cluster], start=1):
            marker_rows.append(
                {
                    'cluster': cluster,
                    'strategy': 'intersection',
                    'rank': rank,
                    'gene': gene,
                }
            )

    markers_df = pd.DataFrame(marker_rows)
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
        if len(candidate_genes) < 3:
            candidate_genes = strategy_one_markers[cluster][:3]
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
            f'5.{sample_id}_marker_genes_dotplot.png', dpi=300, bbox_inches='tight'
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
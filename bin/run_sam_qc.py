#!/usr/bin/env python

import scanpy as sc
import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
from samalg import SAM

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

    sc.tl.rank_genes_groups(sam.adata, groupby='leiden_clusters', method='wilcoxon')
    rank_df = sc.get.rank_genes_groups_df(sam.adata, None)
    top_markers = rank_df.groupby('group').head(2)['names'].tolist()
    external_markers = []
    if marker_genes_file:
        with open(marker_genes_file) as f:
            external_markers = [line.strip() for line in f if line.strip()]
    genes_to_plot = list(dict.fromkeys(top_markers + external_markers))
    genes_to_plot = [g for g in genes_to_plot if g in sam.adata.raw.var_names]
    if genes_to_plot:
        expr = sam.adata.raw.to_adata()[:, genes_to_plot].to_df()
        expr['cluster'] = sam.adata.obs['leiden_clusters'].values
        mean_expr = expr.groupby('cluster').mean()
        mean_expr = mean_expr.reindex(sorted(mean_expr.index, key=lambda x: int(x)))
        g = sns.clustermap(
            mean_expr,
            cmap='viridis',
            figsize=(0.5 * len(genes_to_plot) + 5, 0.5 * mean_expr.shape[0] + 5),
        )
        g.fig.suptitle(f'{sample_id} Marker Gene Expression')
        g.ax_heatmap.set_xlabel('Gene')
        g.ax_heatmap.set_ylabel('Cluster')
        g.savefig(f'5.{sample_id}_marker_genes_heatmap.png')
        plt.close(g.fig)

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

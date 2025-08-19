#!/usr/bin/env python

import scanpy as sc
import scrublet as scr
import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
from samalg import SAM

def run_qc2(
    sample_id,
    input_h5ad,
    min_genes,
    min_cells,
    max_genes,
    max_counts,
    max_mito,
    mito_prefixes,
    mito_gene_list,
    n_hvg,
    n_neighbors,
    n_pcs,
    marker_genes_file=None,
):
    """
    Run single-cell data QC.
    adata - after removing ambient RNA, filter lower limits
    adata_QC1 - mitochondria filtration
    adata_QC2 - doublet removal
    """
    # 1. Load data
    adata = sc.read_h5ad(input_h5ad)
    adata.var_names_make_unique()
    adata.var_names_make_unique()
    sc.pp.normalize_total(adata, target_sum=1e6) #log2 CPM, but input has been processed by SoupX
    sc.pp.log1p(adata)


    # 2. --- Initial filtering and QC calculations ---
    # a. Basic filtering on lower limit
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    print("---mito list---")
    print(mito_gene_list)
    # b. Calculate mitochondrial gene percentage
    if mito_gene_list:
        with open(mito_gene_list) as f:
            mito_genes = {line.strip() for line in f if line.strip()}
        adata.var['mito'] = adata.var_names.isin(mito_genes)
    else:
        adata.var['mito'] = adata.var_names.str.startswith(tuple(mito_prefixes))
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=['mito'], percent_top=None, log1p=False, inplace=True
    )

    adata_QC1 = adata[adata.obs.pct_counts_mito < max_mito * 100, :] # Scanpy calculates this as a percentage


    # c. Filter doublets
    # @Prateek I didn't filter the upper limit here as I think they would be duplicated with function of SCR
    # adata = adata[adata.obs.n_genes_by_counts < max_genes, :]
    # adata = adata[adata.obs.total_counts < max_counts, :]
    raw_matrix = adata_QC1.X.todense()
    scrub = scr.Scrublet(raw_matrix)
    doublet_score, predicted_doublets = scrub.scrub_doublets()
    real_cells = np.logical_not(predicted_doublets)
    adata_QC1.obs['predicted_doublet'] = predicted_doublets

    # # Backup: Cell number dependent doublet expectation (https://kb.10xgenomics.com/hc/en-us/articles/360054599512-What-is-the-cell-multiplet-rate-when-using-the-3-CellPlex-Kit-for-Cell-Multiplexing)
    # rate = adata_QC1.n_obs / 1000 * 0.008
    # scrub = scr.Scrublet(adata_QC1.X, expected_doublet_rate=rate, random_state=0)
    # adata_QC1.obs['doublet_score'], predicted_doublets = scrub.scrub_doublets()
    # scrub.plot_histogram()
    # plt.show()
    # auto_threshold = scrub.threshold_
    # final_threshold = auto_threshold
    # adata_QC1.obs['predicted_doublet'] = adata_QC1.obs['doublet_score'] > final_threshold
    #
    adata_QC2 = adata_QC1[real_cells,:]

    # e. Identify highly variable genes and compute PCA/neighbor graph
    sc.pp.highly_variable_genes(adata_QC2, n_top_genes=n_hvg)
    sc.pp.pca(adata_QC2, n_comps=n_pcs)
    sc.pp.neighbors(adata_QC2, n_neighbors=n_neighbors, n_pcs=n_pcs)

    # d. Create QC plots
    # Fig 0: Mitochodria Removal
    fig0, axes = plt.subplots(1, 2, figsize=(15, 5))
    sc.pl.violin(adata, 'pct_counts_mito', jitter=0.4, ax=axes[0], show=False)
    sc.pl.violin(adata_QC1, 'pct_counts_mito', jitter=0.4, ax=axes[1], show=False)
    fig0.suptitle(f'{sample_id} - Before Filtering - QC1')
    fig0.tight_layout()
    fig0.savefig(f'1.{sample_id}_violin_before_mito_filtering_QC1.png')
    plt.close(fig0)

    # Fig 1: Doublet cells shown on scatter plots
    adata_QC1.obs['predicted_doublet'] = adata_QC1.obs['predicted_doublet'].astype('category')
    custom_palette = ['#DDDDDD', 'red']
    fig1 = plt.figure()
    sc.pl.scatter(
        adata_QC1,  # as marks saved within QC1
        x='total_counts',
        y='n_genes_by_counts',
        color='predicted_doublet',
        palette=custom_palette,
        show=False,
        title=f'{sample_id} - Predicted Doublets Highlighted - QC2'
    )
    fig1.savefig(f'2.{sample_id}_scatter_doublet_highlight_QC2.png')
    plt.close(fig1)

    # Fig 2: Global violin before and after doublet removal
    fig2, axes = plt.subplots(2, 2, figsize=(18, 10))
    fig2.suptitle(f'{sample_id} - QC Metrics: Before vs. After Doublet Removal', fontsize=16)

    axes[0, 0].set_title('Genes per Cell (Before)')
    sc.pl.violin(adata_QC1, 'n_genes_by_counts', jitter=0.4, ax=axes[0, 0], show=False)
    axes[0, 1].set_title('Counts per Cell (Before)')
    sc.pl.violin(adata_QC1, 'total_counts', jitter=0.4, ax=axes[0, 1], show=False)

    axes[1, 0].set_title('Genes per Cell (After)')
    sc.pl.violin(adata_QC2, 'n_genes_by_counts', jitter=0.4, ax=axes[1, 0], show=False)
    axes[1, 1].set_title('Counts per Cell (After)')
    sc.pl.violin(adata_QC2, 'total_counts', jitter=0.4, ax=axes[1, 1], show=False)

    fig2.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig2.savefig(f'3.{sample_id}_violin_comparison_QC2.png')
    plt.close(fig2)


    # e. Save stats for multiqc
    # --- Save pre-filtering statistics for MultiQC use ---
    number_cells = adata.n_obs
    median_genes = np.median(adata.obs['n_genes_by_counts'])
    median_counts = np.median(adata.obs['total_counts'])

    number_cells_QC1 = adata_QC1.n_obs
    median_genes_QC1 = np.median(adata_QC1.obs['n_genes_by_counts'])
    median_counts_QC1 = np.median(adata_QC1.obs['total_counts'])

    number_cells_QC2 = adata_QC2.n_obs
    median_genes_QC2  = np.median(adata_QC2.obs['n_genes_by_counts'])
    median_counts_QC2 = np.median(adata_QC2.obs['total_counts'])


    # 3. CLUSTERING using SAM algorithm
    # ----------------
    # store raw (pre‐normalized) counts in .raw for safe‐keeping
    adata_QC2.raw = adata_QC2.copy()
    # The SAM library no longer accepts ``n_genes`` during initialization.
    # Instead, instantiate ``SAM`` with the AnnData object and pass the
    # desired number of genes and principal components to ``run`` via the
    # ``n_genes`` and ``npcs`` arguments. Using the old ``d`` parameter
    # caused a ``TypeError`` because it is not a valid keyword argument.
    sam = SAM(adata_QC2)
    sam.run(
        projection='umap',
        preprocessing='StandardScaler',
        k=n_neighbors,
        npcs=n_pcs,
        n_genes=n_hvg,
    )
    sam.leiden_clustering(res=0.5) # As a QC pipeline, here the resolution is relatively low

    # Fig 3: UMAP for Leiden clusters generated by SAM
    sc.pl.umap(
        sam.adata,
        color='leiden_clusters',
        show=False,
        title=f'{sample_id} - Leiden Clustering (UMAP) - QC2'
    )
    plt.savefig(f'4.{sample_id}_umap_leiden_QC2.png')
    plt.close()

    # Identify marker genes and plot heatmap
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

    # 4. EXPORT RAW COUNTS MATRIX
    # ----------------------------
    # pull it out as a dense DataFrame
    raw_mat = pd.DataFrame(
        sam.adata.raw.X.todense(),
        index=sam.adata.obs_names,
        columns=sam.adata.raw.var_names
    )
    raw_mat.to_csv(f'{sample_id}_raw_counts.csv')
    sam.adata.obs['leiden_clusters'].to_csv(f'{sample_id}_leiden_clusters.csv', header=True)

    # --- Generate a simple tab-delimited file for MultiQC ---
    # with open(f'{sample_id}_qc_metrics.tsv', 'w') as f:
    #     f.write("metric\tvalue\n")
    #     f.write(f"cells_initial\t{number_cells}\n")
    #     f.write(f"median_genes_initial\t{int(median_genes)}\n")
    #     f.write(f"median_counts_initial\t{int(median_counts)}\n")
    #     f.write(f"cells_after_MT_Removal\t{number_cells_QC1}\n")
    #     f.write(f"median_genes_MT_Removal\t{int(median_genes_QC1)}\n")
    #     f.write(f"median_counts_MT_Removal\t{int(median_counts_QC1)}\n")
    #     f.write(f"cells_after_Doublet_Removal\t{number_cells_QC2}\n")
    #     f.write(f"median_genes_Doublet_Removal\t{int(median_genes_QC2)}\n")
    #     f.write(f"median_counts_Doublet_Removal\t{int(median_counts_QC2)}\n")

    with open(f"{sample_id}_cells.tsv", "w") as f:
        f.write("# plot_type: 'table'\n")
        f.write("# section_name: 'Cells QC Metrics'\n")
        f.write("# description: 'Cell counts at different filtering steps'\n")
        f.write("# pconfig:\n")
        f.write("#     sortRows: false\n")
        f.write("Sample\tcells_initial\tcells_after_MT_Removal\tcells_after_Doublet_Removal\n")
        f.write(f"{sample_id}\t{number_cells}\t{number_cells_QC1}\t{number_cells_QC2}\n")

    with open(f"{sample_id}_counts.tsv", "w") as f:
        f.write("# plot_type: 'table'\n")
        f.write("# section_name: 'Counts QC Metrics'\n")
        f.write("# description: 'Median UMI counts at different filtering steps'\n")
        f.write("# pconfig:\n")
        f.write("#     sortRows: false\n")
        f.write("Sample\tmedian_counts_initial\tmedian_counts_MT_Removal\tmedian_counts_Doublet_Removal\n")
        f.write(f"{sample_id}\t{int(median_counts)}\t{int(median_counts_QC1)}\t{int(median_counts_QC2)}\n")

    with open(f"{sample_id}_genes.tsv", "w") as f:
        f.write("# plot_type: 'table'\n")
        f.write("# section_name: 'Genes QC Metrics'\n")
        f.write("# description: 'Median genes at different filtering steps'\n")
        f.write("# pconfig:\n")
        f.write("#     sortRows: false\n")
        f.write("Sample\tmedian_genes_initial\tmedian_genes_MT_Removal\tmedian_genes_Doublet_Removal\n")
        f.write(f"{sample_id}\t{int(median_genes)}\t{int(median_genes_QC1)}\t{int(median_genes_QC2)}\n")

    # Save processed AnnData object
    sam.adata.write_h5ad(f'{sample_id}_filtered_QC2.h5ad')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run SAM QC on STARsolo output")
    parser.add_argument('--sample_id', required=True, help="Sample identifier")
    parser.add_argument('--input_h5ad', required=True, help="Sample h5ad")

    #parser.add_argument('--matrix_dir', required=True, help="Path to STARsolo matrix directory")
    parser.add_argument('--min_genes', type=int, default=600, help="Min genes per cell")
    parser.add_argument('--min_cells', type=int, default=3, help="Min cells per gene")
    parser.add_argument('--max_genes', type=int, default=6000, help="Max genes per cell")
    parser.add_argument('--max_counts', type=int, default=20000, help="Max UMI counts per cell")
    parser.add_argument('--max_mito', type=float, default=0.2, help="Max mitochondrial percentage (as a fraction, e.g., 0.2 for 20%)")
    parser.add_argument('--mito_prefixes', nargs='+', default=['mt-'], help="Prefix(es) for mitochondrial genes")
    parser.add_argument('--mito_gene_list', type=str, default=None, help="Path to text file listing mitochondrial genes (one per line)")
    parser.add_argument('--n_hvg', type=int, default=3000, help="Number of highly variable genes")
    parser.add_argument('--n_neighbors', type=int, default=15, help="Number of nearest neighbors in PCA space")
    parser.add_argument('--n_pcs', type=int, default=50, help="Number of principal components")
    parser.add_argument('--marker_genes', type=str, default=None, help="Optional path to marker genes list (one gene per line)")

    args = parser.parse_args()

    run_qc2(
        args.sample_id,
        args.input_h5ad,
        args.min_genes,
        args.min_cells,
        args.max_genes,
        args.max_counts,
        args.max_mito,
        args.mito_prefixes,
        args.mito_gene_list,
        args.n_hvg,
        args.n_neighbors,
        args.n_pcs,
        args.marker_genes,
    )

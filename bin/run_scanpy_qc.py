#!/usr/bin/env python

import scanpy as sc
import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt

def run_qc2(sample_id, matrix_dir, min_genes, min_cells, max_genes, max_counts, max_mito, mito_prefixes):
    """
    Run single-cell data QC.
    adata - only lower limit filtration
    adata_QC1 - mitochondria filtration
    adata_QC2 - doublet removal
    """
    # 1. Load data
    adata = sc.read_10x_mtx(matrix_dir, var_names='gene_symbols', cache=True)
    adata.var_names_make_unique()
    
    # 2. --- Initial filtering and QC calculations ---
    # a. Basic filtering on lower limit
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)

    # b. Calculate mitochondrial gene percentage
    adata.var['mito'] = adata.var_names.str.startswith(tuple(mito_prefixes))
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mito'], percent_top=None, log1p=False, inplace=True)
  
    adata_QC1 = adata[adata.obs.pct_counts_mito < max_mito * 100, :] # Scanpy calculates this as a percentage

    
    # c. Filter doublets 
    # @Prateek I didn't filter the upper limit here as I think they would be duplicated with function of SCR
    # adata = adata[adata.obs.n_genes_by_counts < max_genes, :]
    # adata = adata[adata.obs.total_counts < max_counts, :]
    raw_matrix = adata_QC1.X.todense()
    scrub = scr.Scrublet(raw_matrix)
    doublet_score, predicted_doublets = scrub.scrub_doublets()
    real_cells = np.logical_not(predicted_doublets)
    adata_QC2 = adata_QC1[real_cells,:]
    
    
    # d. Create QC plots 
    # Fig 0: Mitochodria Removal
    fig0, axes = plt.subplots(1, 2, figsize=(15, 5))
    sc.pl.violin(adata, 'pct_counts_mito', jitter=0.4, ax=axes[0], show=False)
    sc.pl.violin(adata_QC1, 'pct_counts_mito', jitter=0.4, ax=axes[1], show=False)
    fig0.suptitle(f'{sample_id} - Before Filtering - QC1')
    fig0.tight_layout()
    fig0.savefig(f'{sample_id}_violin_before_filtering_QC1.png')
    plt.close(fig0)
    
    # Fig 1: Doublet cells shown on scatter plots
    fig1 = plt.figure()
    sc.pl.scatter(
        adata, 
        x='total_counts', 
        y='n_genes_by_counts', 
        color='predicted_doublet',  
        show=False,
        title=f'{sample_id} - Predicted Doublets Highlighted - QC2'
    )
    fig1.savefig(f'{sample_id}_scatter_doublet_highlight_QC2.png')
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
    fig2.savefig(f'{sample_id}_violin_comparison_QC2.png')
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

    
    # 3. CLUSTERING @Prateek potential flexibilities for parameters? I don't scale often (annotated), what do you think?
    # ----------------
    # store raw (pre‐normalized) counts in .raw for safe‐keeping
    adata_QC2.raw = adata_QC2.copy()
    # (a) Normalize & log-transform
    sc.pp.normalize_total(adata_QC2, target_sum=1e6)
    sc.pp.log1p(adata_QC2)    
    sc.pp.highly_variable_genes(adata_QC2, n_top_genes=2000, subset=True)
    # sc.pp.scale(adata_QC2, max_value=10)
    sc.tl.pca(adata_QC2, svd_solver='arpack')
    sc.pp.neighbors(adata_QC2, n_neighbors=10, n_pcs=40)
    sc.tl.leiden(adata_QC2, resolution=1.0, key_added='leiden_clusters')


    # Fig 3: UMAP fot Leiden
    fig3 = plt.figure()
    sc.pl.umap(
        adata_QC2,
        color='leiden_clusters',
        show=False,
        title=f'{sample_id} - Leiden Clustering (UMAP) - QC2'
    )
    fig3.savefig(f'{sample_id}_umap_leiden_QC2.png')
    plt.close(fig3)
    
    # 4. EXPORT RAW COUNTS MATRIX
    # ----------------------------
    # pull it out as a dense DataFrame
    raw_mat = pd.DataFrame(
        adata_QC2.raw.X.todense(),
        index=adata_QC2.obs_names,
        columns=adata_QC2.raw.var_names
    )
    raw_mat.to_csv(f'{sample_id}_raw_counts.csv')
    adata_QC2.obs['leiden_clusters'].to_csv(f'{sample_id}_leiden_clusters.csv', header=True)
         
    # --- Generate a simple tab-delimited file for MultiQC ---
    with open(f'{sample_id}_qc_metrics.tsv', 'w') as f:
        f.write("metric\tvalue\n")
        f.write(f"cells_initial\t{number_cells}\n")
        f.write(f"median_genes_initial\t{int(median_genes)}\n")
        f.write(f"median_counts_initial\t{int(median_counts)}\n")
        f.write(f"cells_after_MT_Removal\t{number_cells_QC1}\n")
        f.write(f"median_genes_MT_Removal\t{int(median_genes_QC1)}\n")
        f.write(f"median_counts_MT_Removal\t{int(median_counts_QC1)}\n")
        f.write(f"cells_after_Doublet_Removal\t{number_cells_QC2}\n")
        f.write(f"median_genes_Doublet_Removal\t{int(median_genes_QC2)}\n")
        f.write(f"median_counts_Doublet_Removal\t{int(median_counts_QC2)}\n")
        
    # @Prateek: Do we need to keep more h5ads above?
    adata.write_h5ad(f'{sample_id}_filtered_QC2.h5ad')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run Scanpy QC on STARsolo output")
    parser.add_argument('--sample_id', required=True, help="Sample identifier")
    parser.add_argument('--matrix_dir', required=True, help="Path to STARsolo matrix directory")
    parser.add_argument('--min_genes', type=int, default=600, help="Min genes per cell")
    parser.add_argument('--min_cells', type=int, default=3, help="Min cells per gene")
    parser.add_argument('--max_genes', type=int, default=6000, help="Max genes per cell")
    parser.add_argument('--max_counts', type=int, default=20000, help="Max UMI counts per cell")
    parser.add_argument('--max_mito', type=float, default=0.2, help="Max mitochondrial percentage (as a fraction, e.g., 0.2 for 20%)")
    parser.add_argument('--mito_prefixes', nargs='+', default=['mt-'], help="Prefix(es) for mitochondrial genes")

    args = parser.parse_args()

    run_qc2(args.sample_id, args.matrix_dir, args.min_genes, args.min_cells, args.max_genes, args.max_counts, args.max_mito, args.mito_prefixes)

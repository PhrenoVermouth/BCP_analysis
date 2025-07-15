#!/usr/bin/env python

import scanpy as sc
import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt

def run_qc(sample_id, matrix_dir, min_genes, min_cells, max_genes, max_counts, max_mito, mito_prefixes):
    """
    Run single-cell data QC.
    """
    # 1. Load data
    adata = sc.read_10x_mtx(matrix_dir, var_names='gene_symbols', cache=True)
    adata.var_names_make_unique()

    # --- Initial filtering and QC calculations ---
    # a. Basic filtering
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)

    # b. Calculate mitochondrial gene percentage
    # Note: The prefix 'mt-' is a common default; may need to adjust based on your species.
    adata.var['mito'] = adata.var_names.str.startswith(tuple(mito_prefixes))
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mito'], percent_top=None, log1p=False, inplace=True)

    # c. Create QC plots (before filtering)
    fig1, axes = plt.subplots(1, 3, figsize=(15, 5))
    sc.pl.violin(adata, 'n_genes_by_counts', jitter=0.4, ax=axes[0], show=False)
    sc.pl.violin(adata, 'total_counts', jitter=0.4, ax=axes[1], show=False)
    sc.pl.violin(adata, 'pct_counts_mito', jitter=0.4, ax=axes[2], show=False)
    fig1.suptitle(f'{sample_id} - Before Filtering')
    fig1.tight_layout()
    fig1.savefig(f'{sample_id}_violin_before_filtering.png')
    plt.close(fig1)

    fig2 = plt.figure()
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', color='pct_counts_mito', show=False)
    plt.title(f'{sample_id} - Before Filtering')
    fig2.savefig(f'{sample_id}_scatter_before_filtering.png')
    plt.close(fig2)

    # --- Save pre-filtering statistics for MultiQC use ---
    initial_cells = adata.n_obs
    median_genes = np.median(adata.obs['n_genes_by_counts'])
    median_counts = np.median(adata.obs['total_counts'])
    
    # --- Perform filtering ---
    adata = adata[adata.obs.n_genes_by_counts < max_genes, :]
    adata = adata[adata.obs.total_counts < max_counts, :]
    adata = adata[adata.obs.pct_counts_mito < max_mito * 100, :] # Scanpy calculates this as a percentage

    # --- Save post-filtering statistics ---
    final_cells = adata.n_obs
    
    # --- Generate a simple tab-delimited file for MultiQC ---
    with open(f'{sample_id}_qc_metrics.tsv', 'w') as f:
        f.write("metric\tvalue\n")
        f.write(f"initial_cells\t{initial_cells}\n")
        f.write(f"median_genes_per_cell\t{int(median_genes)}\n")
        f.write(f"median_counts_per_cell\t{int(median_counts)}\n")
        f.write(f"final_cells\t{final_cells}\n")

    # --- Save results ---
    # Save filtered AnnData file
    adata.write_h5ad(f'{sample_id}_filtered.h5ad')

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

    run_qc(args.sample_id, args.matrix_dir, args.min_genes, args.min_cells, args.max_genes, args.max_counts, args.max_mito, args.mito_prefixes)

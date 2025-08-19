#!/usr/bin/env python

import scanpy as sc
import scrublet as scr
import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt
import seaborn as sns

def run_scrublet(sample_id, matrix_dir, min_genes, min_cells, max_mito, mito_prefixes, mito_gene_list=None):
    adata = sc.read_10x_mtx(matrix_dir, var_names='gene_symbols')
    adata.var_names_make_unique()

    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)

    if mito_gene_list:
        with open(mito_gene_list) as f:
            mito_genes = {line.strip() for line in f if line.strip()}
        adata.var['mito'] = adata.var_names.isin(mito_genes)
    else:
        adata.var['mito'] = adata.var_names.str.startswith(tuple(mito_prefixes))
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mito'], percent_top=None, log1p=False, inplace=True)

    adata_QC1 = adata[adata.obs.pct_counts_mito < max_mito, :]

    raw_matrix = adata_QC1.X.todense()
    scrub = scr.Scrublet(raw_matrix)
    doublet_score, predicted_doublets = scrub.scrub_doublets()
    adata_QC1.obs['predicted_doublet'] = predicted_doublets
    adata_QC2 = adata_QC1[~predicted_doublets, :]

    # Fig 0: Mitochondria removal
    fig0, axes = plt.subplots(1, 2, figsize=(15, 5))
    sc.pl.violin(adata, 'pct_counts_mito', jitter=0.4, ax=axes[0], show=False)
    sc.pl.violin(adata_QC1, 'pct_counts_mito', jitter=0.4, ax=axes[1], show=False)
    fig0.suptitle(f'{sample_id} - Mito Filtering - QC1')
    fig0.tight_layout()
    fig0.savefig(f'1.{sample_id}_violin_mito_filtering_QC1.png')
    plt.close(fig0)

    # Fig 1: Doublet scatter
    adata_QC1.obs['predicted_doublet'] = adata_QC1.obs['predicted_doublet'].astype('category')
    custom_palette = ['#DDDDDD', 'red']
    fig1, ax = plt.subplots()
    sc.pl.scatter(
        adata_QC1,
        x='total_counts',
        y='n_genes_by_counts',
        color='predicted_doublet',
        palette=custom_palette,
        ax=ax,
        show=False,
        title=f'{sample_id} - Predicted Doublets Highlighted - QC2'
    )
    fig1.savefig(
        f'2.{sample_id}_scatter_doublet_highlight_QC2.png',
        dpi=300,
        bbox_inches='tight')
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

    # Save whitelist
    adata_QC2.obs_names.to_series().to_csv(f'{sample_id}_whitelist.txt', index=False, header=False)

    # Save metrics for MultiQC
    number_cells = adata.n_obs
    median_genes = int(np.median(adata.obs['n_genes_by_counts']))
    median_counts = int(np.median(adata.obs['total_counts']))

    number_cells_QC1 = adata_QC1.n_obs
    median_genes_QC1 = int(np.median(adata_QC1.obs['n_genes_by_counts']))
    median_counts_QC1 = int(np.median(adata_QC1.obs['total_counts']))

    number_cells_QC2 = adata_QC2.n_obs
    median_genes_QC2 = int(np.median(adata_QC2.obs['n_genes_by_counts']))
    median_counts_QC2 = int(np.median(adata_QC2.obs['total_counts']))

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
        f.write(f"{sample_id}\t{median_counts}\t{median_counts_QC1}\t{median_counts_QC2}\n")

    with open(f"{sample_id}_genes.tsv", "w") as f:
        f.write("# plot_type: 'table'\n")
        f.write("# section_name: 'Genes QC Metrics'\n")
        f.write("# description: 'Median genes at different filtering steps'\n")
        f.write("# pconfig:\n")
        f.write("#     sortRows: false\n")
        f.write("Sample\tmedian_genes_initial\tmedian_genes_MT_Removal\tmedian_genes_Doublet_Removal\n")
        f.write(f"{sample_id}\t{median_genes}\t{median_genes_QC1}\t{median_genes_QC2}\n")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run Scrublet doublet detection")
    parser.add_argument('--sample_id', required=True, help="Sample identifier")
    parser.add_argument('--matrix_dir', required=True, help="Path to raw count matrix directory")
    parser.add_argument('--min_genes', type=int, default=600, help="Min genes per cell")
    parser.add_argument('--min_cells', type=int, default=3, help="Min cells per gene")
    parser.add_argument('--max_mito', type=float, default=0.2, help="Max mitochondrial percentage (fraction)")
    parser.add_argument('--mito_prefixes', nargs='+', default=['mt-'], help="Prefix(es) for mitochondrial genes")
    parser.add_argument('--mito_gene_list', type=str, default=None, help="Path to text file listing mitochondrial genes (one per line)")
    args = parser.parse_args()

    run_scrublet(
        args.sample_id,
        args.matrix_dir,
        args.min_genes,
        args.min_cells,
        args.max_mito,
        args.mito_prefixes,
        args.mito_gene_list,
    )

#!/usr/bin/env python

import argparse
import os
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scrublet as scr

def _load_summary(summary_csv: Optional[str]):
    estimated_cells = None
    sequencing_saturation = ""
    if summary_csv and os.path.exists(summary_csv):
        try:
            summary_df = pd.read_csv(
                summary_csv,
                header=None,
                names=["metric", "value"],
            )
            est_cells = summary_df.loc[
                summary_df["metric"] == "Estimated Number of Cells", "value"
            ]
            if not est_cells.empty:
                estimated_cells = int(float(est_cells.iloc[0]))

            saturation = summary_df.loc[
                summary_df["metric"] == "Sequencing Saturation", "value"
            ]
            if not saturation.empty:
                sequencing_saturation = round(float(saturation.iloc[0]), 4)
        except Exception:
            estimated_cells = None
            sequencing_saturation = ""
    return estimated_cells, sequencing_saturation


def _plot_knee(knee, expected_num_cells: Optional[int], sample_id: str):
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.loglog(knee, range(len(knee)), linewidth=5, color="g")

    if expected_num_cells is not None and expected_num_cells < len(knee):
        ax.axvline(x=knee[expected_num_cells], linewidth=3, color="k")
        ax.axhline(y=expected_num_cells, linewidth=3, color="k")

    ax.set_xlabel("UMI Counts")
    ax.set_ylabel("Set of Barcodes")
    plt.grid(True, which="both")
    fig.savefig(
        f"0.{sample_id}_knee_plot_mqc.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)


def _apply_mingene_filter(adata, min_gene_threshold: int):
    filtered = adata[adata.obs["n_genes_by_counts"] >= min_gene_threshold, :]
    if filtered.n_obs < 5000 and min_gene_threshold > 200:
        min_gene_threshold = 200
        filtered = adata[adata.obs["n_genes_by_counts"] >= min_gene_threshold, :]
    return filtered, min_gene_threshold


def run_scrublet(
    sample_id,
    matrix_dir,
    min_genes,
    min_cells,
    max_mito,
    mito_prefixes,
    mito_gene_list=None,
    summary_csv=None,
):
    adata = sc.read_10x_mtx(matrix_dir, var_names='gene_symbols')
    adata.var_names_make_unique()

    sc.pp.filter_genes(adata, min_cells=min_cells)

    if mito_gene_list:
        with open(mito_gene_list) as f:
            mito_genes = {line.strip() for line in f if line.strip()}
        adata.var['mito'] = adata.var_names.isin(mito_genes)
    else:
        adata.var['mito'] = adata.var_names.str.startswith(tuple(mito_prefixes))
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mito'], percent_top=None, log1p=False, inplace=True)

    knee = np.sort((np.array(adata.X.sum(axis=1))).flatten())[::-1]
    expected_num_cells, sequencing_saturation = _load_summary(summary_csv)
    _plot_knee(knee, expected_num_cells, sample_id)

    adata.write_h5ad(f"{sample_id}_initial.h5ad")

    raw_matrix = adata.X.todense()
    scrub = scr.Scrublet(raw_matrix, expected_doublet_rate=0.06)
    doublet_score, predicted_doublets = scrub.scrub_doublets(
        min_counts=2,
        min_cells=3,
        min_gene_variability_pctl=85,
        n_prin_comps=30,
    )
    adata.obs['predicted_doublet'] = predicted_doublets
    adata_QC1 = adata[~predicted_doublets, :]

    adata_QC_min, min_gene_threshold = _apply_mingene_filter(adata_QC1, min_genes)
    adata_QC2 = adata_QC_min[adata_QC_min.obs.pct_counts_mito < max_mito, :]

    # Fig 0: Doublet scatter
    adata.obs['predicted_doublet'] = adata.obs['predicted_doublet'].astype('category')
    custom_palette = ['#DDDDDD', 'red']
    fig0, ax = plt.subplots()
    sc.pl.scatter(
        adata,
        x='total_counts',
        y='n_genes_by_counts',
        color='predicted_doublet',
        palette=custom_palette,
        ax=ax,
        show=False,
        title=f'{sample_id} - Predicted Doublets Highlighted - QC1',
    )
    fig0.savefig(
        f'1.{sample_id}_scatter_doublet_highlight_QC1_mqc.png',
        dpi=300,
        bbox_inches='tight')
    plt.close(fig0)

    # Fig 1: Global violin before and after doublet removal
    fig1, axes = plt.subplots(2, 2, figsize=(18, 10))
    fig1.suptitle(f'{sample_id} - QC Metrics: Before vs. After Doublet Removal', fontsize=16)

    axes[0, 0].set_title('Genes per Cell (Before)')
    sc.pl.violin(adata, 'n_genes_by_counts', jitter=0.4, ax=axes[0, 0], show=False)
    axes[0, 1].set_title('Counts per Cell (Before)')
    sc.pl.violin(adata, 'total_counts', jitter=0.4, ax=axes[0, 1], show=False)

    axes[1, 0].set_title('Genes per Cell (After Doublet Removal)')
    sc.pl.violin(adata_QC1, 'n_genes_by_counts', jitter=0.4, ax=axes[1, 0], show=False)
    axes[1, 1].set_title('Counts per Cell (After Doublet Removal)')
    sc.pl.violin(adata_QC1, 'total_counts', jitter=0.4, ax=axes[1, 1], show=False)

    fig1.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig1.savefig(f'2.{sample_id}_violin_comparison_QC1_mqc.png')
    plt.close(fig1)

    # Fig 2: Mitochondria removal
    fig2, axes = plt.subplots(1, 2, figsize=(15, 5))
    sc.pl.violin(adata_QC_min, 'pct_counts_mito', jitter=0.4, ax=axes[0], show=False)
    sc.pl.violin(adata_QC2, 'pct_counts_mito', jitter=0.4, ax=axes[1], show=False)
    fig2.suptitle(f'{sample_id} - Mito Filtering - QC2 (min_genes={min_gene_threshold})')
    fig2.tight_layout()
    fig2.savefig(f'3.{sample_id}_violin_mito_filtering_QC2_mqc.png')
    plt.close(fig2)

    # Save whitelist
    adata_QC2.obs_names.to_series().to_csv(f'{sample_id}_whitelist.txt', index=False, header=False)

    # Save metrics for MultiQC
    number_cells = adata.n_obs
    median_genes = int(np.median(adata.obs['n_genes_by_counts']))
    median_counts = int(np.median(adata.obs['total_counts']))

    number_cells_QC1 = adata_QC1.n_obs

    number_cells_QC_min = adata_QC_min.n_obs
    median_genes_QC_min = int(np.median(adata_QC_min.obs['n_genes_by_counts']))
    median_counts_QC_min = int(np.median(adata_QC_min.obs['total_counts']))

    number_cells_QC2 = adata_QC2.n_obs
    median_genes_QC2 = int(np.median(adata_QC2.obs['n_genes_by_counts']))
    median_counts_QC2 = int(np.median(adata_QC2.obs['total_counts']))

    doublet_fraction = 0
    if number_cells > 0:
        doublet_fraction = round((number_cells - number_cells_QC1) / number_cells, 4)

    with open(f"{sample_id}_total_metaqc_partial.tsv", "w") as f:
        f.write("# plot_type: 'table'\n")
        f.write("# section_name: 'Total_metaQC'\n")
        f.write("# description: 'Combined QC metrics before ambient RNA correction'\n")
        f.write("# pconfig:\n")
        f.write("#     sortRows: false\n")
        f.write(
            "Sample\tInitial_n_cells\tInitial_median_genes\tInitial_median_counts\tDoublet_fraction\tQC_db_min_n_cells\tQC_db_min_median_genes\tQC_db_min_median_counts\tQC_db_min_mt_n_cells\tQC_db_min_mt_median_genes\tQC_db_min_mt_median_counts\tSequencing_saturation\tRho\n"
        )
        f.write(
            f"{sample_id}\t{number_cells}\t{median_genes}\t{median_counts}\t{doublet_fraction}\t{number_cells_QC_min}\t{median_genes_QC_min}\t{median_counts_QC_min}\t{number_cells_QC2}\t{median_genes_QC2}\t{median_counts_QC2}\t{sequencing_saturation}\t\n"
        )

    adata_QC1.write_h5ad(f"{sample_id}_rmdoublet.h5ad")
    adata_QC_min.write_h5ad(f"{sample_id}_rmdoublet_mingene{min_gene_threshold}.h5ad")
    adata_QC2.write_h5ad(f"{sample_id}_rmMT{max_mito}.h5ad")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run Scrublet doublet detection")
    parser.add_argument('--sample_id', required=True, help="Sample identifier")
    parser.add_argument('--matrix_dir', required=True, help="Path to raw count matrix directory")
    parser.add_argument('--min_genes', type=int, default=500, help="Min genes per cell before relaxation")
    parser.add_argument('--min_cells', type=int, default=3, help="Min cells per gene")
    parser.add_argument('--max_mito', type=float, default=0.2, help="Max mitochondrial percentage (fraction)")
    parser.add_argument('--mito_prefixes', nargs='+', default=['mt-'], help="Prefix(es) for mitochondrial genes")
    parser.add_argument('--mito_gene_list', type=str, default=None, help="Path to text file listing mitochondrial genes (one per line)")
    parser.add_argument('--summary_csv', type=str, default=None, help="Path to STARsolo Summary.csv for sequencing saturation")
    args = parser.parse_args()

    run_scrublet(
        args.sample_id,
        args.matrix_dir,
        args.min_genes,
        args.min_cells,
        args.max_mito,
        args.mito_prefixes,
        args.mito_gene_list,
        args.summary_csv,
    )

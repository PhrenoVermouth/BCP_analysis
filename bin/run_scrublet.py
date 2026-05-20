#!/usr/bin/env python

import argparse
import os
import sys
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scrublet as scr

# Fallback rule: if STAR's filtered matrix has more than TRIGGER cells, the
# barcode list is likely contaminated by empty droplets. Downstream then keeps
# only the top TARGET cells ranked by UMI counts. Surfaced via [FALLBACK] log
# lines and read back by bin/nf_debug_agent.py.
CELL_COUNT_FALLBACK_TRIGGER = 70000
CELL_COUNT_FALLBACK_TARGET = 40000

# Hard-coded threshold used when scrublet's auto-detected threshold > 0.4.
SCRUBLET_THRESHOLD_FALLBACK = 0.4


def _log_fallback(message: str) -> None:
    """Emit a parseable fallback marker. Picked up by nf_debug_agent.py."""
    print(f"[FALLBACK] {message}", flush=True)
    sys.stderr.write(f"[FALLBACK] {message}\n")
    sys.stderr.flush()


def _calculate_expected_doublet_rate(expected_num_cells: Optional[int]) -> float:
    """Scale expected doublet rate with expected_num_cells using piecewise linear rules."""
    if expected_num_cells is None:
        return 0.06
    if expected_num_cells <= 0:
        return round(0.01, 2)
    if expected_num_cells <= 20000:
        return round(0.01 + (0.08 - 0.01) * (expected_num_cells / 20000), 2)
    if expected_num_cells <= 40000:
        return round(0.08 + (0.2 - 0.08) * ((expected_num_cells - 20000) / 20000), 2)
    return round(0.2, 2)

    # if expected_num_cells <= 0:
    #     return 0.01
    # if expected_num_cells <= 20000:
    #     return 0.01 + (0.08 - 0.01) * (expected_num_cells / 20000)
    # if expected_num_cells <= 40000:
    #     return 0.08 + (0.2 - 0.08) * ((expected_num_cells - 20000) / 20000)
    # return 0.2

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


def _plot_knee(
    knee,
    expected_num_cells: Optional[int],
    sample_id: str,
    fallback_n_cells: Optional[int] = None,
):
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.loglog(knee, range(len(knee)), linewidth=5, color="g")

    legend_needed = False
    if expected_num_cells is not None and expected_num_cells < len(knee):
        ax.axvline(
            x=knee[expected_num_cells],
            linewidth=3,
            color="k",
            label=f"STARsolo expected = {expected_num_cells}",
        )
        ax.axhline(y=expected_num_cells, linewidth=3, color="k")
        legend_needed = True

    if fallback_n_cells is not None and fallback_n_cells < len(knee):
        ax.axvline(
            x=knee[fallback_n_cells],
            linewidth=2,
            color="red",
            linestyle="--",
            label=f"top-N fallback = {fallback_n_cells}",
        )
        ax.axhline(
            y=fallback_n_cells,
            linewidth=2,
            color="red",
            linestyle="--",
        )
        legend_needed = True

    if legend_needed:
        ax.legend(loc="lower left", fontsize=9)

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


def _load_matrix(matrix_dir):
    adata = sc.read_10x_mtx(matrix_dir, var_names='gene_symbols')
    adata.var_names_make_unique()
    return adata


def _save_scrublet_debug_plots(scrub, sample_id: str, tag: str = "debug"):
    """Save scrublet histogram and raw score CSV for debugging.
    Called immediately after scrub_doublets() so plots exist even if the run later fails.
    """
    prefix = f"DEBUG_{sample_id}_doublet_histogram_{tag}"

    # --- Save raw doublet scores to CSV for manual threshold selection ---
    try:
        scores = scrub.doublet_scores_obs_
        if scores is not None:
            pd.Series(scores, name="doublet_score").to_csv(
                f"{prefix}_scores.csv", index=False
            )
    except Exception:
        pass  # Scores may not exist if scrub_doublets crashed early

    # --- Save histogram plot ---
    try:
        histogram_plot = scrub.plot_histogram()
        histogram_figure = None
        if hasattr(histogram_plot, "get_figure"):
            histogram_figure = histogram_plot.get_figure()
        elif hasattr(histogram_plot, "savefig"):
            histogram_figure = histogram_plot
        elif isinstance(histogram_plot, (list, tuple)):
            for plot_obj in histogram_plot:
                if hasattr(plot_obj, "get_figure"):
                    histogram_figure = plot_obj.get_figure()
                    break
                if hasattr(plot_obj, "savefig"):
                    histogram_figure = plot_obj
                    break
        if histogram_figure is None:
            histogram_figure = plt.gcf()

        histogram_figure.savefig(f"{prefix}.png", dpi=200, bbox_inches="tight")
        plt.close(histogram_figure)
    except Exception as plot_err:
        # Fallback: manually plot observed + simulated side by side for threshold judgment
        try:
            obs = scrub.doublet_scores_obs_
            sim = scrub.doublet_scores_sim_

            fig, axes = plt.subplots(1, 2, figsize=(12, 4))
            fig.suptitle(f"{sample_id} — scrublet debug ({tag})", fontsize=12)

            # Left: observed on log-y scale (same as scrublet's default)
            if obs is not None and len(obs) > 0:
                axes[0].hist(obs, bins=50, color="steelblue", edgecolor="white")
                axes[0].set_yscale("log")
                axes[0].set_xlabel("Doublet Score")
                axes[0].set_ylabel("Cell Count (log)")
                axes[0].set_title("Observed transcriptomes")
                # Mark candidate thresholds for reference
                for t in [0.2, 0.3, 0.4]:
                    axes[0].axvline(t, color="red", linestyle="--", alpha=0.5, linewidth=1, label=f"t={t}")
                axes[0].legend(fontsize=8)

            # Right: simulated as probability density
            if sim is not None and len(sim) > 0:
                axes[1].hist(sim, bins=50, color="darkorange", edgecolor="white", density=True)
                axes[1].set_xlabel("Doublet Score")
                axes[1].set_ylabel("Prob. density")
                axes[1].set_title("Simulated doublets")
                sim_median = float(np.median(sim))
                axes[1].axvline(sim_median, color="black", linestyle="-", linewidth=1.5,
                                label=f"sim median={sim_median:.2f}")
                axes[1].legend(fontsize=8)

            fig.tight_layout()
            fig.savefig(f"{prefix}_fallback.png", dpi=200, bbox_inches="tight")
            plt.close(fig)
        except Exception:
            pass  # Nothing more we can do


def run_scrublet(
    sample_id,
    matrix_dir,
    min_genes,
    min_cells,
    max_mito,
    mito_prefixes,
    mito_gene_list=None,
    summary_csv=None,
    knee_matrix_dir=None,
    manual_threshold=None,
):
    adata = _load_matrix(matrix_dir)

    # Compute the knee curve from the unfiltered matrix. When no separate
    # --knee_matrix_dir was given, reuse the just-loaded adata BEFORE
    # filter_genes mutates it; loading the same .mtx twice would be wasted I/O
    # and produce an identical knee curve since filter_genes only drops genes.
    if knee_matrix_dir:
        knee_data = _load_matrix(knee_matrix_dir)
        knee = np.sort(np.array(knee_data.X.sum(axis=1)).flatten())[::-1]
    else:
        knee = np.sort(np.array(adata.X.sum(axis=1)).flatten())[::-1]

    sc.pp.filter_genes(adata, min_cells=min_cells)

    if mito_gene_list:
        with open(mito_gene_list) as f:
            mito_genes = {line.strip() for line in f if line.strip()}
        adata.var['mito'] = adata.var_names.isin(mito_genes)
    else:
        adata.var['mito'] = adata.var_names.str.startswith(tuple(mito_prefixes))
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mito'], percent_top=None, log1p=False, inplace=True)
    expected_num_cells, sequencing_saturation = _load_summary(summary_csv)

    # Decide cell-count fallback before drawing the knee plot so the red
    # dashed reference line is only added when the fallback actually applies.
    n_cells_initial = int(adata.n_obs)
    cell_count_fallback_applied = n_cells_initial > CELL_COUNT_FALLBACK_TRIGGER

    _plot_knee(
        knee,
        expected_num_cells,
        sample_id,
        fallback_n_cells=CELL_COUNT_FALLBACK_TARGET if cell_count_fallback_applied else None,
    )

    adata.write_h5ad(f"{sample_id}_initial.h5ad")

    if cell_count_fallback_applied:
        # Keep the top-N cells by total UMI counts (computed by calculate_qc_metrics).
        umi_counts = adata.obs["total_counts"].to_numpy()
        top_local = np.argsort(umi_counts)[::-1][:CELL_COUNT_FALLBACK_TARGET]
        top_local.sort()  # preserve original barcode order for stability
        kept_min_umi = float(umi_counts[top_local].min())
        _log_fallback(
            f"type=cell_count sample={sample_id} "
            f"original_cells={n_cells_initial} keep_top={CELL_COUNT_FALLBACK_TARGET} "
            f"min_umi_kept={kept_min_umi:.0f} "
            f"reason='STAR filtered cells > {CELL_COUNT_FALLBACK_TRIGGER}'"
        )
        adata = adata[top_local, :].copy()

######################################

    adata_QC_min, min_gene_threshold = _apply_mingene_filter(adata, min_genes)

    raw_matrix = adata_QC_min.X.todense()
    expected_doublet_rate = _calculate_expected_doublet_rate(adata_QC_min.n_obs) #Update 260101
    #expected_doublet_rate = _calculate_expected_doublet_rate(expected_num_cells)
    scrub = scr.Scrublet(raw_matrix, expected_doublet_rate=expected_doublet_rate)

    # --- Run scrub_doublets, but save histogram + scores before raising any error ---
    try:
        result = scrub.scrub_doublets(
            min_counts=2,
            min_cells=3,
            min_gene_variability_pctl=85,
            n_prin_comps=30,
        )
    except Exception as e:
        # Even if scrub_doublets itself throws, try to save whatever diagnostic we can
        _save_scrublet_debug_plots(scrub, sample_id, tag="crash")
        raise RuntimeError(f"scrub_doublets raised an exception: {e}") from e

    # Always save the histogram immediately after scrub_doublets, before any validation
    _save_scrublet_debug_plots(scrub, sample_id, tag="initial")

    if result is None or not isinstance(result, tuple) or len(result) != 2:
        raise RuntimeError(f"scrub_doublets returned unexpected result: {result}")

    doublet_score, predicted_doublets = result

    # Snapshot scrublet's auto-detected threshold BEFORE any call_doublets()
    # below may overwrite scrub.threshold_. None means auto-detection failed.
    auto_threshold = getattr(scrub, "threshold_", None)

    if manual_threshold is not None:
        predicted_doublets = scrub.call_doublets(threshold=manual_threshold)

    # Fallback A: scrublet failed to auto-detect a threshold (returned None
    # for predicted_doublets and no manual threshold was provided). Force the
    # hard-coded threshold so the pipeline doesn't crash.
    none_threshold_fallback_applied = False
    if predicted_doublets is None and manual_threshold is None:
        none_threshold_fallback_applied = True
        _log_fallback(
            f"type=none_threshold sample={sample_id} "
            f"fallback_threshold={SCRUBLET_THRESHOLD_FALLBACK} "
            f"reason='scrublet returned None for predicted_doublets (auto-threshold detection failed)'"
        )
        predicted_doublets = scrub.call_doublets(threshold=SCRUBLET_THRESHOLD_FALLBACK)

    if predicted_doublets is None:
        raise RuntimeError(
            f"call_doublets returned None even after fallback to "
            f"threshold={SCRUBLET_THRESHOLD_FALLBACK}. "
            f"Check 'DEBUG_{sample_id}_doublet_histogram_initial.png'."
        )

    adata_QC_min.obs['predicted_doublet'] = predicted_doublets
    adata_QC1 = adata_QC_min[~predicted_doublets, :]

    adata_QC2 = adata_QC1[adata_QC1.obs.pct_counts_mito < max_mito, :]


    # Save metrics for MultiQC
    # adata: STAR-filtered; adata_min: min200/500; adata_QC1: rm doublet; adata_QC2: rm high MT
    number_cells = adata.n_obs
    median_genes = int(np.median(adata.obs['n_genes_by_counts']))
    median_counts = int(np.median(adata.obs['total_counts']))

    number_cells_QC_min = adata_QC_min.n_obs
    median_genes_QC_min = int(np.median(adata_QC_min.obs['n_genes_by_counts']))
    median_counts_QC_min = int(np.median(adata_QC_min.obs['total_counts']))

    number_cells_QC1 = adata_QC1.n_obs

    number_cells_QC2 = adata_QC2.n_obs
    median_genes_QC2 = int(np.median(adata_QC2.obs['n_genes_by_counts']))
    median_counts_QC2 = int(np.median(adata_QC2.obs['total_counts']))

    doublet_fraction = 0
    if number_cells > 0:
        doublet_fraction = round((number_cells_QC_min - number_cells_QC1) / number_cells_QC_min, 4)

######################################## Fallback B: scrublet auto-threshold > 0.4 is too lenient
    threshold_fallback_applied = False
    if (
        manual_threshold is None
        and not none_threshold_fallback_applied
        and auto_threshold is not None
        and auto_threshold < SCRUBLET_THRESHOLD_FALLBACK
    ):
        threshold_fallback_applied = True
        _log_fallback(
            f"type=threshold sample={sample_id} "
            f"auto_threshold={float(auto_threshold):.4f} "
            f"fallback_threshold={SCRUBLET_THRESHOLD_FALLBACK} "
            f"reason='scrublet auto threshold > {SCRUBLET_THRESHOLD_FALLBACK}'"
        )
        predicted_doublets = scrub.call_doublets(threshold=SCRUBLET_THRESHOLD_FALLBACK)

        if predicted_doublets is None:
            raise RuntimeError(
                f"call_doublets returned None for predicted_doublets after "
                f"threshold={SCRUBLET_THRESHOLD_FALLBACK} fallback. "
                f"Check 'DEBUG_{sample_id}_doublet_histogram_initial.png'."
            )

        adata_QC_min.obs['predicted_doublet'] = predicted_doublets
        adata_QC1 = adata_QC_min[~predicted_doublets, :]
        adata_QC2 = adata_QC1[adata_QC1.obs.pct_counts_mito < max_mito, :]

        number_cells_QC1 = adata_QC1.n_obs
        number_cells_QC2 = adata_QC2.n_obs
        median_genes_QC2 = int(np.median(adata_QC2.obs['n_genes_by_counts']))
        median_counts_QC2 = int(np.median(adata_QC2.obs['total_counts']))

        doublet_fraction = 0
        if number_cells > 0:
            doublet_fraction = round((number_cells_QC_min - number_cells_QC1) / number_cells_QC_min, 4)

###################################### Archive 1228
#    raw_matrix = adata.X.todense()
#    expected_doublet_rate = _calculate_expected_doublet_rate(expected_num_cells)
#    scrub = scr.Scrublet(raw_matrix, expected_doublet_rate=expected_doublet_rate)
#    doublet_score, predicted_doublets = scrub.scrub_doublets(
#        min_counts=2,
#        min_cells=3,
#        min_gene_variability_pctl=85,
#        n_prin_comps=30,
#    )
#    adata.obs['predicted_doublet'] = predicted_doublets
#    adata_QC1 = adata[~predicted_doublets, :]

#   adata_QC_min, min_gene_threshold = _apply_mingene_filter(adata_QC1, min_genes)
#   adata_QC2 = adata_QC_min[adata_QC_min.obs.pct_counts_mito < max_mito, :]
######################################
    histogram_plot = scrub.plot_histogram()
    histogram_figure = None
    if hasattr(histogram_plot, "get_figure"):
        histogram_figure = histogram_plot.get_figure()
    elif hasattr(histogram_plot, "savefig"):
        histogram_figure = histogram_plot
    elif isinstance(histogram_plot, (list, tuple)):
        for plot_obj in histogram_plot:
            if hasattr(plot_obj, "get_figure"):
                histogram_figure = plot_obj.get_figure()
                break
            if hasattr(plot_obj, "savefig"):
                histogram_figure = plot_obj
                break
    if histogram_figure is None:
        histogram_figure = plt.gcf()

    # When any threshold fallback fired, annotate each histogram panel:
    #   black solid  = original auto-detected threshold (what scrublet picked)
    #                  -- omitted if auto-detection failed (None case)
    #   red dashed   = fallback threshold actually used (0.4)
    # Red is reserved for the fallback value across all QC plots.
    any_threshold_fallback = threshold_fallback_applied or none_threshold_fallback_applied
    if any_threshold_fallback:
        for ax in histogram_figure.axes:
            if auto_threshold is not None:
                ax.axvline(
                    auto_threshold,
                    color="black",
                    linestyle="-",
                    linewidth=1.2,
                    label=f"auto threshold = {auto_threshold:.2f}",
                )
            ax.axvline(
                SCRUBLET_THRESHOLD_FALLBACK,
                color="red",
                linestyle="--",
                linewidth=1.2,
                label=f"fallback threshold = {SCRUBLET_THRESHOLD_FALLBACK:.2f}",
            )
            ax.legend(fontsize=8, loc="upper right")

    histogram_figure.savefig(
        f"4.{sample_id}_doublet_score_histogram_QC1_mqc.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(histogram_figure)

    # Fig 0: Doublet scatter
    adata_QC_min.obs['predicted_doublet'] = adata_QC_min.obs['predicted_doublet'].astype('category')
    custom_palette = ['#DDDDDD', 'red']
    fig0, ax = plt.subplots()
    sc.pl.scatter(
        adata_QC_min,
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
    sc.pl.violin(adata_QC_min, 'n_genes_by_counts', jitter=0.4, ax=axes[0, 0], show=False)
    axes[0, 1].set_title('Counts per Cell (Before)')
    sc.pl.violin(adata_QC_min, 'total_counts', jitter=0.4, ax=axes[0, 1], show=False)

    axes[1, 0].set_title('Genes per Cell (After Doublet Removal)')
    sc.pl.violin(adata_QC1, 'n_genes_by_counts', jitter=0.4, ax=axes[1, 0], show=False)
    axes[1, 1].set_title('Counts per Cell (After Doublet Removal)')
    sc.pl.violin(adata_QC1, 'total_counts', jitter=0.4, ax=axes[1, 1], show=False)

    fig1.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig1.savefig(f'2.{sample_id}_violin_comparison_QC1_mqc.png')
    plt.close(fig1)

    # Fig 2: Mitochondria removal (vertical orientation for side-by-side panel placement)
    fig2, axes = plt.subplots(2, 1, figsize=(7, 12))
    sc.pl.violin(adata_QC1, 'pct_counts_mito', jitter=0.4, ax=axes[0], show=False)
    axes[0].set_title('pct_counts_mito (Before MT filtering)')
    sc.pl.violin(adata_QC2, 'pct_counts_mito', jitter=0.4, ax=axes[1], show=False)
    axes[1].set_title('pct_counts_mito (After MT filtering)')
    fig2.suptitle(f'{sample_id} - Mito Filtering - QC2 (min_genes={min_gene_threshold})')
    fig2.tight_layout(rect=[0, 0.03, 1, 0.97])
    fig2.savefig(f'3.{sample_id}_violin_mito_filtering_QC2_mqc.png')
    plt.close(fig2)

    # Save whitelist
    adata_QC2.obs_names.to_series().to_csv(f'{sample_id}_whitelist.txt', index=False, header=False)

    with open(f"{sample_id}_total_metaqc_partial.tsv", "w") as f:
        f.write("# plot_type: 'table'\n")
        f.write("# section_name: 'Total_metaQC'\n")
        f.write("# description: 'Combined QC metrics before ambient RNA correction'\n")
        f.write("# pconfig:\n")
        f.write("#     sortRows: false\n")
        f.write(
            "Sample\tInitial_n_cells\tInitial_median_genes\tInitial_median_counts\tQC_min_n_cells\tQC_min_median_genes\tQC_min_median_counts\tDoublet_fraction\tQC_min_db_mt_n_cells\tQC_min_db_mt_median_genes\tQC_min_db_mt_median_counts\tSequencing_saturation\tRho\n"
        )
        f.write(
            f"{sample_id}\t{number_cells}\t{median_genes}\t{median_counts}\t{number_cells_QC_min}\t{median_genes_QC_min}\t{median_counts_QC_min}\t{doublet_fraction}\t{number_cells_QC2}\t{median_genes_QC2}\t{median_counts_QC2}\t{sequencing_saturation}\t\n"
        )

    adata_QC_min.write_h5ad(f"{sample_id}_rmdoublet_mingene{min_gene_threshold}.h5ad")
    adata_QC1.write_h5ad(f"{sample_id}_rmdoublet.h5ad")
    adata_QC2.write_h5ad(f"{sample_id}_rmMT{max_mito}.h5ad")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run Scrublet doublet detection")
    parser.add_argument('--sample_id', required=True, help="Sample identifier")
    parser.add_argument('--matrix_dir', required=True, help="Path to filtered count matrix directory")
    parser.add_argument('--knee_matrix_dir', required=False, help="Path to raw count matrix directory for knee plot")
    parser.add_argument('--min_genes', type=int, default=500, help="Min genes per cell before relaxation")
    parser.add_argument('--min_cells', type=int, default=3, help="Min cells per gene")
    parser.add_argument('--max_mito', type=float, default=0.2, help="Max mitochondrial percentage (fraction)")
    parser.add_argument('--mito_prefixes', nargs='+', default=['mt-'], help="Prefix(es) for mitochondrial genes")
    parser.add_argument('--mito_gene_list', type=str, default=None, help="Path to text file listing mitochondrial genes (one per line)")
    parser.add_argument('--summary_csv', type=str, default=None, help="Path to STARsolo Summary.csv for sequencing saturation")
    parser.add_argument('--manual_threshold', type=float, default=None, help="Optional manual Scrublet threshold override")
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
        args.knee_matrix_dir,
        args.manual_threshold,
    )

#!/usr/bin/env python

import argparse
import time
from typing import Dict, List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy as sp
from numba import njit, prange
from samalg import SAM
from tqdm import tqdm


# ============================================================
# Helpers for marker detection
# ============================================================
def enforce_numeric_cluster_order(adata, cluster_key="leiden_clusters"):
    """
    Ensure cluster labels are strings and ordered categoricals in numeric order,
    so Scanpy plots show 0,1,2,3,... rather than 0,1,10,11,...
    """
    adata.obs[cluster_key] = adata.obs[cluster_key].astype(str)

    cluster_order = sorted(
        adata.obs[cluster_key].unique().tolist(),
        key=lambda x: int(x)  # numeric sort
    )

    adata.obs[cluster_key] = pd.Categorical(
        adata.obs[cluster_key],
        categories=cluster_order,
        ordered=True
    )
    return cluster_order


def get_X_labels_genes(adata, cluster_key="leiden_clusters"):
    X = adata.X
    labels = np.asarray(adata.obs[cluster_key].values)  # will be strings
    genes = np.asarray(adata.var_names)
    return X, labels, genes


def get_gene_vector(Xmat, gene_idx):
    """Return 1D ndarray for gene_idx across all cells; handles sparse/dense."""
    col = Xmat[:, gene_idx]
    if hasattr(col, "toarray"):  # sparse column
        return col.toarray().ravel()
    return np.asarray(col).ravel()


def preselect_topK_candidates_scanpy_wilcoxon(
    adata,
    *,
    cluster_key="leiden_clusters",
    topK=200,
    min_cells_in_cluster=10,
    use_raw=False,
    layer=None,
    corr_method="benjamini-hochberg",
):
    # valid clusters only
    vc = adata.obs[cluster_key].value_counts()
    valid_clusters = vc[vc >= min_cells_in_cluster].index.tolist()

    # run once for all groups
    sc.tl.rank_genes_groups(
        adata,
        groupby=cluster_key,
        method="wilcoxon",
        use_raw=use_raw,
        layer=layer,
        corr_method=corr_method,
        n_genes=topK,
        pts=False,
    )

    rg = adata.uns["rank_genes_groups"]
    names = rg["names"]  # structured array: names[group][rank]

    var_index = {g: i for i, g in enumerate(adata.var_names)}

    candidates = {}
    for cl in valid_clusters:
        gene_names = list(names[cl])
        idx = [var_index[g] for g in gene_names if g in var_index]
        candidates[cl] = np.asarray(idx, dtype=int)

    return candidates


@njit(parallel=True)
def bsMeanDiff(A, B, nBS=500, alpha=0.1):
    diffs = np.zeros(nBS)
    np.random.seed(0)
    for i in prange(nBS):
        a = np.random.choice(A, len(A))
        b = np.random.choice(B, len(B))
        diffs[i] = a.mean() - b.mean()
    return (
        np.quantile(diffs, alpha / 2),
        diffs.mean(),
        np.quantile(diffs, 1 - alpha / 2),
    )


def refine_candidates_for_one_cluster(
    X, labels, genes, cl, gene_idx,
    *,
    nBS=500,
    alpha=0.1,
    expr_thresh=0.25,
    min_pct_in=0.4,
    max_pct_out=0.2,
    min_ratio=4.0,
    max_p=1e-3,
    require_ci_positive=True,
):
    in_mask = (labels == cl)
    out_mask = ~in_mask

    n_in = int(in_mask.sum())
    n_out = int(out_mask.sum())
    if n_in == 0 or n_out == 0:
        return pd.DataFrame()

    gene_idx = np.asarray(gene_idx, dtype=int)
    m = gene_idx.size
    if m == 0:
        return pd.DataFrame()

    ci_lo = np.zeros(m, dtype=float)
    bs_mu = np.zeros(m, dtype=float)
    ci_hi = np.zeros(m, dtype=float)
    p = np.zeros(m, dtype=float)

    mean_in = np.zeros(m, dtype=float)
    mean_out = np.zeros(m, dtype=float)
    pct_in = np.zeros(m, dtype=float)
    pct_out = np.zeros(m, dtype=float)
    ratio = np.zeros(m, dtype=float)

    for k, j in enumerate(gene_idx):
        x = get_gene_vector(X, j)
        x_in = x[in_mask]
        x_out = x[out_mask]

        mean_in[k] = x_in.mean()
        mean_out[k] = x_out.mean()

        pct_in[k] = np.mean(x_in > expr_thresh)
        pct_out[k] = np.mean(x_out > expr_thresh)

        if mean_out[k] > 0:
            ratio[k] = mean_in[k] / mean_out[k]
        else:
            ratio[k] = np.inf if mean_in[k] > 0 else np.nan

        ci_lo[k], bs_mu[k], ci_hi[k] = bsMeanDiff(x_in, x_out, nBS=nBS, alpha=alpha)

        try:
            p[k] = sp.stats.mannwhitneyu(x_in, x_out, alternative="two-sided").pvalue
        except ValueError:
            p[k] = 1.0

    res = pd.DataFrame({
        "cluster": cl,
        "gene": genes[gene_idx],
        "gene_idx": gene_idx,
        "n_in": n_in,
        "n_out": n_out,
        "mean_in": mean_in,
        "mean_out": mean_out,
        "mean_diff_in_minus_out": mean_in - mean_out,
        "ratio_mean_in_over_out": ratio,
        "pct_in": pct_in,
        "pct_out": pct_out,
        "bs_ci_lo": ci_lo,
        "bs_mean_diff": bs_mu,
        "bs_ci_hi": ci_hi,
        "p_value": p,
        "-log10_p": -np.log10(np.maximum(p, 1e-300)),
    })

    pass_mask = (
        (res["pct_in"] >= min_pct_in) &
        (res["pct_out"] <= max_pct_out) &
        (res["ratio_mean_in_over_out"] >= min_ratio) &
        (res["p_value"] <= max_p)
    )
    if require_ci_positive:
        pass_mask &= (res["bs_ci_lo"] > 0)

    res["passes_filter"] = pass_mask.values

    # rank by bs_mean_diff (fallback ordering)
    res = res.sort_values("bs_mean_diff", ascending=False)
    return res


def run_qc2(
    sample_id,
    input_h5ad,
    whitelist,
    n_hvg,
    n_neighbors,
    n_pcs,
    marker_genes_file=None,
    *,
    cluster_key="leiden_clusters",
    topK=200,
    min_cells_in_cluster=10,
    expr_thresh=0.25,
    min_pct_in=0.4,
    max_pct_out=0.2,
    min_ratio=4.0,
    max_p=1e-3,
    require_ci_positive=True,
    nBS=500,
    alpha=0.1,
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
    ## 20251225 added
    sam.preprocess_data()
    sam.run(
        projection='umap',
        preprocessing='StandardScaler',
        weight_mode = 'rms',
        k=n_neighbors,
        npcs=n_pcs,
        n_genes=n_hvg,
    )
    ###
    sam.leiden_clustering(res=0.5)

    cluster_order = enforce_numeric_cluster_order(sam.adata, cluster_key=cluster_key)
    cluster_labels = sam.adata.obs[cluster_key].astype(str)
    clusters_sorted = cluster_order

    cluster_stats = (
        sam.adata.obs.assign(**{cluster_key: cluster_labels})
        .groupby(cluster_key)[['n_counts', 'n_genes']]
        .median()
        .loc[clusters_sorted]
    )

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw={'width_ratios': [2, 1]})
    
    sc.pl.umap(sam.adata,color=[cluster_key],legend_loc='on data',s = 40, ax=axes[0],show=False, title=f'{sample_id} - Leiden Clustering (UMAP) - QC2')

    y_pos = np.arange(len(clusters_sorted))
    bar_height = 0.35
    y_counts = y_pos - bar_height / 2
    y_genes = y_pos + bar_height / 2
    ax_counts = axes[1]
    ax_genes = ax_counts.twiny()

    counts_bars = ax_counts.barh(
        y_counts,
        cluster_stats['n_counts'],
        height=bar_height,
        color='steelblue',
        label='n_counts',
        align='center'
    )
    genes_bars = ax_genes.barh(
        y_genes,
        cluster_stats['n_genes'],
        height=bar_height,
        color='orange',
        alpha=0.7,
        label='n_genes',
        align='center'
    )

    ax_counts.set_yticks(y_pos)
    ax_counts.set_yticklabels(clusters_sorted)
    ax_counts.invert_yaxis()
    ax_counts.set_xlabel('Median n_counts')
    ax_counts.set_title('Median n_counts and n_genes by Leiden cluster')

    ax_genes.set_xlabel('Median n_genes')
    ax_genes.grid(False)
    ax_counts.legend(handles=[counts_bars, genes_bars], loc='lower right')

    ### Added 12/23 waiting to be tested
    ax_counts.xaxis.set_ticks_position('top')
    ax_counts.xaxis.set_label_position('top')
    ax_genes.xaxis.set_ticks_position('bottom')
    ax_genes.xaxis.set_label_position('bottom')

    plt.tight_layout()
    plt.savefig(
        f'4.{sample_id}_umap_leiden_QC2_mqc.png', dpi=300, bbox_inches='tight'
    )
    plt.close(fig)

    # ============================================================
    # Marker detection pipeline: Scanpy Wilcoxon pre-screen + refined filtering
    # ============================================================
    cluster_order = enforce_numeric_cluster_order(sam.adata, cluster_key=cluster_key)

    X, labels, genes = get_X_labels_genes(sam.adata, cluster_key=cluster_key)
    clusters = np.asarray(cluster_order, dtype=object)

    t_total_start = time.perf_counter()

    # ---- prescreen timing ----
    t_prescreen_start = time.perf_counter()

    candidates = preselect_topK_candidates_scanpy_wilcoxon(
        sam.adata,
        cluster_key=cluster_key,
        topK=topK,
        min_cells_in_cluster=min_cells_in_cluster,
        use_raw=False,
        layer=None,
    )

    t_prescreen_end = time.perf_counter()
    print(f"[TIMER] Prescreen (Scanpy Wilcoxon): {t_prescreen_end - t_prescreen_start:.2f} s")

    # ---- refine timing ----
    t_refine_start = time.perf_counter()

    all_ranked_candidates: Dict[str, pd.DataFrame] = {}
    all_strict_markers: Dict[str, pd.DataFrame] = {}
    top5_rows = []
    cluster_times: Dict[str, float] = {}

    for cl in tqdm(clusters, desc="refine clusters", unit="cluster"):
        if cl not in candidates or candidates[cl].size == 0:
            continue

        t0 = time.perf_counter()

        res_cand = refine_candidates_for_one_cluster(
            X, labels, genes, cl, candidates[cl],
            nBS=nBS,
            alpha=alpha,
            expr_thresh=expr_thresh,
            min_pct_in=min_pct_in,
            max_pct_out=max_pct_out,
            min_ratio=min_ratio,
            max_p=max_p,
            require_ci_positive=require_ci_positive,
        )

        t1 = time.perf_counter()
        cluster_times[cl] = t1 - t0

        all_ranked_candidates[cl] = res_cand
        all_strict_markers[cl] = res_cand.loc[res_cand["passes_filter"]].copy()

        strict_ranked = (
            res_cand.loc[res_cand["passes_filter"]]
            .sort_values(
                ["pct_out", "ratio_mean_in_over_out", "p_value", "bs_mean_diff"],
                ascending=[True, False, True, False],
            )
        )

        top5 = strict_ranked.head(5).copy()

        if len(top5) < 5:
            need = 5 - len(top5)
            filler = (
                res_cand.sort_values("bs_mean_diff", ascending=False)
                .loc[~res_cand["gene"].isin(top5["gene"])]
                .head(need)
            )
            top5 = pd.concat([top5, filler], ignore_index=True)

        top5["rank_in_cluster"] = np.arange(1, len(top5) + 1)
        top5["passed_strict_filter"] = top5["passes_filter"].values
        top5_rows.append(top5)

    t_refine_end = time.perf_counter()
    print(f"[TIMER] Refinement (bootstrap on candidates): {t_refine_end - t_refine_start:.2f} s")

    if len(top5_rows) > 0:
        top5_markers_by_cluster = pd.concat(top5_rows, ignore_index=True)

        top5_markers_by_cluster["_cluster_num"] = top5_markers_by_cluster["cluster"].astype(int)
        top5_markers_by_cluster = (
            top5_markers_by_cluster
            .sort_values(["_cluster_num", "rank_in_cluster"])
            .drop(columns="_cluster_num")
        )

        top5_markers_by_cluster = top5_markers_by_cluster[
            [
                "cluster", "rank_in_cluster", "gene",
                "bs_mean_diff", "bs_ci_lo", "bs_ci_hi",
                "p_value", "-log10_p",
                "pct_in", "pct_out",
                "mean_in", "mean_out", "ratio_mean_in_over_out",
                "passes_filter", "passed_strict_filter",
            ]
        ]
    else:
        top5_markers_by_cluster = pd.DataFrame()

    t_total_end = time.perf_counter()

    # ---- timing summary ----
    print("\n========= TIMING SUMMARY =========")
    print(f"Total runtime: {t_total_end - t_total_start:.2f} s")
    print(f"Prescreen runtime: {t_prescreen_end - t_prescreen_start:.2f} s")
    print(f"Refinement runtime: {t_refine_end - t_refine_start:.2f} s")

    if len(cluster_times) > 0:
        ct_times = np.array(list(cluster_times.values()))
        print(f"Avg refine time / cluster: {ct_times.mean():.2f} s")
        print(f"Min / Max refine time: {ct_times.min():.2f} / {ct_times.max():.2f} s")
    print("==================================\n")

    if len(top5_rows) > 0:
        top5_markers_by_cluster.to_csv(f"{sample_id}_markers.csv", index=False)

    sam.adata.uns['marker_genes'] = {
        'top5_markers_by_cluster': top5_markers_by_cluster.to_dict(orient='list') if len(top5_rows) > 0 else {},
        'all_ranked_candidates': {
            cl: df.reset_index(drop=True).to_dict(orient='list') for cl, df in all_ranked_candidates.items()
        },
        'strict_markers': {
            cl: df.reset_index(drop=True).to_dict(orient='list') for cl, df in all_strict_markers.items()
        },
        'timing': {
            'total_seconds': t_total_end - t_total_start,
            'prescreen_seconds': t_prescreen_end - t_prescreen_start,
            'refine_seconds': t_refine_end - t_refine_start,
            'per_cluster_seconds': cluster_times,
        },
    }

    genes_for_dotplot: List[str] = []
    if len(top5_rows) > 0:
        for cluster in cluster_order:
            cluster_top = top5_markers_by_cluster[top5_markers_by_cluster["cluster"] == cluster]
            for gene in cluster_top["gene"]:
                if gene not in genes_for_dotplot:
                    genes_for_dotplot.append(gene)

    if genes_for_dotplot:
        sc.pl.dotplot(
            sam.adata,
            genes_for_dotplot,
            cluster_key,
            standard_scale='var',
            dendrogram=False,
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
    sam.adata.obs[cluster_key].to_csv(f'{sample_id}_{cluster_key}.csv', header=True)
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
    parser.add_argument('--cluster_key', type=str, default='leiden_clusters', help="Obs column to use for clustering")
    parser.add_argument('--topK', type=int, default=200, help="Number of top genes per cluster from Wilcoxon prescreen")
    parser.add_argument('--min_cells_in_cluster', type=int, default=10, help="Minimum cells required for a cluster to be processed")
    parser.add_argument('--expr_thresh', type=float, default=0.25, help="Expression threshold for pct_in/pct_out calculations")
    parser.add_argument('--min_pct_in', type=float, default=0.4, help="Minimum pct_in for refined marker filtering")
    parser.add_argument('--max_pct_out', type=float, default=0.2, help="Maximum pct_out for refined marker filtering")
    parser.add_argument('--min_ratio', type=float, default=4.0, help="Minimum mean-in/mean-out ratio")
    parser.add_argument('--max_p', type=float, default=1e-3, help="Maximum p-value for refined marker filtering")
    parser.add_argument('--require_ci_positive', action='store_true', default=True, help="Require bootstrap CI to be strictly positive")
    parser.add_argument('--no_require_ci_positive', action='store_false', dest='require_ci_positive', help="Disable bootstrap CI positivity requirement")
    parser.add_argument('--nBS', type=int, default=500, help="Number of bootstrap iterations")
    parser.add_argument('--alpha', type=float, default=0.1, help="Alpha level for bootstrap confidence intervals")
    args = parser.parse_args()

    run_qc2(
        args.sample_id,
        args.input_h5ad,
        args.whitelist,
        args.n_hvg,
        args.n_neighbors,
        args.n_pcs,
        args.marker_genes,
        cluster_key=args.cluster_key,
        topK=args.topK,
        min_cells_in_cluster=args.min_cells_in_cluster,
        expr_thresh=args.expr_thresh,
        min_pct_in=args.min_pct_in,
        max_pct_out=args.max_pct_out,
        min_ratio=args.min_ratio,
        max_p=args.max_p,
        require_ci_positive=args.require_ci_positive,
        nBS=args.nBS,
        alpha=args.alpha,
    )

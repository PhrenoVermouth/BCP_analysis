#!/usr/bin/env python

import argparse
import time
from typing import Dict, List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sparse
from samalg import SAM


# ============================================================
# Helpers for marker detection (vectorized expression-based filtering)
# ============================================================
def enforce_numeric_cluster_order(adata, cluster_key="leiden_clusters"):
    """
    Ensure cluster labels are strings and ordered categoricals in numeric order,
    so Scanpy plots show 0,1,2,3,... rather than 0,1,10,11,...
    """
    labels_str = adata.obs[cluster_key].astype(str).values
    try:
        cats_sorted = [str(x) for x in np.sort(pd.to_numeric(pd.unique(labels_str)))]
    except Exception:
        cats_sorted = sorted(pd.unique(labels_str))

    adata.obs[cluster_key] = pd.Categorical(
        adata.obs[cluster_key].astype(str),
        categories=cats_sorted,
        ordered=True,
    )
    return cats_sorted


def get_cluster_stats(adata, cluster_key: str, expr_thresh: float):
    """
    Pre-compute mean and fraction matrices for all clusters (core speedup step).

    Returns:
        means_df: DataFrame (index=clusters, columns=genes) - mean expression per cluster
        fracs_df: DataFrame (index=clusters, columns=genes) - fraction expressed per cluster
        sizes: Series (index=clusters) - number of cells per cluster
        global_mean: ndarray - global mean expression per gene
        global_frac: ndarray - global fraction expressed per gene
    """
    # Use CSR format for fast row slicing (cells)
    X = adata.X
    if sparse.issparse(X) and not sparse.isspmatrix_csr(X):
        X = X.tocsr()

    # Binary matrix for fraction calculation
    if sparse.issparse(X):
        X_bin = X > expr_thresh
    else:
        X_bin = X > expr_thresh

    clusters = adata.obs[cluster_key].astype("category")
    unique_clusters = clusters.cat.categories

    # Initialize result containers
    means_df = pd.DataFrame(
        index=unique_clusters, columns=adata.var_names, dtype=np.float32
    )
    fracs_df = pd.DataFrame(
        index=unique_clusters, columns=adata.var_names, dtype=np.float32
    )
    sizes = clusters.value_counts()

    # Iterate over clusters (only ~20 iterations, much faster than per-gene)
    for ct in unique_clusters:
        mask = (clusters == ct).values
        if not np.any(mask):
            continue

        X_ct = X[mask]
        X_bin_ct = X_bin[mask]

        # Mean (axis=0 for column mean)
        m = X_ct.mean(axis=0)
        means_df.loc[ct] = np.asarray(m).flatten().astype(np.float32)

        # Fraction expressed
        f = X_bin_ct.mean(axis=0)
        fracs_df.loc[ct] = np.asarray(f).flatten().astype(np.float32)

    # Global statistics
    global_mean = np.asarray(X.mean(axis=0)).flatten()
    global_frac = np.asarray((X > expr_thresh).mean(axis=0)).flatten()

    return means_df, fracs_df, sizes, global_mean, global_frac


def vectorized_marker_scoring(
    adata,
    cluster_key: str,
    means_df: pd.DataFrame,
    fracs_df: pd.DataFrame,
    sizes: pd.Series,
    global_mean: np.ndarray,
    global_frac: np.ndarray,
    *,
    expr_thresh: float,
    min_frac_in: float,
    max_clusters: int,
    diff_thresh: float,
    frac_gain_thresh: float,
    test_top: int,
    top_n: int,
    penalize_bg_expr: bool,
) -> tuple:
    """
    Vectorized marker scoring logic - loops over clusters, not genes.
    """
    total_cells = adata.n_obs

    # Compute specificity constraint: how many clusters express each gene
    n_clusters_expressed = (fracs_df >= min_frac_in).sum(axis=0)

    # Get Wilcoxon results
    rank_res = sc.get.rank_genes_groups_df(adata, group=None)

    all_rows: List[pd.DataFrame] = []
    top_genes_by_cluster: Dict[str, List[str]] = {}

    # Loop over clusters (not genes!)
    for ct in means_df.index:
        # Get Wilcoxon candidates for this cluster
        ct_rank_df = rank_res[rank_res["group"] == ct].head(test_top)
        if ct_rank_df.empty:
            continue

        candidate_genes = ct_rank_df["names"].values
        pvals = ct_rank_df["pvals"].values
        pval_map = pd.Series(pvals, index=candidate_genes)

        # Filter to genes that exist in our data
        valid_genes = [g for g in candidate_genes if g in means_df.columns]
        if not valid_genes:
            continue

        # Get gene indices for lookup
        gene_indices = [adata.var_names.get_loc(g) for g in valid_genes]

        # Extract target cluster statistics (vectorized)
        ct_mean = means_df.loc[ct, valid_genes].values.astype(np.float64)
        ct_frac = fracs_df.loc[ct, valid_genes].values.astype(np.float64)

        # Calculate background statistics
        ct_n = sizes[ct]
        bg_n = total_cells - ct_n

        # Global stats for these genes
        g_mean_sub = global_mean[gene_indices]
        g_frac_sub = global_frac[gene_indices]

        # Derive background mean: bg_mean = (total_sum - ct_sum) / bg_n
        ct_sum = ct_mean * ct_n
        total_sum = g_mean_sub * total_cells
        bg_mean = np.maximum((total_sum - ct_sum) / bg_n, 0)

        # Derive background fraction
        ct_count_expr = ct_frac * ct_n
        total_count_expr = g_frac_sub * total_cells
        bg_frac = np.maximum((total_count_expr - ct_count_expr) / bg_n, 0)

        # ---- Apply vectorized filters ----
        # 1. Min fraction constraint
        mask_fct = ct_frac >= min_frac_in

        # 2. Max clusters constraint
        n_cl_vals = n_clusters_expressed.loc[valid_genes].values
        mask_spec = n_cl_vals <= max_clusters

        # 3. Diff & gain thresholds
        mean_diff = ct_mean - bg_mean
        safe_fct = np.maximum(ct_frac, 1e-12)
        frac_gain = (ct_frac - bg_frac) / safe_fct
        mask_diff = (mean_diff > diff_thresh) | (frac_gain > frac_gain_thresh)

        # Combined mask
        final_mask = mask_fct & mask_spec & mask_diff

        if not np.any(final_mask):
            continue

        # ---- Compute scores (vectorized) ----
        sel_idx = np.where(final_mask)[0]
        sel_genes = np.array(valid_genes)[sel_idx]
        sel_mean_ct = ct_mean[sel_idx]
        sel_mean_bg = bg_mean[sel_idx]
        sel_mean_diff = mean_diff[sel_idx]
        sel_frac_gain = frac_gain[sel_idx]
        sel_fct = ct_frac[sel_idx]
        sel_fbg = bg_frac[sel_idx]
        sel_n_cl = n_cl_vals[sel_idx]

        # Get p-values
        sel_p = pval_map.loc[sel_genes].values.astype(np.float64)
        sel_p = np.maximum(sel_p, 1e-300)

        # Expression term
        expr_term = np.log1p(sel_mean_ct)
        if penalize_bg_expr:
            expr_term -= 0.5 * np.log1p(sel_mean_bg)

        # Composite score
        score = (
            np.clip(expr_term, 0, None)
            * np.clip(sel_frac_gain, 0, None)
            * np.clip(sel_mean_diff, 0, None)
            * (-np.log10(sel_p))
        )

        # Build result DataFrame
        df_ct = pd.DataFrame(
            {
                "gene": sel_genes,
                "cluster": ct,
                "p": sel_p,
                "mean_ct": sel_mean_ct,
                "mean_bg": sel_mean_bg,
                "mean_diff": sel_mean_diff,
                "fct": sel_fct,
                "fbg": sel_fbg,
                "frac_gain": sel_frac_gain,
                "n_clusters_expressed": sel_n_cl,
                "score": score,
            }
        )

        # Sort and take top N
        df_ct = df_ct.sort_values(["score", "mean_ct"], ascending=False).head(top_n)
        df_ct["rank_in_cluster"] = np.arange(1, len(df_ct) + 1)
        all_rows.append(df_ct)
        top_genes_by_cluster[ct] = df_ct["gene"].tolist()

    if not all_rows:
        return pd.DataFrame(), {}

    final_df = pd.concat(all_rows, ignore_index=True)
    return final_df, top_genes_by_cluster


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
    # New marker detection parameters
    expr_thresh=0.25,       # threshold for "expressed"
    min_frac_in=0.20,       # gene must be expressed in >= this fraction of cells in target cluster
    max_clusters=5,         # gene can be expressed (>= min_frac_in) in at most this many clusters
    diff_thresh=1.0,        # mean_ct - mean_bg threshold (prefilter)
    frac_gain_thresh=0.7,   # (fct - fbg)/fct threshold (prefilter)
    test_top=5000,          # candidates per cluster to test
    top_n=5,                # top markers per cluster to keep
    penalize_bg_expr=False, # whether to penalize background expression in score
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
    ## 20251225 added, 1230 changed to None
    sam.preprocess_data()
    sam.run(
        projection='umap',
        preprocessing=None,
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
    # Marker detection pipeline: Vectorized expression-based filtering
    # ============================================================
    adata = sam.adata
    cluster_order = enforce_numeric_cluster_order(adata, cluster_key=cluster_key)
    clusters_str = cluster_order

    t_total_start = time.perf_counter()

    # ---- Step 1: Pre-compute cluster statistics (vectorized) ----
    t_stats_start = time.perf_counter()
    print("Computing cluster statistics (vectorized)...",flush=True)
    means_df, fracs_df, sizes, global_mean, global_frac = get_cluster_stats(
        adata, cluster_key, expr_thresh
    )
    t_stats_end = time.perf_counter()
    print(f"[TIMER] Cluster stats: {t_stats_end - t_stats_start:.2f} s",flush=True)

    # ---- Step 2: Scanpy Wilcoxon precandidates ----
    t_prescreen_start = time.perf_counter()
    sc.tl.rank_genes_groups(
        adata,
        groupby=cluster_key,
        method="wilcoxon",
        n_genes=min(test_top, adata.n_vars),
        use_raw=False,
        pts=True,
    )
    t_prescreen_end = time.perf_counter()
    print(f"[TIMER] Wilcoxon prescreen: {t_prescreen_end - t_prescreen_start:.2f} s")

    # ---- Step 3: Vectorized marker scoring ----
    t_refine_start = time.perf_counter()
    print("Scoring markers (vectorized)...",flush=True)
    markers_df_all, top_genes_by_cluster = vectorized_marker_scoring(
        adata,
        cluster_key,
        means_df,
        fracs_df,
        sizes,
        global_mean,
        global_frac,
        expr_thresh=expr_thresh,
        min_frac_in=min_frac_in,
        max_clusters=max_clusters,
        diff_thresh=diff_thresh,
        frac_gain_thresh=frac_gain_thresh,
        test_top=test_top,
        top_n=top_n,
        penalize_bg_expr=penalize_bg_expr,
    )
    t_refine_end = time.perf_counter()
    t_total_end = time.perf_counter()

    # ---- Timing summary ----
    print("\n========= TIMING SUMMARY =========",flush=True)
    print(f"Total runtime: {t_total_end - t_total_start:.2f} s")
    print(f"  - Cluster stats: {t_stats_end - t_stats_start:.2f} s")
    print(f"  - Wilcoxon prescreen: {t_prescreen_end - t_prescreen_start:.2f} s")
    print(f"  - Vectorized scoring: {t_refine_end - t_refine_start:.2f} s")
    print("==================================\n")

    print(f"Clusters total: {len(clusters_str)}")
    print(f"Clusters with markers: {len(top_genes_by_cluster)}")
    print(f"Constraints: min_frac_in={min_frac_in}, max_clusters={max_clusters}, expr_thresh={expr_thresh}")

    # Save results
    if not markers_df_all.empty:
        markers_df_all.to_csv(f"{sample_id}_markers.csv", index=False)

        # Build top markers table
        top_markers_rows = []
        for ct in clusters_str:
            if ct in top_genes_by_cluster:
                ct_df = markers_df_all[markers_df_all["cluster"] == ct]
                for rank_i, gene in enumerate(top_genes_by_cluster[ct], 1):
                    gene_row = ct_df[ct_df["gene"] == gene]
                    if len(gene_row) > 0:
                        row = gene_row.iloc[0].to_dict()
                        row["rank_in_cluster"] = rank_i
                        top_markers_rows.append(row)

        if top_markers_rows:
            top_markers_df = pd.DataFrame(top_markers_rows)
            top_markers_df.to_csv(f"{sample_id}_top_markers.csv", index=False)

    sam.adata.uns['marker_genes'] = {
        'top_genes_by_cluster': top_genes_by_cluster,
        'all_markers': markers_df_all.to_dict(orient='list') if not markers_df_all.empty else {},
        'parameters': {
            'expr_thresh': expr_thresh,
            'min_frac_in': min_frac_in,
            'max_clusters': max_clusters,
            'diff_thresh': diff_thresh,
            'frac_gain_thresh': frac_gain_thresh,
            'test_top': test_top,
            'top_n': top_n,
        },
        'timing': {
            'total_seconds': t_total_end - t_total_start,
            'stats_seconds': t_stats_end - t_stats_start,
            'prescreen_seconds': t_prescreen_end - t_prescreen_start,
            'scoring_seconds': t_refine_end - t_refine_start,
        },
    }

    # ---- Dotplot (grouped by cluster) ----
    var_names_dict = {ct: top_genes_by_cluster[ct] for ct in clusters_str if ct in top_genes_by_cluster}

    if not var_names_dict:
        print("No clusters had any passing markers. Debug:")
        print("  categories:", clusters_str[:10])
        print("  rank_genes_groups keys:", list(adata.uns.get("rank_genes_groups", {}).keys()))
    else:
        sc.pl.dotplot(
            sam.adata,
            var_names=var_names_dict,
            groupby=cluster_key,
            standard_scale="var",
            dot_max=0.8,
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
    parser.add_argument('--n_neighbors', type=int, default=20, help="Number of nearest neighbors in PCA space")
    parser.add_argument('--n_pcs', type=int, default=150, help="Number of principal components")
    parser.add_argument('--marker_genes', type=str, default=None, help="Optional path to marker genes list (one gene per line)")
    parser.add_argument('--cluster_key', type=str, default='leiden_clusters', help="Obs column to use for clustering")
    # New marker detection parameters
    parser.add_argument('--expr_thresh', type=float, default=0.25, help="Expression threshold for considering a gene 'expressed'")
    parser.add_argument('--min_frac_in', type=float, default=0.20, help="Gene must be expressed in >= this fraction of cells in target cluster")
    parser.add_argument('--max_clusters', type=int, default=5, help="Gene can be expressed (>= min_frac_in) in at most this many clusters")
    parser.add_argument('--diff_thresh', type=float, default=1.0, help="mean_ct - mean_bg threshold (prefilter)")
    parser.add_argument('--frac_gain_thresh', type=float, default=0.7, help="(fct - fbg)/fct threshold (prefilter)")
    parser.add_argument('--test_top', type=int, default=5000, help="Candidates per cluster to test from Wilcoxon prescreen")
    parser.add_argument('--top_n', type=int, default=5, help="Top markers per cluster to keep")
    parser.add_argument('--penalize_bg_expr', action='store_true', default=False, help="Penalize background expression in score calculation")
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
        expr_thresh=args.expr_thresh,
        min_frac_in=args.min_frac_in,
        max_clusters=args.max_clusters,
        diff_thresh=args.diff_thresh,
        frac_gain_thresh=args.frac_gain_thresh,
        test_top=args.test_top,
        top_n=args.top_n,
        penalize_bg_expr=args.penalize_bg_expr,
    )

#!/usr/bin/env python
"""
SAM-based marker calling for multiome integration.

Called from within run_archr_multiome.R via system().
Uses the same vectorized marker scoring algorithm as run_sam_qc.py.

Input:  RNA counts h5ad (X = raw counts) + cluster assignment CSV
Output: marker CSV + SAM-processed h5ad (X = SAM-processed, .raw = counts)
"""

import argparse
import time
from typing import Dict, List

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sparse
from samalg import SAM


# ============================================================
# Helpers (same algorithm as run_sam_qc.py)
# ============================================================
def enforce_numeric_cluster_order(adata, cluster_key="cluster"):
    labels_str = adata.obs[cluster_key].astype(str).values
    unique_labels = list(pd.unique(labels_str))
    try:
        # Try pure numeric first (1, 2, 10)
        cats_sorted = [str(x) for x in np.sort(pd.to_numeric(unique_labels))]
    except (ValueError, TypeError):
        # Try natural sort for prefixed labels (C1, C2, C10)
        import re
        def _nat_key(s):
            parts = re.split(r'(\d+)', s)
            return [int(p) if p.isdigit() else p.lower() for p in parts]
        cats_sorted = sorted(unique_labels, key=_nat_key)
    adata.obs[cluster_key] = pd.Categorical(
        adata.obs[cluster_key].astype(str),
        categories=cats_sorted,
        ordered=True,
    )
    return cats_sorted


def get_cluster_stats(adata, cluster_key, expr_thresh):
    X = adata.X
    if sparse.issparse(X) and not sparse.isspmatrix_csr(X):
        X = X.tocsr()
    X_bin = X > expr_thresh

    clusters = adata.obs[cluster_key].astype("category")
    unique_clusters = clusters.cat.categories
    means_df = pd.DataFrame(index=unique_clusters, columns=adata.var_names, dtype=np.float32)
    fracs_df = pd.DataFrame(index=unique_clusters, columns=adata.var_names, dtype=np.float32)
    sizes = clusters.value_counts()

    for ct in unique_clusters:
        mask = (clusters == ct).values
        if not np.any(mask):
            continue
        X_ct = X[mask]
        X_bin_ct = X_bin[mask]
        means_df.loc[ct] = np.asarray(X_ct.mean(axis=0)).flatten().astype(np.float32)
        fracs_df.loc[ct] = np.asarray(X_bin_ct.mean(axis=0)).flatten().astype(np.float32)

    global_mean = np.asarray(X.mean(axis=0)).flatten()
    global_frac = np.asarray((X > expr_thresh).mean(axis=0)).flatten()
    return means_df, fracs_df, sizes, global_mean, global_frac


def vectorized_marker_scoring(
    adata, cluster_key, means_df, fracs_df, sizes, global_mean, global_frac,
    *, expr_thresh, min_frac_in, max_clusters, diff_thresh, frac_gain_thresh,
    test_top, top_n, penalize_bg_expr,
):
    total_cells = adata.n_obs
    n_clusters_expressed = (fracs_df >= min_frac_in).sum(axis=0)
    rank_res = sc.get.rank_genes_groups_df(adata, group=None)

    all_rows: List[pd.DataFrame] = []
    top_genes_by_cluster: Dict[str, List[str]] = {}

    for ct in means_df.index:
        ct_rank_df = rank_res[rank_res["group"] == ct].head(test_top)
        if ct_rank_df.empty:
            continue
        candidate_genes = ct_rank_df["names"].values
        pvals = ct_rank_df["pvals"].values
        pval_map = pd.Series(pvals, index=candidate_genes)
        valid_genes = [g for g in candidate_genes if g in means_df.columns]
        if not valid_genes:
            continue

        gene_indices = [adata.var_names.get_loc(g) for g in valid_genes]
        ct_mean = means_df.loc[ct, valid_genes].values.astype(np.float64)
        ct_frac = fracs_df.loc[ct, valid_genes].values.astype(np.float64)
        ct_n = sizes[ct]
        bg_n = total_cells - ct_n
        g_mean_sub = global_mean[gene_indices]
        g_frac_sub = global_frac[gene_indices]
        ct_sum = ct_mean * ct_n
        total_sum = g_mean_sub * total_cells
        bg_mean = np.maximum((total_sum - ct_sum) / bg_n, 0)
        ct_count_expr = ct_frac * ct_n
        total_count_expr = g_frac_sub * total_cells
        bg_frac = np.maximum((total_count_expr - ct_count_expr) / bg_n, 0)

        mask_fct = ct_frac >= min_frac_in
        n_cl_vals = n_clusters_expressed.loc[valid_genes].values
        mask_spec = n_cl_vals <= max_clusters
        mean_diff = ct_mean - bg_mean
        safe_fct = np.maximum(ct_frac, 1e-12)
        frac_gain = (ct_frac - bg_frac) / safe_fct
        mask_diff = (mean_diff > diff_thresh) | (frac_gain > frac_gain_thresh)
        final_mask = mask_fct & mask_spec & mask_diff

        if not np.any(final_mask):
            continue

        sel_idx = np.where(final_mask)[0]
        sel_genes = np.array(valid_genes)[sel_idx]
        sel_mean_ct = ct_mean[sel_idx]
        sel_mean_bg = bg_mean[sel_idx]
        sel_mean_diff = mean_diff[sel_idx]
        sel_frac_gain = frac_gain[sel_idx]
        sel_fct = ct_frac[sel_idx]
        sel_fbg = bg_frac[sel_idx]
        sel_n_cl = n_cl_vals[sel_idx]
        sel_p = pval_map.loc[sel_genes].values.astype(np.float64)
        sel_p = np.maximum(sel_p, 1e-300)

        expr_term = np.log1p(sel_mean_ct)
        if penalize_bg_expr:
            expr_term -= 0.5 * np.log1p(sel_mean_bg)

        score = (
            np.clip(expr_term, 0, None)
            * np.clip(sel_frac_gain, 0, None)
            * np.clip(sel_mean_diff, 0, None)
            * (-np.log10(sel_p))
        )

        df_ct = pd.DataFrame({
            "gene": sel_genes, "cluster": ct, "p": sel_p,
            "mean_ct": sel_mean_ct, "mean_bg": sel_mean_bg, "mean_diff": sel_mean_diff,
            "fct": sel_fct, "fbg": sel_fbg, "frac_gain": sel_frac_gain,
            "n_clusters_expressed": sel_n_cl, "score": score,
        })
        df_ct = df_ct.sort_values(["score", "mean_ct"], ascending=False).head(top_n)
        df_ct["rank_in_cluster"] = np.arange(1, len(df_ct) + 1)
        all_rows.append(df_ct)
        top_genes_by_cluster[ct] = df_ct["gene"].tolist()

    if not all_rows:
        return pd.DataFrame(), {}
    return pd.concat(all_rows, ignore_index=True), top_genes_by_cluster


# ============================================================
# Main
# ============================================================
def main():
    parser = argparse.ArgumentParser(description="SAM marker calling for multiome")
    parser.add_argument("--input", required=True, help="RNA counts h5ad (X = raw counts)")
    parser.add_argument("--clusters", required=True, help="CSV with barcode,cluster columns")
    parser.add_argument("--sample_id", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--n_hvg", type=int, default=3000)
    parser.add_argument("--n_neighbors", type=int, default=20)
    parser.add_argument("--n_pcs", type=int, default=150)
    parser.add_argument("--top_n", type=int, default=5)
    parser.add_argument("--expr_thresh", type=float, default=0.25)
    parser.add_argument("--min_frac_in", type=float, default=0.20)
    parser.add_argument("--max_clusters", type=int, default=2)
    parser.add_argument("--diff_thresh", type=float, default=1.0)
    parser.add_argument("--frac_gain_thresh", type=float, default=0.7)
    parser.add_argument("--test_top", type=int, default=5000)
    args = parser.parse_args()

    import os
    os.makedirs(args.outdir, exist_ok=True)

    # 1. Load h5ad (X = raw counts)
    print(f"[SAM] Loading RNA counts from {args.input}", flush=True)
    adata = sc.read_h5ad(args.input)
    adata.var_names_make_unique()

    # 2. Load cluster assignments from ArchR
    print("[SAM] Loading cluster assignments...", flush=True)
    cl_df = pd.read_csv(args.clusters, index_col=0)
    # Align to adata barcodes
    common = adata.obs_names.intersection(cl_df.index)
    adata = adata[common, :].copy()
    adata.obs["cluster"] = cl_df.loc[common, "cluster"].values
    cluster_key = "cluster"

    # 3. SAM preprocessing (saves counts to .raw before transforming)
    print("[SAM] Running SAM preprocessing...", flush=True)
    sc.pp.filter_cells(adata, min_genes=0)
    adata.obs["n_counts"] = np.asarray(adata.X.sum(axis=1)).flatten()
    adata.raw = adata.copy()  # preserve raw counts in .raw

    sam = SAM(adata)
    sam.preprocess_data()
    sam.run(
        projection="umap",
        weight_mode="rms",
        k=args.n_neighbors,
        npcs=args.n_pcs,
        n_genes=args.n_hvg,
    )

    # Restore cluster labels (SAM may reorder obs)
    adata = sam.adata
    if "cluster" not in adata.obs.columns:
        adata.obs["cluster"] = cl_df.loc[adata.obs_names, "cluster"].values

    # 4. Marker calling with provided clusters (no re-clustering)
    print("[SAM] Running marker calling with provided cluster labels...", flush=True)
    cluster_order = enforce_numeric_cluster_order(adata, cluster_key=cluster_key)

    t_start = time.perf_counter()

    means_df, fracs_df, sizes, global_mean, global_frac = get_cluster_stats(
        adata, cluster_key, args.expr_thresh
    )

    sc.tl.rank_genes_groups(
        adata, groupby=cluster_key, method="wilcoxon",
        n_genes=min(args.test_top, adata.n_vars), use_raw=False, pts=True,
    )

    markers_df, top_genes_by_cluster = vectorized_marker_scoring(
        adata, cluster_key, means_df, fracs_df, sizes, global_mean, global_frac,
        expr_thresh=args.expr_thresh, min_frac_in=args.min_frac_in,
        max_clusters=args.max_clusters, diff_thresh=args.diff_thresh,
        frac_gain_thresh=args.frac_gain_thresh, test_top=args.test_top,
        top_n=args.top_n, penalize_bg_expr=False,
    )

    t_end = time.perf_counter()
    print(f"[SAM] Marker calling completed in {t_end - t_start:.2f}s", flush=True)
    print(f"[SAM] Clusters with markers: {len(top_genes_by_cluster)}/{len(cluster_order)}", flush=True)

    # 5. Save outputs
    marker_csv = os.path.join(args.outdir, f"{args.sample_id}_sam_markers.csv")
    if not markers_df.empty:
        markers_df.to_csv(marker_csv, index=False)
        print(f"[SAM] Marker table saved: {marker_csv}", flush=True)
    else:
        # Write empty CSV with header so R can still read it
        pd.DataFrame(columns=["gene", "cluster", "score", "rank_in_cluster"]).to_csv(
            marker_csv, index=False
        )
        print("[SAM] WARNING: No markers found, wrote empty table", flush=True)

    # Save SAM-processed h5ad (.raw = counts, X = SAM-processed)
    h5ad_out = os.path.join(args.outdir, f"{args.sample_id}_sam_processed.h5ad")
    adata.write_h5ad(h5ad_out)
    print(f"[SAM] Processed h5ad saved: {h5ad_out}", flush=True)

    print("[SAM] Done.", flush=True)


if __name__ == "__main__":
    main()

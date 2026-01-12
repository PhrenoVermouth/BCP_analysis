#!/usr/bin/env python

"""
Attach RNA velocity layers (spliced/unspliced) to an AnnData object without
altering the count matrix used for existing analyses.
"""

import argparse
import anndata as ad
import pandas as pd
import scanpy as sc
import scipy.sparse as sp


def load_velocity_layers(velo_dir: str):
    """Load spliced/unspliced matrices and their metadata from a Velocyto dir."""
    spliced = sc.read_mtx(f"{velo_dir}/spliced.mtx.gz").T
    unspliced = sc.read_mtx(f"{velo_dir}/unspliced.mtx.gz").T

    barcodes = pd.read_csv(f"{velo_dir}/barcodes.tsv.gz", header=None, sep="\t")[0].values
    genes = pd.read_csv(f"{velo_dir}/features.tsv.gz", header=None, sep="\t")[0].values

    velo_adata = ad.AnnData(
        X=spliced.X,
        obs=pd.DataFrame(index=barcodes),
        var=pd.DataFrame(index=genes),
        layers={
            "spliced": spliced.X,
            "unspliced": unspliced.X,
        },
    )
    velo_adata.var_names_make_unique()
    velo_adata.obs_names_make_unique()
    return velo_adata


def add_velocity_layers(base_h5ad: str, velo_dir: str, output: str):
    base = sc.read_h5ad(base_h5ad)
    base.var_names_make_unique()
    base.obs_names_make_unique()

    velo = load_velocity_layers(velo_dir)

    common_cells = base.obs_names.intersection(velo.obs_names)
    common_genes = base.var_names.intersection(velo.var_names)

    if len(common_cells) == 0 or len(common_genes) == 0:
        raise ValueError("No overlapping cells or genes found between counts and velocity matrices.")

    base_subset = base[common_cells, common_genes].copy()
    velo_subset = velo[common_cells, common_genes].copy()

    base_subset.layers["spliced"] = velo_subset.layers["spliced"].copy()
    base_subset.layers["unspliced"] = velo_subset.layers["unspliced"].copy()

    base_subset.write_h5ad(output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Attach Velocyto layers to an AnnData file.")
    parser.add_argument("--counts_h5ad", required=True, help="Path to the counts-based h5ad (unchanged).")
    parser.add_argument("--velo_dir", required=True, help="Path to the Velocyto matrix directory (raw or filtered).")
    parser.add_argument("--output", required=True, help="Output h5ad path with velocity layers attached.")
    args = parser.parse_args()

    add_velocity_layers(args.counts_h5ad, args.velo_dir, args.output)

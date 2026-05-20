#!/bin/bash
# Script to download 10x Genomics E18 Mouse Brain Multiome Test Dataset

mkdir -p e18_mouse_brain_test
cd e18_mouse_brain_test

echo "Downloading GEX (RNA) FASTQs..."
curl -O https://cf.10xgenomics.com/samples/cell-arc/1.0.0/e18_mouse_brain_fresh_5k/e18_mouse_brain_fresh_5k_gex_fastqs.tar

echo "Downloading ATAC FASTQs..."
curl -O https://cf.10xgenomics.com/samples/cell-arc/1.0.0/e18_mouse_brain_fresh_5k/e18_mouse_brain_fresh_5k_atac_fastqs.tar

echo "Unpacking FASTQs..."
tar -xvf e18_mouse_brain_fresh_5k_gex_fastqs.tar
tar -xvf e18_mouse_brain_fresh_5k_atac_fastqs.tar

# Download mm10 reference (if you don't already have one configured for ATAC and STAR)
# echo "Downloading mm10 ARC reference (approx 20GB)..."
# curl -O https://cf.10xgenomics.com/supp/cell-arc/refdata-cellranger-arc-mm10-2020-A-2.0.0.tar.gz
# tar -xzvf refdata-cellranger-arc-mm10-2020-A-2.0.0.tar.gz

echo "Test files downloaded to: $(pwd)"

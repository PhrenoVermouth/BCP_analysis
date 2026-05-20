#!/bin/bash
TAR="/home/gyang/SHERLOCK_software/BCP_analysis/new_resources_multiome/cellranger-arc-2.1.0.tar.gz"
OUTDIR="/home/gyang/SHERLOCK_software/BCP_analysis/new_resources_multiome"

echo "Extracting ATAC whitelist..."
tar -xzf $TAR cellranger-arc-2.1.0/lib/python/atac/barcodes/737K-arc-v1.txt.gz -O > $OUTDIR/multiome_atac.txt.gz

echo "Extracting RNA whitelist..."
tar -xzf $TAR cellranger-arc-2.1.0/lib/python/cellranger/barcodes/737K-arc-v1.txt.gz -O > $OUTDIR/multiome_rna.txt.gz

echo "Unzipping..."
gunzip -f $OUTDIR/multiome_atac.txt.gz
gunzip -f $OUTDIR/multiome_rna.txt.gz

echo "Generating translation CSV..."
paste -d ',' $OUTDIR/multiome_atac.txt $OUTDIR/multiome_rna.txt > $OUTDIR/barcode_translation.csv

echo "Done."

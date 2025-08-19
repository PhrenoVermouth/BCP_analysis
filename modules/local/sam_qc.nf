// modules/local/sam_qc.nf

process SAM_QC {
    tag "$meta.id"
    publishDir "$params.outdir/sam_qc/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(corrected_h5ad)

    output:
    tuple val(meta), path("*.png"), emit: qc_plots
    tuple val(meta), path("*.h5ad"), emit: filtered_adata
    path "*_cells.tsv", emit: qc_cells_metrics
    path "*_counts.tsv", emit: qc_counts_metrics
    path "*_genes.tsv", emit: qc_genes_metrics

    script:
    def mito_prefixes = params.mito_prefixes.join(' ')
    def marker_arg = params.marker_genes ? "\\\n        --marker_genes ${params.marker_genes}" : ""
    """
    run_sam_qc.py \\
        --sample_id ${meta.id} \\
        --input_h5ad ${corrected_h5ad} \\
        --min_genes ${params.min_genes_per_cell} \\
        --min_cells ${params.min_cells_per_gene} \\
        --max_genes ${params.max_genes_per_cell} \\
        --max_counts ${params.max_counts_per_cell} \\
        --max_mito ${params.max_mito} \\
        --mito_prefixes ${mito_prefixes} \\
        --n_hvg ${params.n_hvg} \\
        --n_neighbors ${params.n_neighbors} \\
        --n_pcs ${params.n_pcs}  ${marker_arg}
    """
}

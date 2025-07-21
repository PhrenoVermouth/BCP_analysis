// modules/local/scanpy_qc.nf

process SCANPY_QC {
    tag "$meta.id"
    publishDir "$params.outdir/scanpy_qc/${meta.id}", mode: 'copy'
    
    input:
    tuple val(meta), path(corrected_h5ad)

    output:
    tuple val(meta), path("*.png"), emit: qc_plots
    tuple val(meta), path("*.h5ad"), emit: filtered_adata
    path "*_qc_metrics.tsv", emit: qc_metrics

    script:
    def mito_prefixes = params.mito_prefixes.join(' ')
    """
    run_scanpy_qc.py \\
        --sample_id ${meta.id} \\
        --input_h5ad ${corrected_h5ad} \\
        --min_genes ${params.min_genes_per_cell} \\
        --min_cells ${params.min_cells_per_gene} \\
        --max_genes ${params.max_genes_per_cell} \\
        --max_counts ${params.max_counts_per_cell} \\
        --max_mito ${params.max_mito_percent} \\
        --mito_prefixes ${mito_prefixes}
    """
}

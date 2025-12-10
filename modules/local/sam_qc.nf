// modules/local/sam_qc.nf

process SAM_QC {
    tag "$meta.id"
    publishDir "$params.outdir/sam_qc/${meta.id}", mode: 'copy'
    errorStrategy 'ignore'
    input:
    tuple val(meta), path(corrected_h5ad), path(whitelist)

    output:
    tuple val(meta), path("*.png"), emit: qc_plots
    tuple val(meta), path("*.h5ad"), emit: filtered_adata

    script:
    def marker_arg = params.marker_genes ? "\\\n        --marker_genes ${params.marker_genes}" : ""
    """
    run_sam_qc.py \
        --sample_id ${meta.id} \
        --input_h5ad ${corrected_h5ad} \
        --whitelist ${whitelist} \
        --n_hvg ${params.n_hvg} \
        --n_neighbors ${params.n_neighbors} \
        --n_pcs ${params.n_pcs}  ${marker_arg}
    """
}

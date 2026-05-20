// modules/local/qc_panel.nf

process QC_PANEL {
    tag "$meta.id"
    publishDir "$params.outdir/qc_panel/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(knee_plot), path(histogram), path(violin_qc1), path(mito_qc2), path(soupx_combined), path(umap), path(dotplot)

    output:
    tuple val(meta), path("${meta.id}_qc_panel_mqc.png"), emit: panel

    script:
    """
    assemble_qc_panel.R \
        --sample_id ${meta.id} \
        --knee_plot ${knee_plot} \
        --histogram ${histogram} \
        --violin_qc1 ${violin_qc1} \
        --mito_qc2 ${mito_qc2} \
        --soupx_combined ${soupx_combined} \
        --umap ${umap} \
        --dotplot ${dotplot} \
        --output ${meta.id}_qc_panel_mqc.png
    """
}

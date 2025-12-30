// modules/local/scrublet.nf

process SCRUBLET {
    tag "$meta.id"
    publishDir "$params.outdir/scrublet/${meta.id}", mode: 'copy'
    errorStrategy 'ignore'
    input:
    tuple val(meta), path(gzipped_dir)

    output:
    tuple val(meta), path("${meta.id}_whitelist.txt"), emit: whitelist
    tuple val(meta), path("*.png"), emit: qc_plots
    tuple val(meta), path("*_total_metaqc_partial.tsv"), emit: metaqc_partial
    tuple val(meta), path("*.h5ad"), emit: scrublet_h5ads

    script:
    def mito_prefixes_str = params.mito_prefixes instanceof List ? params.mito_prefixes.join(' ') : params.mito_prefixes
    def mito_list_arg = params.mito_gene_list ? "\\\n        --mito_gene_list ${params.mito_gene_list}" : ""
    def mito_max = meta.containsKey('max_mito') ? meta.max_mito : params.max_mito
    def summary_csv = "${gzipped_dir}/GeneFull/Summary.csv"
    def knee_matrix_dir = "${gzipped_dir}/GeneFull/raw"
    """
    run_scrublet.py \
        --sample_id ${meta.id} \
        --matrix_dir ${gzipped_dir}/GeneFull/filtered \
        --knee_matrix_dir ${knee_matrix_dir} \
        --min_genes ${params.min_genes_per_cell} \
        --min_cells ${params.min_cells_per_gene} \
        --max_mito ${mito_max} \
        --summary_csv ${summary_csv} \
        --mito_prefixes ${mito_prefixes_str} ${mito_list_arg}
    """
}

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
    path "*_cells_mqc.tsv", emit: qc_cells_metrics
    path "*_counts_mqc.tsv", emit: qc_counts_metrics
    path "*_genes_mqc.tsv", emit: qc_genes_metrics

    script:
    def mito_prefixes = params.mito_prefixes.join(' ')
    def mito_list_arg = params.mito_gene_list ? "\\\n        --mito_gene_list ${params.mito_gene_list}" : ""
    def summary_csv = "${gzipped_dir}/GeneFull/Summary.csv"
    """
    run_scrublet.py \
        --sample_id ${meta.id} \
        --matrix_dir ${gzipped_dir}/GeneFull/raw \
        --min_genes ${params.min_genes_per_cell} \
        --min_cells ${params.min_cells_per_gene} \
        --max_mito ${params.max_mito} \
        --summary_csv ${summary_csv} \
        --mito_prefixes ${mito_prefixes} ${mito_list_arg}
    """
}

// modules/local/soupx.nf

process SOUPX {
    tag "$meta.id"
    errorStrategy 'ignore'
    publishDir "$params.outdir/soupx/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(gzipped_dir), path(whitelist)
    output:
    tuple val(meta), path("${meta.id}_pre_soupx.h5ad"), emit: pre_h5ad
    tuple val(meta), path("${meta.id}_rm_ambient.h5ad"), emit: ambient_h5ad
    tuple val(meta), path("${meta.id}_soupx_rho.tsv"), emit: rho
    path("0.${meta.id}_ambient_RNA_removed_mqc.png"), emit: ambient_plot
    path("0.${meta.id}_soupx_contamination_estimation_mqc.png"), emit: contamination_plot
    tuple val(meta), path("0.${meta.id}_soupx_combined_mqc.png"), emit: combined_plot

    script:
    """
    Rscript ${baseDir}/bin/run_soupx.R \
        --raw_dir ${gzipped_dir}/GeneFull/raw \
        --filtered_dir ${gzipped_dir}/GeneFull/filtered \
        --whitelist ${whitelist} \
        --sample_id ${meta.id}

    """
}

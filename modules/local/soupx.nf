// modules/local/soupx.nf

process SOUPX {
    tag "$meta.id"
    errorStrategy 'ignore'
    publishDir "$params.outdir/soupx/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(gzipped_dir)
    output:
    tuple val(meta), path("${meta.id}_corrected.h5ad"), emit: corrected_h5ad
    tuple val(meta), path("${meta.id}_pre_soupx.h5ad"), emit: pre_h5ad
    path("0.${meta.id}_ambient_RNA_removed.png"), emit: ambient_plot
    path("0.${meta.id}_soupx_contamination_estimation.png"), emit: contamination_plot
    script:
    """
    Rscript ${baseDir}/bin/run_soupx.R \
        --raw_dir ${gzipped_dir}/GeneFull/raw \
        --filtered_dir ${gzipped_dir}/GeneFull/filtered \
        --sample_id ${meta.id}

    """
}

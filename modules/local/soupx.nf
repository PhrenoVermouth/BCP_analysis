// modules/local/soupx.nf

process SOUPX {
    tag "$meta.id"
    publishDir "$params.outdir/soupx/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(gzipped_solo_dir) 
    output:
    tuple val(meta), path("${meta.id}_corrected.h5ad"), emit: corrected_h5ad

    script:
    """
    Rscript ${baseDir}/bin/run_soupx.R \\
        --raw_dir ${gzipped_solo_dir}/GeneFull/raw \\
        --filtered_dir ${gzipped_solo_dir}/GeneFull/filtered \\
        --sample_id ${meta.id}
    """
}

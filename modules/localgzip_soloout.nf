// modules/local/gzip_soloout.nf

process GZIP_SOLO_OUTPUT {
    tag "$meta.id"
    publishDir "$params.outdir/gzipped_matrix/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(solo_out_dir) // 输入来自STAR_SOLO的顶层输出目录

    output:
    tuple val(meta), path(solo_out_dir), emit: gzipped_dir 

    script:
    """
    find ${solo_out_dir}/GeneFull -type f \\( -name "*.tsv" -o -name "*.mtx" \\) -print0 | xargs -0 -r gzip
    """
}

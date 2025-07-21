// modules/local/gzip_solout.nf

process GZIP_SOLO_OUTPUT {
    tag "$meta.id"
    publishDir "$params.outdir/gzipped_matrix/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(solo_out_dir) 

    output:
    tuple val(meta), path(solo_out_dir), emit: gzipped_dir 

    script:
    """
    gzip ${solo_out_dir}/GeneFull/raw/barcodes.tsv
    gzip ${solo_out_dir}/GeneFull/raw/features.tsv
    gzip ${solo_out_dir}/GeneFull/raw/matrix.mtx
    """
}

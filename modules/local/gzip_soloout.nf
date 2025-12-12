// modules/local/gzip_soloout.nf

process GZIP_SOLO_OUTPUT {
    tag "$meta.id"
    publishDir "$params.outdir/gzipped_matrix/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(solo_out_dir) 

    output:
    tuple val(meta), path(solo_out_dir), emit: gzipped_dir

    script:
    def feature_dirs = []
    if (params.run_mode.toLowerCase() == 'genefull') {
        feature_dirs << "${solo_out_dir}/GeneFull"
    } else {
        feature_dirs << "${solo_out_dir}/Velocyto"
    }
    def dirs_to_zip = feature_dirs.collect { "\"${it}\"" }.join(' ')
    """
    for dir in ${dirs_to_zip}; do
        if [[ -d \${dir} ]]; then
            find \${dir} -type f \\( -name "*.tsv" -o -name "*.mtx" \\) -print0 | xargs -0 -r gzip
        fi
    done
    """
}

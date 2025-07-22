// modules/local/multiqc.nf

process MULTIQC {
    tag "MultiQC Report"
    publishDir "$params.outdir/multiqc", mode: 'copy'

    input:
    path files
    path config

    output:
    path "multiqc_report.html"
    path "multiqc_data"

    script:
    """
    multiqc . -c ${config} -o . -n multiqc_report.html
    """
}

// modules/local/multiqc.nf

process MULTIQC {
    tag "MultiQC Report"
    publishDir "$params.outdir/multiqc", mode: 'copy'

    input:
    path files

    output:
    path "multiqc_report.html"

    script:
    """
    multiqc .
    """
}

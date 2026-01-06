// modules/local/multiqc.nf

process MULTIQC {
    tag "MultiQC Report"
    publishDir "$params.outdir/multiqc", mode: 'copy'

    input:
    path files collect: true
    path config

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data", emit: data

    script:
    """
    multiqc . -f  -c ${config} -o . -n multiqc_report.html
    mv multiqc_report_data multiqc_data
    """
}

// modules/local/multiqc.nf

def reportName = params.run_mode?.toLowerCase() == 'velocity'
    ? 'multiqc_report_velocyto.html'
    : 'multiqc_report.html'

process MULTIQC {
    tag "MultiQC Report"
    publishDir "$params.outdir/multiqc", mode: 'copy'

    input:
    path files
    path config

    output:
    path reportName, emit: report
    path "multiqc_data", emit: data

    script:
    """
    multiqc . -f  -c ${config} -o . -n ${reportName}
    mv multiqc_report_data multiqc_data
    """
}

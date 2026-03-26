// modules/local/debug_analysis.nf

process DEBUG_ANALYSIS {
    tag "Debug and AI failure analysis"
    publishDir "$params.outdir/pipeline_info", mode: 'copy'

    input:
    path multiqc_report
    val trace_file
    val work_dir
    val base_url
    val api_key

    output:
    path "debug_payload.json"
    path "debug_report.md"
    path "debug_report.html"
    path "ai_analysis.md"
    path "ai_analysis.html"

    script:
    def apiKeyExport = api_key ? """export OPENAI_API_KEY='${api_key.replace("'", "'\"'\"'")}'""" : "true"
    """
    set -euo pipefail

    TRACE_PATH="${trace_file}"
    WORK_DIR="${work_dir}"
    BASE_URL="${base_url}"

    if [ ! -f "\$TRACE_PATH" ]; then
      echo "Trace file not found: \$TRACE_PATH" >&2
      exit 1
    fi

    ${apiKeyExport}

    if [ -z \"\${OPENAI_API_KEY:-}\" ]; then
      echo \"OPENAI_API_KEY is not set. Pass --debug_api_key or export OPENAI_API_KEY before running.\" >&2
      exit 1
    fi

    python ${projectDir}/bin/nf_debug_agent.py --trace "\$TRACE_PATH" --work-dir "\$WORK_DIR" --report debug_report.md --html debug_report.html > debug_payload.json

    python ${projectDir}/bin/analyze_failures.py debug_payload.json --html ai_analysis.html --base-url "\$BASE_URL" > ai_analysis.md
    """
}

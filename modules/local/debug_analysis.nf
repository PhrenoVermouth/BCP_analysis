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
    path "debug_report_*.md"
    path "debug_report_*.html"
    path "ai_analysis_*.md"
    path "ai_analysis_*.html"

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

    timestamp="$(date +%Y%m%d_%H%M%S)"
    debug_report_md="debug_report_${timestamp}.md"
    debug_report_html="debug_report_${timestamp}.html"
    ai_analysis_md="ai_analysis_${timestamp}.md"
    ai_analysis_html="ai_analysis_${timestamp}.html"

    python ${projectDir}/bin/nf_debug_agent.py --trace "\$TRACE_PATH" --work-dir "\$WORK_DIR" --report "\$debug_report_md" --html "\$debug_report_html" > debug_payload.json

    python ${projectDir}/bin/analyze_failures.py debug_payload.json --html "\$ai_analysis_html" --base-url "\$BASE_URL" > "\$ai_analysis_md"
    """
}

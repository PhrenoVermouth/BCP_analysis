from __future__ import annotations
import sys
import argparse
import csv
import json
import re
import math
import html
from pathlib import Path
from typing import Dict, List, Optional, Sequence

# ==========================================
# Part 1: Argument Parsing & Setup
# ==========================================

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Collect forensic information for failed Nextflow tasks."
    )
    parser.add_argument(
        "--trace",
        required=True,
        help="Path to the Nextflow trace execution report file (txt/tsv).",
    )
    parser.add_argument(
        "--work-dir",
        required=True,
        help="Root path to the Nextflow work directory.",
    )
    parser.add_argument(
        "--report",
        default="debug_report.md",
        help="Output path for Markdown report (default: debug_report.md).",
    )
    parser.add_argument(
        "--html",
        default="debug_report.html",
        help="Output path for HTML report (default: debug_report.html).",
    )
    return parser.parse_args()


def sniff_delimiter(trace_path: Path) -> str:
    try:
        with trace_path.open("r", encoding="utf-8", errors="replace") as handle:
            sample = handle.read(4096)
            dialect = csv.Sniffer().sniff(sample, delimiters=["\t", ",", ";"])
            return dialect.delimiter
    except Exception:
        return "\t"


def read_trace_rows(trace_path: Path) -> List[Dict[str, str]]:
    delimiter = sniff_delimiter(trace_path)
    rows: List[Dict[str, str]] = []
    with trace_path.open("r", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        for row in reader:
            rows.append(row)
    return rows


def build_task_dir(work_dir: Path, task_hash: str) -> Path:
    clean_hash = task_hash.replace("/", "")
    if len(clean_hash) < 2:
        return work_dir / clean_hash

    bucket = clean_hash[:2]
    remainder = clean_hash[2:]
    bucket_dir = work_dir / bucket

    if not bucket_dir.exists():
        return bucket_dir / remainder

    try:
        candidates = list(bucket_dir.glob(f"{remainder}*"))
        if candidates:
            return candidates[0]
    except Exception:
        pass

    try:
        candidates_full = list(bucket_dir.glob(f"{clean_hash}*"))
        if candidates_full:
            return candidates_full[0]
    except Exception:
        pass

    return bucket_dir / remainder


def tail_lines(lines: Sequence[str], count: int) -> str:
    return "".join(lines[-count:]) if lines else ""

def format_size_human(size_bytes: int) -> str:
    if size_bytes == 0:
        return "0 B"
    size_name = ("B", "KB", "MB", "GB", "TB")
    i = int(math.floor(math.log(size_bytes, 1024)))
    p = math.pow(1024, i)
    s = round(size_bytes / p, 2)
    return f"{s} {size_name[i]}"

# ==========================================
# Part 2: Forensics & Collection
# ==========================================

def read_command_script(task_dir: Path) -> Optional[str]:
    command_path = task_dir / ".command.sh"
    if not command_path.exists():
        return None
    try:
        return command_path.read_text(encoding="utf-8", errors="replace")
    except Exception:
        return None


def read_log_content(task_dir: Path) -> Dict[str, Optional[str]]:
    err_path = task_dir / ".command.err"
    log_path = task_dir / ".command.log"

    if err_path.exists():
        try:
            size = err_path.stat().st_size
        except Exception:
            size = 0
        if size >= 50:
            try:
                lines = err_path.read_text(encoding="utf-8", errors="replace").splitlines(keepends=True)
                if len(lines) > 200:
                    content = tail_lines(lines, 200)
                else:
                    content = "".join(lines)
                return {"source": ".command.err", "content": content}
            except Exception:
                pass

    if log_path.exists():
        try:
            lines = log_path.read_text(encoding="utf-8", errors="replace").splitlines(keepends=True)
            limited_lines = lines[-100:] if len(lines) > 100 else lines
            content = "".join(limited_lines)
            return {"source": ".command.log", "content": content}
        except Exception:
            pass

    if err_path.exists():
        try:
            return {"source": ".command.err", "content": err_path.read_text(encoding="utf-8", errors="replace")}
        except Exception:
            pass

    return {"source": None, "content": None}


# Marker emitted by bin/run_scrublet.py when a fallback rule activates.
# Format example:
#   [FALLBACK] type=cell_count sample=Sample1 original_cells=51234 keep_top=30000 reason='...'
#   [FALLBACK] type=threshold sample=Sample1 auto_threshold=0.5123 fallback_threshold=0.4 reason='...'
FALLBACK_LINE_RE = re.compile(r"\[FALLBACK\]\s*(.+)")
FALLBACK_KV_RE = re.compile(r"(\w+)=('[^']*'|\"[^\"]*\"|\S+)")


def _parse_fallback_details(details: str) -> Dict[str, str]:
    parsed: Dict[str, str] = {}
    for key, value in FALLBACK_KV_RE.findall(details):
        if (value.startswith("'") and value.endswith("'")) or (
            value.startswith('"') and value.endswith('"')
        ):
            value = value[1:-1]
        parsed[key] = value
    return parsed


def scan_for_fallback_events(task_dir: Path) -> List[Dict[str, str]]:
    """Read .command.log / .command.err in full and pull out [FALLBACK] lines.

    Returns a list of dicts with at least: details (raw), and any key=value
    fields parsed from the marker (e.g. type, sample, original_cells, ...).
    """
    if not task_dir.exists():
        return []

    events: List[Dict[str, str]] = []
    seen_raw: set = set()
    for fname in (".command.log", ".command.err", ".command.out"):
        path = task_dir / fname
        if not path.exists():
            continue
        try:
            content = path.read_text(encoding="utf-8", errors="replace")
        except Exception:
            continue
        for match in FALLBACK_LINE_RE.finditer(content):
            details = match.group(1).strip()
            if details in seen_raw:
                continue
            seen_raw.add(details)
            event = _parse_fallback_details(details)
            event["details"] = details
            event["source"] = fname
            events.append(event)
    return events


def read_solo_summary(task_dir: Path) -> Dict[str, str]:
    summaries = {}
    for path in task_dir.iterdir():
        if path.name.endswith(".Solo.out") and path.is_dir():
            summary_csv = path / "GeneFull" / "Summary.csv"
            if not summary_csv.exists():
                summary_csv = path / "Summary.csv"
            if summary_csv.exists():
                try:
                    content = summary_csv.read_text(encoding="utf-8", errors="replace")
                    summaries[path.name] = content
                except Exception:
                    pass
    return summaries


def scan_input_files(task_dir: Path) -> List[Dict[str, object]]:
    file_list: List[Dict[str, object]] = []
    found_directories: List[str] = []

    for path in task_dir.iterdir():
        if path.name.startswith(".command") or path.name.startswith(".exit"):
            continue

        if path.is_dir():
            found_directories.append(path.name)
            continue

        try:
            stat_result = path.stat()
            size = stat_result.st_size
            file_list.append({
                "path": str(path.name),
                "size_bytes": size,
                "size_human": format_size_human(size)
            })
        except Exception as exc:
            file_list.append({"path": str(path.name), "message": f"unable to stat: {exc}"})
            continue

    if "_STARtmp" in found_directories:
        file_list.append({
            "path": "_STARtmp",
            "size_bytes": 0,
            "size_human": "DIR",
            "message": "⚠️ CRASH INDICATOR: Temporary directory exists"
        })

    return file_list


def detect_related_scripts(
    command_content: Optional[str],
    task_dir: Path,
    script_library: Dict[str, str]
) -> List[str]:
    if not command_content:
        return []

    pattern = re.compile(
        r"(?:python(?:3)?|Rscript|bash|sh|perl).*?\s+([a-zA-Z0-9_./-]+\.(?:py|R|r|sh|pl|bash))",
        re.IGNORECASE,
    )

    referenced_scripts: List[str] = []

    for match in pattern.finditer(command_content):
        script_path_str = match.group(1)
        path_obj = Path(script_path_str)
        candidates = []

        if path_obj.is_absolute():
            candidates.append(path_obj)

        candidates.append(task_dir / path_obj.name)

        if not path_obj.is_absolute():
            candidates.append(task_dir / path_obj)

        found_path = None
        for p in candidates:
            if p.exists() and p.is_file():
                found_path = p
                break

        if not found_path:
            continue

        script_name = found_path.name
        if script_name not in script_library:
            try:
                content = found_path.read_text(encoding="utf-8", errors="replace")
                script_library[script_name] = content
            except Exception:
                script_library[script_name] = f"[Error reading script file at {found_path}]"

        if script_name not in referenced_scripts:
            referenced_scripts.append(script_name)

    return referenced_scripts


def collect_task_report(
    row: Dict[str, str],
    work_dir: Path,
    script_library: Dict[str, str]
) -> Dict[str, object]:
    task_hash = row.get("hash") or row.get("task_hash") or ""
    task_dir = build_task_dir(work_dir, task_hash)

    report: Dict[str, object] = {
        "process_name": row.get("process") or row.get("name"),
        "task_id": row.get("task_id"),
        "hash": task_hash,
        "status": row.get("status"),
        "exit_code": row.get("exit") or row.get("exit_code"),
        "work_dir": str(task_dir),
    }

    try:
        if not task_dir.exists():
            report["error"] = f"Directory not found: {task_dir}"
            return report

        command_script = read_command_script(task_dir)
        report["command_script"] = command_script
        report["logs"] = read_log_content(task_dir)
        report["file_list"] = scan_input_files(task_dir)
        report["related_scripts"] = detect_related_scripts(
            command_script, task_dir, script_library
        )
        report["solo_summary"] = read_solo_summary(task_dir)
        report["fallback_events"] = scan_for_fallback_events(task_dir)

    except Exception as exc:
        report["error"] = f"failed to parse task directory {task_dir}: {exc}"

    return report

# ==========================================
# Part 3: Report Generation (Markdown & HTML)
# ==========================================

def _fallback_table_columns(events: List[Dict[str, object]]) -> List[str]:
    """Pick a stable set of columns for the fallback events table."""
    preferred = [
        "type", "sample", "process", "status",
        "original_cells", "keep_top", "min_umi_kept",
        "auto_threshold", "fallback_threshold",
        "reason", "hash",
    ]
    present = {k for ev in events for k in ev.keys()}
    cols = [c for c in preferred if c in present]
    # Append any leftover keys (sorted) we did not anticipate
    extras = sorted(present - set(cols) - {"details", "source", "task_id", "work_dir"})
    cols.extend(extras)
    return cols


def generate_markdown_report(payload: Dict[str, object], output_path: str):
    failures = payload.get("failures", [])
    fallback_events = payload.get("fallback_events", [])
    script_library = payload.get("script_library", {})

    md_lines = []
    md_lines.append(f"# 🧬 Nextflow Failure Analysis Report")
    md_lines.append(
        f"**Total Failures:** {len(failures)}   "
        f"**Fallback Events:** {len(fallback_events)}\n"
    )

    if fallback_events:
        md_lines.append("## 🟡 Fallback Events")
        md_lines.append(
            "Automatic fallbacks triggered inside tasks (emitted via `[FALLBACK]` log lines)."
        )
        cols = _fallback_table_columns(fallback_events)
        header = "| " + " | ".join(cols) + " |"
        sep = "| " + " | ".join(["---"] * len(cols)) + " |"
        md_lines.append(header)
        md_lines.append(sep)
        for ev in fallback_events:
            cells = []
            for c in cols:
                val = ev.get(c, "")
                if c == "hash" and isinstance(val, str) and val:
                    val = f"`{val[:8]}…`"
                cells.append(str(val))
            md_lines.append("| " + " | ".join(cells) + " |")
        md_lines.append("")

    md_lines.append("## 🔍 Failed Tasks Summary")

    for i, fail in enumerate(failures, 1):
        process = fail.get('process_name', 'Unknown')
        exit_code = fail.get('exit_code', 'N/A')
        task_hash = fail.get('hash', '')
        icon = "❌"
        if str(exit_code) == "137": icon = "💀 (OOM)"

        md_lines.append(f"### {i}. {process} (Exit: {exit_code}) {icon}")
        md_lines.append(f"- **Task Hash:** `{task_hash}`")
        md_lines.append(f"- **Work Dir:** `{fail.get('work_dir')}`")

        if fail.get('error'):
             md_lines.append(f"\n> **⚠️ Parser Error:** {fail.get('error')}")

        task_fallbacks = fail.get("fallback_events", [])
        if task_fallbacks:
            md_lines.append("\n**🟡 Fallback Events for this task:**")
            for ev in task_fallbacks:
                md_lines.append(f"- `[FALLBACK] {ev.get('details', '')}`")

        solo_sums = fail.get("solo_summary", {})
        if solo_sums:
             md_lines.append("\n**📊 STARsolo Metrics (Summary.csv):**")
             for name, content in solo_sums.items():
                 md_lines.append(f"- *Source: `{name}`*")
                 md_lines.append("```csv")
                 md_lines.append(content.strip())
                 md_lines.append("```")

        files = fail.get("file_list", [])
        if files:
            md_lines.append("\n**📂 Task Directory Files:**")
            md_lines.append("| File Name | Size | Notes |")
            md_lines.append("| --- | --- | --- |")
            for f in files:
                name = f.get('path')
                size_str = f.get('size_human', 'N/A')
                msg = f.get('message', '')
                if f.get('size_bytes') == 0 and not msg:
                    msg = "⚠️ Empty"
                md_lines.append(f"| `{name}` | {size_str} | {msg} |")

        scripts = fail.get("related_scripts", [])
        if scripts:
            md_lines.append(f"\n**📜 Script(s) Involved:** `{', '.join(scripts)}`")

        logs = fail.get("logs", {})
        if logs and logs.get("content"):
            source = logs.get("source", "Log")
            md_lines.append(f"\n<details><summary><b>📄 View Error Log ({source})</b></summary>\n")
            md_lines.append("```text")
            md_lines.append(logs['content'])
            md_lines.append("```\n</details>\n")

        md_lines.append("---\n")

    if script_library:
        md_lines.append("## 📚 Source Code Library (Reference)")
        for name, content in script_library.items():
            md_lines.append(f"<details><summary><b>📝 {name} (Click to expand)</b></summary>\n")
            lang = "bash"
            if name.endswith(".R") or name.endswith(".r"): lang = "r"
            if name.endswith(".py"): lang = "python"
            md_lines.append(f"``` {lang}")
            md_lines.append(content)
            md_lines.append("```\n</details>\n")

    with open(output_path, "w", encoding="utf-8") as f:
        f.write("\n".join(md_lines))
    print(f"Human-readable MD report saved to: {output_path}", file=sys.stderr)


def generate_html_report(payload: Dict[str, object], output_path: str):
    failures = payload.get("failures", [])
    fallback_events = payload.get("fallback_events", [])
    script_library = payload.get("script_library", {})

    style = """
    <style>
        body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif; padding: 20px; max-width: 1200px; margin: auto; background-color: #f5f7fa; color: #333; }
        .header { margin-bottom: 30px; border-bottom: 2px solid #eaeaea; padding-bottom: 20px; text-align: center; }
        .card { background: white; border: 1px solid #e1e4e8; margin-bottom: 20px; padding: 25px; border-radius: 8px; box-shadow: 0 4px 6px rgba(0,0,0,0.05); }
        .card.oom { border-left: 6px solid #24292e; }
        .card.error { border-left: 6px solid #d73a49; }
        h1, h2, h3 { margin-top: 0; color: #24292e; }
        table { border-collapse: collapse; width: 100%; margin-top: 15px; font-size: 0.9em; }
        th, td { border: 1px solid #e1e4e8; padding: 10px; text-align: left; }
        th { background-color: #f6f8fa; font-weight: 600; }
        tr:hover { background-color: #f1f8ff; }
        .badge { display: inline-block; padding: 4px 8px; border-radius: 2em; font-size: 0.85em; font-weight: bold; border: 1px solid transparent; }
        .badge-red { background: #ffeef0; color: #b31d28; border-color: #f9dbe1; }
        .badge-gray { background: #eff3f6; color: #24292e; border-color: #e1e4e8; }
        .badge-yellow { background: #fff8c5; color: #735c0f; border-color: #ffea7f; }
        .card.fallback { border-left: 6px solid #d4a72c; }
        pre { background: #24292e; color: #f8f8f2; padding: 15px; border-radius: 6px; overflow-x: auto; font-family: "SFMono-Regular", Consolas, "Liberation Mono", Menlo, monospace; font-size: 0.85em; line-height: 1.45; }
        details { margin-top: 10px; border: 1px solid #e1e4e8; border-radius: 6px; }
        summary { background: #f6f8fa; padding: 12px; cursor: pointer; font-weight: 600; outline: none; list-style: none; display: flex; align-items: center; }
        summary::after { content: '+'; margin-left: auto; font-weight: bold; }
        details[open] summary::after { content: '-'; }
        details[open] summary { border-bottom: 1px solid #e1e4e8; }
        .content { padding: 15px; }
        .path { font-family: monospace; color: #0366d6; word-break: break-all; background: rgba(27,31,35,0.05); padding: 2px 4px; border-radius: 3px; }
    </style>
    """

    html_content = [f"<!DOCTYPE html><html><head><meta charset='utf-8'><title>Nextflow Debug Forensics</title>{style}</head><body>"]
    html_content.append("<div class='header'><h1>🧬 Nextflow Failure Forensics</h1>")
    html_content.append(
        f"<p><strong>Total Failures:</strong> "
        f"<span style='font-size:1.5em; font-weight:bold; color:#d73a49;'>{len(failures)}</span>"
        f"&nbsp;&nbsp;<strong>Fallback Events:</strong> "
        f"<span style='font-size:1.5em; font-weight:bold; color:#b08800;'>{len(fallback_events)}</span></p></div>"
    )

    if fallback_events:
        cols = _fallback_table_columns(fallback_events)
        html_content.append("<div class='card fallback'>")
        html_content.append("<h2>🟡 Fallback Events</h2>")
        html_content.append(
            "<p style='color:#586069;'>Automatic fallbacks triggered inside tasks "
            "(emitted via <code>[FALLBACK]</code> log lines by "
            "<code>bin/run_scrublet.py</code>).</p>"
        )
        html_content.append("<table><thead><tr>")
        for c in cols:
            html_content.append(f"<th>{html.escape(c)}</th>")
        html_content.append("</tr></thead><tbody>")
        for ev in fallback_events:
            html_content.append("<tr>")
            for c in cols:
                raw_val = ev.get(c, "")
                val = html.escape(str(raw_val))
                if c == "hash" and isinstance(raw_val, str) and raw_val:
                    val = f"<code>{html.escape(raw_val[:8])}…</code>"
                if c == "type":
                    val = f"<span class='badge badge-yellow'>{val}</span>"
                html_content.append(f"<td>{val}</td>")
            html_content.append("</tr>")
        html_content.append("</tbody></table>")
        html_content.append("</div>")

    for i, fail in enumerate(failures, 1):
        process = html.escape(fail.get('process_name', 'Unknown'))
        exit_code = str(fail.get('exit_code', 'N/A'))
        task_hash = fail.get('hash', '')

        card_class = "error"
        status_text = f"Exit: {exit_code}"
        if exit_code == "137":
            card_class = "oom"
            status_text += " (OOM 💀)"

        html_content.append(f"<div class='card {card_class}'>")
        html_content.append(f"<div style='display:flex; justify-content:space-between; align-items:center;'><h3>{i}. {process}</h3> <span class='badge badge-gray'>{status_text}</span></div>")
        html_content.append(f"<p><strong>Task Hash:</strong> <code>{task_hash}</code></p>")
        html_content.append(f"<p><strong>Work Dir:</strong> <span class='path'>{fail.get('work_dir')}</span></p>")

        if fail.get('error'):
            html_content.append(f"<div style='color:#d73a49; font-weight:bold; background:#ffeef0; padding:10px; border-radius:4px;'>⚠️ {html.escape(fail.get('error'))}</div>")

        task_fallbacks = fail.get("fallback_events", [])
        if task_fallbacks:
            html_content.append("<h4>🟡 Fallback Events for this task</h4>")
            html_content.append("<ul style='margin-top:5px;'>")
            for ev in task_fallbacks:
                details = html.escape(str(ev.get("details", "")))
                html_content.append(f"<li><code>[FALLBACK] {details}</code></li>")
            html_content.append("</ul>")

        solo_sums = fail.get("solo_summary", {})
        if solo_sums:
            html_content.append("<h4>📊 STARsolo Metrics (Summary.csv)</h4>")
            for name, content in solo_sums.items():
                 html_content.append(f"<p style='margin-bottom:5px; font-weight:600; color:#586069;'>Source: {html.escape(name)}</p>")
                 html_content.append(f"<pre style='background:#f6f8fa; color:#24292e; border:1px solid #e1e4e8;'>{html.escape(content.strip())}</pre>")

        files = fail.get("file_list", [])
        if files:
            html_content.append("<h4>📂 Task Directory Files</h4>")
            html_content.append("<table><thead><tr><th>File Name</th><th>Size</th><th>Notes</th></tr></thead><tbody>")
            for f in files:
                name = html.escape(f.get('path', ''))
                size_str = f.get('size_human', 'N/A')
                msg = f.get('message', '')

                note_html = ""
                if f.get('size_bytes') == 0 and not msg:
                    note_html = "<span class='badge badge-red'>⚠️ Empty</span>"
                elif msg:
                    note_html = f"<span class='badge badge-red'>{html.escape(msg)}</span>"

                html_content.append(f"<tr><td>{name}</td><td style='white-space:nowrap;'>{size_str}</td><td>{note_html}</td></tr>")
            html_content.append("</tbody></table>")

        scripts = fail.get("related_scripts", [])
        if scripts:
            script_str = ", ".join(scripts)
            html_content.append(f"<p style='margin-top:20px;'><strong>📜 Related Scripts:</strong> <code>{script_str}</code></p>")

        logs = fail.get("logs", {})
        if logs and logs.get("content"):
            source = logs.get("source", "Log")
            html_content.append(f"<details><summary>📄 View Error Log ({source})</summary><div class='content'>")
            html_content.append(f"<pre>{html.escape(logs['content'])}</pre></div></details>")

        html_content.append("</div>")

    if script_library:
        html_content.append("<h2>📚 Source Code Library</h2>")
        for name, content in script_library.items():
            html_content.append(f"<div class='card'><details><summary>📝 {html.escape(name)} (Click to expand)</summary><div class='content'>")
            html_content.append(f"<pre>{html.escape(content)}</pre></div></details></div>")

    html_content.append("</body></html>")

    with open(output_path, "w", encoding="utf-8") as f:
        f.write("\n".join(html_content))
    print(f"Human-readable HTML report saved to: {output_path}", file=sys.stderr)


# ==========================================
# Part 4: Main Execution
# ==========================================

def main() -> None:
    args = parse_args()
    trace_path = Path(args.trace).expanduser().resolve()
    work_dir = Path(args.work_dir).expanduser().resolve()

    rows = read_trace_rows(trace_path)
    failed_rows = [row for row in rows if (row.get("status") or "").upper() == "FAILED"]

    script_library: Dict[str, str] = {}
    failures: List[Dict[str, object]] = []

    import sys
    print(f"Analyzing {len(failed_rows)} failed tasks...", file=sys.stderr)

    for row in failed_rows:
        try:
            failures.append(collect_task_report(row, work_dir, script_library))
        except Exception as exc:
            failures.append({
                "hash": row.get("hash"),
                "error": f"unexpected error: {exc}",
            })

    # Scan every task (not just FAILED) for [FALLBACK] markers emitted by
    # bin/run_scrublet.py. This surfaces fallbacks that occurred inside
    # tasks that still completed successfully.
    print(f"Scanning {len(rows)} tasks for fallback events...", file=sys.stderr)
    fallback_events: List[Dict[str, object]] = []
    for row in rows:
        task_hash = row.get("hash") or row.get("task_hash") or ""
        if not task_hash:
            continue
        task_dir = build_task_dir(work_dir, task_hash)
        try:
            events = scan_for_fallback_events(task_dir)
        except Exception:
            continue
        for ev in events:
            enriched: Dict[str, object] = dict(ev)
            enriched["process"] = row.get("process") or row.get("name", "")
            enriched["task_id"] = row.get("task_id", "")
            enriched["status"] = row.get("status", "")
            enriched["hash"] = task_hash
            enriched["work_dir"] = str(task_dir)
            fallback_events.append(enriched)

    final_payload = {
        "summary": f"Found {len(failures)} failed tasks; {len(fallback_events)} fallback events.",
        "failures": failures,
        "fallback_events": fallback_events,
        "script_library": script_library,
    }

    generate_markdown_report(final_payload, args.report)
    generate_html_report(final_payload, args.html)

    # stdout stays as clean JSON so downstream tools can pipe it directly.
    print(json.dumps(final_payload, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

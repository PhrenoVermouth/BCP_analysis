# BCP Analysis Pipeline

A Nextflow-based pipeline for single-cell RNA sequencing data preprocessing, specifically designed for BCP (Billion Cell Program) analysis. The pipeline performs comprehensive single-cell data processing including alignment, ambient RNA removal, doublet detection, and basic quality control metrics.

<img width="1436" height="1727" alt="pipeline_dag_2026-01-06_22-07-46" src="https://github.com/user-attachments/assets/e7fb4529-e5d6-42ad-831e-51797479684a" />


## Features

- **Single-cell alignment** using STAR
- **Ambient RNA removal** for cleaner expression profiles
- **Doublet detection** with automatic fallbacks for pathological samples (see below)
- **Comprehensive QC metrics** with MultiQC reporting
- **Embedded QC plots** included in MultiQC report
- **Automated preprocessing** with minimal manual intervention
- **Failure forensics**: when tasks fail or fallbacks trigger, `nf_debug_agent` collects logs, scripts, and STARsolo summaries into a single Markdown/HTML report

## Future Development

Parameters will be optimized and made configurable for different species and organ types to enhance pipeline flexibility and accuracy.

## Environment Configuration

### Current Setup
The pipeline currently requires manual conda environment configuration for testing and development purposes. 

### Planned Updates
Singularity container support will be implemented upon completion of the development phase to ensure better reproducibility and easier deployment across different computing environments.

## Installation

```bash
# Clone the repository
git clone https://github.com/PhrenoVermouth/BCP_analysis.git
cd BCP_analysis

# Create conda environment
mamba env create -f bin/environment.yml
```

## Usage

```bash
# Run the pipeline
nextflow run ~/BCP_analysis/main.nf -profile standard
```
### Run modes

- **Genefull only** (default):
  ```bash
  nextflow run ~/BCP_analysis/main.nf --run_mode genefull
  ```
- **Velocity using prior GeneFull outputs**: rerun in the same project directory after a completed GeneFull run. The pipeline
  will reuse `${outdir}/soupx/<sample>/<sample>_corrected.h5ad` by default.
  ```bash
  nextflow run ~/BCP_analysis/main.nf --run_mode velocity
  ```
  If GeneFull outputs live elsewhere, optionally add a `counts_h5ad` column in `samples.csv` to point to custom `.h5ad`
  locations.
- **Multiome (RNA + ATAC)**: runs Chromap-based ATAC alongside the RNA pipeline. Requires the multiome whitelist /
  translation files configured in `nextflow.config` (`multiome_atac_whitelist`, `multiome_rna_whitelist`,
  `multiome_barcode_translation`).
  ```bash
  nextflow run ~/BCP_analysis/main.nf --run_mode multiome
  ```

To skip SoupX ambient-RNA correction entirely (useful when SoupX struggles on a dataset and you want to feed the raw
Scrublet-cleaned matrix directly downstream), pass `--bypass_soupX true`.

### Scrublet threshold control

By default Scrublet uses its automatic threshold. Three escape hatches are available, listed in increasing specificity:

1. **Global manual override** for all samples:
   ```bash
   nextflow run ~/BCP_analysis/main.nf --scrublet_manual_threshold 0.25
   ```
   This is forwarded to `run_scrublet.py --manual_threshold`.

2. **Per-sample overrides** via `--scrublet_threshold_map /path/to/file`. Same `sample1, sample2 = value` format as
   `--mito_max_map`. Listed samples use the mapped threshold; unlisted samples fall back to
   `--scrublet_manual_threshold` (or pure automatic detection if that is also unset).

3. **Automatic fallbacks** (see next section) — applied transparently when Scrublet's auto-detection misbehaves; no
   action required.

### Automatic fallbacks in `run_scrublet.py`

The pipeline now self-corrects three pathological situations and emits a `[FALLBACK] …` log line for each one (these
lines are picked up automatically by `nf_debug_agent` and surfaced in the debug report):

| Trigger | Action | Marker |
| --- | --- | --- |
| STAR-filtered cell count > **50,000** | Keep only the top **30,000** barcodes by total UMI counts before doublet detection. Knee plot adds a red dashed reference line at the new boundary; the original black STARsolo line is retained. | `type=cell_count` |
| Scrublet's auto-detected threshold > **0.4** | Force `call_doublets(threshold=0.4)`. Doublet histogram annotates the original auto threshold with a black solid line and the fallback with a red dashed line. | `type=threshold` |
| `scrub_doublets()` returns `None` for `predicted_doublets` (auto-detection failed) | Force `call_doublets(threshold=0.4)` instead of raising. Histogram shows only the red dashed fallback line (no auto threshold to mark). | `type=none_threshold` |

Manual overrides (`--scrublet_manual_threshold` / `--scrublet_threshold_map`) always take precedence over the
threshold fallbacks; the cell-count fallback runs regardless.

### Debugging failed tasks

By default the pipeline runs `bin/nf_debug_agent.py` after the workflow completes (disable with `--debug_off true`).
It reads `trace.txt`, walks each task's work directory, and produces:

- `debug_report.md` — Markdown summary (failed tasks + fallback events table)
- `debug_report.html` — HTML version with collapsible sections and code highlighting
- A JSON payload on stdout for downstream tooling

The **Fallback Events** section lists every `[FALLBACK]` marker the agent finds — including ones from *successful*
tasks — so you can audit how many samples relied on the safety nets without digging through individual `.command.log`
files.

`nf_debug_agent` can also be run standalone for an arbitrary trace + work directory:
```bash
python bin/nf_debug_agent.py \
    --trace path/to/execution_trace_TIMESTAMP.txt \
    --work-dir path/to/work \
    --report debug_report.md \
    --html debug_report.html
```

### Mitochondrial filtering

- **Global threshold** (default): `--max_mito` sets a single mitochondrial percentage cutoff for all samples (default: `0.2`).
- **Per-sample overrides**: provide `--mito_max_map /path/to/file` where each line maps one or more sample IDs to a cutoff using
  the format `sample1, sample2 = value`. Blank lines and lines starting with `#` are ignored. Any sample not listed will fall back
  to `--max_mito`.

  Example (`resource/AC.mito`):
  ```text
  efm, em, fatfm, fatm, fbfm, fbm, hbfm, hbm, kdfm, kdm, lvfm, lvm, mbfm, mbm = 0.2
  hfm, hm, ifm, im, pcf, pcm, skfm, skm = 0.6
  lf, lm, smf, smm, spm, spfm = 0.4
  ```

## Output

The pipeline generates, under `${outdir}` (default `results/`):

- `starsolo/<sample>/` — STAR alignment, Solo.out matrices, knee plot
- `scrublet/<sample>/` — doublet-removed matrices (`*_rmdoublet.h5ad`, `*_rmMT*.h5ad`), per-sample QC PNGs, scrublet
  debug histograms, and `*_whitelist.txt`
- `soupx/<sample>/` — ambient-RNA-corrected matrices (skipped when `--bypass_soupX true`)
- `sam_qc/<sample>/` — SAM-based QC outputs and the final filtered `*_filtered_QC2.h5ad`
- `multiqc/` — MultiQC HTML report consolidating all sample-level metrics and embedded QC plots
- `debug_report.{md,html}` — failure + fallback forensics from `nf_debug_agent` (unless `--debug_off true`)

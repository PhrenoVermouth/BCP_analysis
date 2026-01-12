# BCP Analysis Pipeline

A Nextflow-based pipeline for single-cell RNA sequencing data preprocessing, specifically designed for BCP (Billion Cell Program) analysis. The pipeline performs comprehensive single-cell data processing including alignment, ambient RNA removal, doublet detection, and basic quality control metrics.

<img width="1436" height="1727" alt="pipeline_dag_2026-01-06_22-07-46" src="https://github.com/user-attachments/assets/e7fb4529-e5d6-42ad-831e-51797479684a" />


## Features

- **Single-cell alignment** using STAR
- **Ambient RNA removal** for cleaner expression profiles
- **Doublet detection** to filter multiplets
- **Comprehensive QC metrics** with MultiQC reporting
- **Embedded QC plots** included in MultiQC report
- **Automated preprocessing** with minimal manual intervention

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

The pipeline generates:
- Processed single-cell count matrices
- Quality control reports via MultiQC
- Doublet detection results
- Ambient RNA removal metrics

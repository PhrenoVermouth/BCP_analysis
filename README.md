# BCP Analysis Pipeline

A Nextflow-based pipeline for single-cell RNA sequencing data preprocessing, specifically designed for BCP (B Cell Precursor) analysis. The pipeline performs comprehensive single-cell data processing including alignment, ambient RNA removal, doublet detection, and basic quality control metrics.

<img width="1718" height="1706" alt="Drawing" src="https://github.com/user-attachments/assets/393f0ea7-708e-4563-94fe-72350aa9ce32" />

## Features

- **Single-cell alignment** using STAR
- **Ambient RNA removal** for cleaner expression profiles
- **Doublet detection** to filter multiplets
- **Comprehensive QC metrics** with MultiQC reporting
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
git clone [<repository-url>](https://github.com/PhrenoVermouth/BCP_analysis.git)
cd BCP_analysis

# Create conda environment
mamba env create -f bin/environment.yml
```

## Usage

```bash
# Run the pipeline
nextflow run ~/BCP_analysis/main.nf -profile standard
```

## Output

The pipeline generates:
- Processed single-cell count matrices
- Quality control reports via MultiQC
- Doublet detection results
- Ambient RNA removal metrics

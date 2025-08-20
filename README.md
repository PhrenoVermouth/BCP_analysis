![graphviz](https://github.com/user-attachments/assets/b18166e1-f707-41a7-b86c-576b95cca189)# BCP Analysis Pipeline

A Nextflow-based pipeline for single-cell RNA sequencing data preprocessing, specifically designed for BCP (Billion Cell Program) analysis. The pipeline performs comprehensive single-cell data processing including alignment, ambient RNA removal, doublet detection, and basic quality control metrics.

<img width="556" height="1049" alt="graphviz" src="https://github.com/user-attachments/assets/9045e386-faab-452a-99ac-b9738fee5bdc" />

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

## Output

The pipeline generates:
- Processed single-cell count matrices
- Quality control reports via MultiQC
- Doublet detection results
- Ambient RNA removal metrics

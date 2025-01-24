# NGS Analysis Project

A comprehensive workflow for:

1. **Data Download** (from Zenodo)
2. **Raw Data Alignment** (HiFi and RNA-seq)
3. **Variant Calling** (bcftools)
4. **Bulk RNA-seq Analysis** (edgeR and more)

This repository contains several shell scripts to automate each step. Below you will find **requirements**, **file overviews**, and **usage** instructions.

---

## Table of Contents
- [Requirements](#requirements)
- [Usage Summary](#usage-summary)

---
## Requirements

### System & Command-Line Tools

- **bash** (for running `.sh` scripts)
- **grep**, **find**, **cut**, **paste** (used in merging or auxiliary tasks)

### NGS & QC Tools

- **fastqc** ≥ 0.11  
- **multiqc** ≥ 1.9  
- **samtools** ≥ 1.9  
- **minimap2** ≥ 2.17  
- **qualimap** ≥ 2.2.1  
- **STAR** ≥ 2.7  
- **bcftools** ≥ 1.9  
- **gffread** for GFF to GTF conversion (if used)

### R & Packages

- **R** (≥ 3.6 recommended)
- **edgeR**, **limma**, **dplyr**, **pheatmap**, **VennDiagram**, **mixOmics**  
  (BiocManager can install edgeR, limma, etc.)

**Installation** example with conda:

```bash
conda create -n ngs-env -c bioconda -c conda-forge \
    fastqc multiqc samtools minimap2 qualimap star bcftools gffread
conda activate ngs-env
```

Then in R:
  install.packages("BiocManager")
  BiocManager::install(c("edgeR","limma","pheatmap","VennDiagram","mixOmics"))

## Usage Summary

You can run everything in sequence with:
```bash
make all
```

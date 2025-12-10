# Bulk RNA-seq Analysis Pipeline

This repository contains a customized bulk RNA-seq analysis workflow based on  
[erilu/bulk-rnaseq-analysis](https://github.com/erilu/bulk-rnaseq-analysis),  
with modifications for running on an 8 GB system using WSL/Ubuntu.

---

##  Workflow Overview
1. Raw data & preparation (Fetching SRA and converting → FASTQ
2. Quality check (FastQC & MultiQC)
3. Genome indexing & Alignment (HISAT2)
4. Read quantification (featureCounts & Count Matrix)
5. Differential expression analysis (DESeq2)
6. Visualization (MA plot, Volcano plot, PCA)

---

## 📁 Repository Structure

⚠️ Note for low-RAM systems
Some alignment steps (HISAT2) may fail on 8 GB RAM.  
Increasing swap space (16 GB) is recommended.

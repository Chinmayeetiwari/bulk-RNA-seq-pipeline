# Bulk RNA-seq Analysis Pipeline

This repository documents my complete learning journey in Bulk RNA-seq analysis, starting from raw sequencing files and building a fully functional analysis pipeline step-by-step.

As a biotechnology student exploring computational workflows, I created this project to understand the practical implementation of RNA-seq from data conversion and quality assessment to alignment, quantification, and differential gene expression analysis using DESeq2.

This repository contains a customized bulk RNA-seq analysis workflow based on [erilu/bulk-rnaseq-analysis](https://github.com/erilu/bulk-rnaseq-analysis) The pipeline is optimized for running on an 8 GB RAM system (WSL Ubuntu), which makes it accessible for students and beginners who may not have access to high-performance computing.

This project serves both as my learning record and as a reference workflow for others beginning RNA-seq analysis.

---
##  Workflow Overview
1. Raw data & preparation (Fetching SRA and converting to FASTQ)
2. Quality check (FastQC & MultiQC)
3. Genome indexing & Alignment (HISAT2)
4. Read quantification (featureCounts & Count Matrix)
5. Differential expression analysis (DESeq2)
6. Visualization (MA plot, Volcano plot, PCA)

---
**Raw data & preparation (Fetching SRA and converting to FASTQ)**

## 📁 Repository Structure

⚠️ Note for low-RAM systems
Some alignment steps (HISAT2) may fail on 8 GB RAM.  
Increasing swap space (16 GB) is recommended.

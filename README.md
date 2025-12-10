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
## Raw data & preparation (Fetching SRA and converting to FASTQ)**

**Samples from GEO (GSE106305)**

I retrieved RNA-seq samples from GEO accession **GSE106305**, containing LNCaP and PC3 prostate cancer cell lines under normoxia and hypoxia conditions. More information about sequencing and sample prep can be obtained from Gene Expression Omnibus (GEO).

The samples used in this study are listed below:

| Sample Name                                    | GSM Identifier | SRA Identifier (SRX) | SRA Runs (SRR — download these) |
|------------------------------------------------|----------------|----------------------|------------------------------------------------|
| LNCaP_RNA-Seq_Empty_Vector_Normoxia_rep1       | GSM3145509     | SRX4096735           | SRR7179504, SRR7179505, SRR7179506, SRR7179507 |
| LNCaP_RNA-Seq_Empty_Vector_Normoxia_rep2       | GSM3145510     | SRX4096736           | SRR7179508, SRR7179509, SRR7179510, SRR7179511 |
| LNCaP_RNA-Seq_Empty_Vector_Hypoxia_rep1        | GSM3145513     | SRX4096739           | SRR7179520, SRR7179521, SRR7179522, SRR7179523 |
| LNCaP_RNA-Seq_Empty_Vector_Hypoxia_rep2        | GSM3145514     | SRX4096740           | SRR7179524, SRR7179525, SRR7179526, SRR7179527 |
| PC3_RNA-Seq_siCtrl_Normoxia_rep1               | GSM3145517     | SRX4096743           | SRR7179536 |
| PC3_RNA-Seq_siCtrl_Normoxia_rep2               | GSM3145518     | SRX4096744           | SRR7179537 |
| PC3_RNA-Seq_siCtrl_Hypoxia_rep1                | GSM3145521     | SRX4096747           | SRR7179540 |
| PC3_RNA-Seq_siCtrl_Hypoxia_rep2                | GSM3145522     | SRX4096748           | SRR7179541 |

**Fetching Data from SRA**
1. Install SRA Toolkit
To download sequencing data from the NCBI Sequence Read Archive (SRA), we use the SRA Toolkit.
It provides tools like prefetch (downloads .sra files) and fastq-dump / fasterq-dump (converts to FASTQ format).

```
sudo apt install sra-toolkit
```
2. Downloading One SRR File (Manual Method)

Use prefetch to download; this will download the .sra file into your current directory.

```
prefetch SRR7179504
```
Convert ``` .sra ``` to fastq files
   
```
fastq-dump --outdir fastq --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip SRR7179504.sra
```
This creates: 
```
fastq/SRR7179504_1.fastq.gz
fastq/SRR7179504_2.fastq.gz
```
3. Automating Downloads for Multiple SRR Ids
   
For automation, first check ```python --version``` then create a Python script that reads a list of SRR IDs and runs prefetch + fastq-dump command for each sample downloaded and logs progress.

To run this script in bash
```
python3 bulk_prefetch+fastq_dump.py
```

fastqc fastq/*.fastq.gz -o fastqc_results/ --threads 8
## 📁 Repository Structure

⚠️ Note for low-RAM systems
Some alignment steps (HISAT2) may fail on 8 GB RAM.  
Increasing swap space (16 GB) is recommended.

# Bulk RNA-seq Analysis Pipeline

This repository documents my complete learning journey in Bulk RNA-seq analysis, starting from raw sequencing files and building a fully functional analysis pipeline step-by-step.

As a biotechnology student exploring computational workflows, I created this project to understand the practical implementation of RNA-seq from data conversion and quality assessment to alignment, quantification, and differential gene expression analysis using DESeq2.

This repository contains a customized bulk RNA-seq analysis workflow based on [erilu/bulk-rnaseq-analysis](https://github.com/erilu/bulk-rnaseq-analysis) The pipeline is optimized for running on an 8 GB RAM system (WSL Ubuntu), which makes it accessible for students and beginners who may not have access to high-performance computing.

This project serves both as my learning record and as a reference workflow for others beginning RNA-seq analysis.

---
##  Workflow Overview
1. Raw data retrieval & Fastq Conversion
2. Quality Control (FastQC & MultiQC)
3. Genome indexing & Alignment (HISAT2)
4. Read quantification (featureCounts & Count Matrix)
5. Differential expression analysis (DESeq2)
6. Visualization (MA plot, Volcano plot, PCA)

---
##  Raw data retrieval & Fastq Conversion

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

** Fetching Data from SRA**
1.Install SRA Toolkit
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
3. Convert ``` .sra ``` to fastq files
   
```
fastq-dump --outdir fastq --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip SRR7179504.sra
```
This creates: 
```
fastq/SRR7179504_1.fastq.gz
fastq/SRR7179504_2.fastq.gz
```
Automating Downloads for Multiple SRR Ids
   
For automation, first check ```python --version``` then create a Python script that reads a list of SRR IDs and runs prefetch + fastq-dump command for each sample downloaded and logs progress.

To run this script in bash
```
python3 bulk_prefetch+fastq_dump.py
```

You can find the automation script here: **[`scripts/download_sra_fastq_dump.py`](scripts/download_sra_fastq_dump.py)`**
 
 ## Quality Control (FastQC & MultiQC)
**FastQC**<br>
FastQC is a widely used tool for assessing the quality of raw sequencing data (FASTQ files). It generates visual reports that summarize key quality metrics, including per-base sequence quality, GC content, adapter contamination, sequence length distribution, overrepresented sequences, and duplication levels. These metrics provide an overview of the integrity of the sequencing data and help determine whether the reads are suitable for downstream alignment or if trimming and preprocessing steps are required.

Run FastQC
```
fastqc fastq/*.fastq.gz -o fastqc_results/ --threads 8
```
This command processes all FASTQ files in the fastq/ directory and saves the reports (.html and .zip) in fastqc_results/

**MultiQC**<br>
When you have many samples, FastQC produces many individual reports. MultiQC solves this by gathering all FastQC results and summarizing them into a single, easy-to-read report.

MultiQC is extremely useful for:

-Comparing quality across replicates<br>
-Checking batch effects<br>
-Spotting problematic samples

Run MultiQC
```
multiqc fastqc_results/ -o multiqc_report/
```
This command:<br>

Scans the fastqc_results directory and generates a unified multiqc_report.html and saves it inside multiqc_report folder.

**Read Trimming (Trimmomatic)**

If FastQC reports poor quality scores, adapter contamination, or low-quality trailing bases, trimming should be performed before proceeding with alignment. Trimmomatic is used for quality and adapter trimming. The following command was used to trim a single-end sample:
```
java -jar Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 4 \
fastq/SRR7179504.fastq.gz fastq/SRR7179504_trimmed.fastq.gz \
TRAILING:10 -phred33
```
After trimming, FastQC should be rerun to verify improvement in read quality:

## Genome Indexing and Alignment

HISAT2 was used for genome alignment, and the GRCh38 pre-built index was downloaded from the official HISAT2 repository. The index was extracted and used for aligning the processed FASTQ files. HISAT2 and Samtools must be installed before running alignment.

Download GRCh38 Index
```
wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
tar -xvzf grch38_genome.tar.gz
```

Install Required Tools
```
sudo apt install hisat2
sudo apt install samtools
```
Manual Alignment Command (Single Sample)
```hisat2 -q -x grch38/genome -U fastq/sample1.fastq.gz \
| samtools sort -o sample1.bam \
| samtools index sample1.bam
```
Automated Alignment Script

An automated version of this workflow is provided in: [`scripts/hisat2a.sh`] <br>This script processes multiple FASTQ files and generates sorted and indexed BAM files for all samples.

⚠️ Note for low-RAM systems
Some alignment steps (HISAT2) may fail on 8 GB RAM.  
Increasing swap space (16 GB) is recommended.

## Read quantification (featureCounts & Count Matrix) 
After alignment, read quantification was performed using featureCounts from the Subread package. This step assigns aligned reads to genomic features (genes) based on the GTF annotation file and produces a table of raw counts for each sample. Gene-level counts are required for downstream differential expression analysis. The latest annotation version (Homo_sapiens.GRCh38.114.gtf, Ensembl release) was downloaded for this workflow.

Install featureCounts
```
sudo apt install subread
```
Manual Command for Quantifying One BAM File
```
featureCounts -S 2 -a Homo_sapiens.GRCh38.114.gtf \
-o quants/featurecounts.txt sample.bam
```

Automated FeatureCounts Script

An automated script for generating counts for multiple BAM files is available in: [`scripts/featurecount.sh`] <br> This script processes all aligned and sorted BAM files and outputs a single featureCounts table. 
Run it with:

```
chmod +x scripts/featurecount.sh
./scripts/featurecount.sh
```

**Count Matrix Generation**

The raw featureCounts output is reformatted into a clean count matrix for downstream analysis. This step is handled by:[`scripts/countmatrix.sh`]

Run the script to generate the consolidated count matrix:
```
chmod +x scripts/countmatrix.sh
./scripts/countmatrix.sh
```
The resulting matrix contains gene identifiers as rows and samples as columns, and is used directly in the DESeq2 workflow.

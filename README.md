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

## **Concatenating FASTQ files**
Some samples in this dataset (LNCaP) were split across multiple SRA runs, resulting in several FASTQ files per sample after conversion. These runs must be concatenated into a single FASTQ file for each biological sample. The PC3 samples consist of only one SRA run each, so they are simply renamed.
Concatenate FASTQ files for LNCaP samples
```
cat SRR7179504_pass.fastq.gz SRR7179505_pass.fastq.gz SRR7179506_pass.fastq.gz SRR7179507_pass.fastq.gz > LNCAP_Normoxia_S1.fastq.gz
cat SRR7179508_pass.fastq.gz SRR7179509_pass.fastq.gz SRR7179510_pass.fastq.gz SRR7179511_pass.fastq.gz > LNCAP_Normoxia_S2.fastq.gz
cat SRR7179520_pass.fastq.gz SRR7179521_pass.fastq.gz SRR7179522_pass.fastq.gz SRR7179523_pass.fastq.gz > LNCAP_Hypoxia_S1.fastq.gz
cat SRR7179524_pass.fastq.gz SRR7179525_pass.fastq.gz SRR7179526_pass.fastq.gz SRR7179527_pass.fastq.gz > LNCAP_Hypoxia_S2.fastq.gz
```
Rename PC3 FASTQ files
```
mv SRR7179536_pass.fastq.gz PC3_Normoxia_S1.fastq.gz
mv SRR7179537_pass.fastq.gz PC3_Normoxia_S2.fastq.gz
mv SRR7179540_pass.fastq.gz PC3_Hypoxia_S1.fastq.gz
mv SRR7179541_pass.fastq.gz PC3_Hypoxia_S2.fastq.gz
```

After concatenation and renaming, the directory should contain eight finalized FASTQ files (four LNCaP and four PC3). The individual SRR FASTQ files may be removed:
```
rm SRR*
```
These FASTQ files are now ready for quality control and downstream alignment.
 
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

## **Differential expression analysis**
After generating the final gene-level count matrix, the next step in the workflow is differential expression analysis. This analysis identifies genes that show significant changes in expression between conditions such as normoxia and hypoxia, or between different cell lines such as LNCaP and PC3. DESeq2 was used for this analysis because it provides robust normalization, variance estimation, and statistical testing for RNA-seq count data.

The input to DESeq2 is the processed count matrix produced after running featureCounts and the count matrix formatting script. Using this matrix along with a sample metadata file, DESeq2 performs normalization, estimates dispersion, fits the statistical model, and generates results tables that include fold changes, adjusted p-values, and significance indicators. The analysis also produces exploratory visualizations such as PCA plots and sample clustering, which help assess sample relationships and overall dataset quality.

**Installing and Loading Required Packages**

Before beginning differential expression analysis, the necessary Bioconductor packages must be installed and loaded. DESeq2 is the primary package used for normalization, dispersion estimation, and statistical testing of count-based RNA-seq data. If BiocManager is not already installed, it is added first, followed by the installation of DESeq2. Once installed, the library is loaded into the R session.
```
getwd()
setwd

install.packages("BiocManager")
BiocManager::install(version = "devel")
BiocManager::install("DESeq2")

library(DESeq2)
```

**Preparing the Count Matrix and Metadata
**
The next step involves importing the raw gene count matrix and creating the metadata (sample information) required by DESeq2. The count matrix contains gene identifiers as rows and sample names as columns. Metadata includes variables describing each sample — here, the experimental condition for each replicate (e.g., normoxia vs. hypoxia for both cell lines). The sample order is matched correctly by sorting the column names and assigning conditions in the same order.
```
library(dplyr)
library(tibble)

raw_counts <- read.csv("GSE106305_counts_matrix.csv", header = TRUE, row.names = "Geneid", stringsAsFactors = FALSE)
head(raw_counts)
raw_counts <- raw_counts[,sort(colnames(raw_counts))]
colSums(raw_counts)

condition <- c(rep("LNCAP_Hypoxia", 2), rep("LNCAP_Normoxia", 2), rep("PC3_Hypoxia", 2), rep("PC3_Normoxia", 2))
print(condition)

my_colData <- as.data.frame(condition)
rownames(my_colData) <- colnames(raw_counts)
head(my_colData)
```
**Creating the DESeq2 Dataset Object**

DESeq2 requires its data to be stored in a special object called a DESeqDataSet, which combines the count matrix with the associated sample metadata. The design formula specifies the variable of interest — in this case, the experimental condition. After creating the dataset, the counts can be inspected to confirm that the data has been loaded correctly.
```dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                              colData = my_colData,
                              design = ~condition)
dds
head(counts(dds))
dim(counts(dds))
```
**Basic Count Matrix Inspection**

A quick check is performed to understand the distribution of zeros across genes. Many genes may not be expressed in all samples, and summarizing the number of zero-count genes helps evaluate dataset sparsity.
```
count_matrix <- counts(dds)
dim(count_matrix)
zero_counts_per_gene <- rowSums(count_matrix == 0)
count_matrix <- as.data.frame(count_matrix)
zero_summary <- table(zero_counts_per_gene)
print(zero_summary)
```
**Adding Gene Annotation and Filtering Genes**

After loading the count matrix and creating the DESeq2 dataset, gene annotation is incorporated to associate each Ensembl gene ID with its corresponding gene symbol and biotype. The annotation file is imported and merged with the count matrix using a common gene identifier. To ensure proper matching, version numbers are removed from the Ensembl IDs in both the annotation file and the count matrix before merging. The merged table includes the gene ID, gene symbol, biotype, and all sample columns.
```library(data.table)
annotation_file <- "GRCh38annotation.csv"
annotation <- fread(annotation_file, stringsAsFactors = FALSE)
counts_gse <- read.csv("GSE106305_counts_matrix.csv",
                     header = TRUE,
                     stringsAsFactors = FALSE)

counts_gse$Geneid <- sub("\\..*$", "", counts_gse$Geneid)
annotation$Geneid <- sub("\\..*$", "", annotation$Geneid)

annotated_counts <- left_join(counts_gse, annotation, by = "Geneid") %>%
  select(Geneid, Genesymbol, Genebiotype, 
         LNCAP_Hypoxia_S1, LNCAP_Hypoxia_S2, LNCAP_Normoxia_S1, LNCAP_Normoxia_S2, 
         PC3_Hypoxia_S1, PC3_Hypoxia_S2, PC3_Normoxia_S1, PC3_Normoxia_S2)
```
**Filtering Genes by Biotype**

Only genes belonging to biologically relevant categories (such as protein-coding and immune gene classes) are retained. All other biotypes are removed to focus the analysis on meaningful transcriptional signals. The filtered count table is saved for reference.
```
biotypes_to_keep <- c("protein_coding", "IG_J_gene", "IG_V_gene", "IG_C_gene", "IG_D_gene",
                      "TR_D_gene", "TR_C_gene", "TR_V_gene", "TR_J_gene")

filtered_counts <- annotated_counts %>%
  filter(Genebiotype %in% biotypes_to_keep)

filtered_counts$Geneid <- sub("\\..*$", "", filtered_counts$Geneid)
head(filtered_counts, n = 3)

output_file <- "9biotype_count_matrix.csv"
fwrite(filtered_counts, file = output_file, sep = ",", row.names = FALSE)

zero_counts1 <- rowSums(filtered_counts[, 4:11] == 0)
zero_summary2 <- table(zero_counts1)
print(zero_summary2)
```
**Filtering Low-Abundance Genes**

To reduce noise, genes with zero counts in most samples are removed. Genes expressed in fewer than seven samples are excluded, and the filtered table is exported. The DESeq2 object is then subset accordingly so that only the retained genes proceed to further analysis.
```
keep_genes <- zero_counts1 < 7
filtered_counts_nozero <- filtered_counts[keep_genes, ]
cat("Number of genes after filtering (zeros in <7 samples):", nrow(filtered_counts_nozero), "\n")

new_zero_counts <- rowSums(filtered_counts_nozero[, 4:11] == 0)
cat("New zero counts distribution:\n")
print(table(new_zero_counts))

output_file <- "filtered_biotype_nozero_count_matrix.csv"
fwrite(filtered_counts_nozero, file = output_file, sep = ",", row.names = FALSE)

head(filtered_counts_nozero, n = 3)

dds_filtered <- dds[rownames(dds) %in% filtered_counts_nozero$Geneid, ]
cat("Dimensions of filtered DESeqDataSet:", dim(dds_filtered), "\n")

removed_genes <- filtered_counts[!keep_genes, ]
cat("Biotype distribution of removed genes:\n")
print(table(removed_genes$Genebiotype))
```
**Final Gene Filtering and Creating the Filtered DESeq2 Object**

After annotating genes and applying biotype filtering, an additional filtering step is performed to remove genes with very low expression across most samples. Genes with zero counts in seven or more samples are excluded, ensuring that downstream statistical analysis focuses on genes with sufficient signal. The resulting filtered count matrix is exported, and the DESeq2 dataset is subset to include only the retained genes.
```
print(colnames(filtered_counts))
zero_counts <- rowSums(filtered_counts[, 4:11] == 0)
zero_summary <- table(zero_counts)
print(zero_summary)

keep_genes <- zero_counts < 7
filtered_counts_nozero <- filtered_counts[keep_genes, ]
cat("Number of genes after filtering (zeros in <7 samples):", nrow(filtered_counts_nozero), "\n")

new_zero_counts <- rowSums(filtered_counts_nozero[, 4:11] == 0)
print(table(new_zero_counts))

output_file <- "filtered_biotype_6.csv"
fwrite(filtered_counts_nozero, file = output_file, sep = ",", row.names = FALSE)

head(filtered_counts_nozero, n = 3)

dds_filtered <- dds[rownames(dds) %in% filtered_counts_nozero$Geneid, ]
```
**Exploring Gene Biotype Composition After Filtering**

To understand the impact of filtering on gene categories, the distribution of biotypes in the remaining dataset is summarized and visualized. This bar plot shows the proportion of each retained biotype and helps verify that biologically relevant gene classes are preserved.
```   
library(dplyr)

biotype_counts <- filtered_counts_nozero %>%
  dplyr::count(Genebiotype) %>%
  dplyr::mutate(Proportion = n / sum(n),
                Percentage = Proportion * 100) %>%
  dplyr::rename(Biotype = Genebiotype)

print(biotype_counts)

library(ggplot2)
p <- ggplot(biotype_counts, aes(x = reorder(Biotype, -Proportion), y = Proportion, fill = Biotype)) +
  geom_bar(stat = "identity") +
  labs(title = "Proportion of Genes by Biotype",
       x = "Gene Biotypes",
       y = "Proportion") +
  scale_y_continuous(labels = scales::percent_format(scale = 100)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  scale_fill_brewer(palette = "Set2")

p

output_plot <- "genebiotype_proportions1.png"
ggsave(output_plot, plot = p, width = 8, height = 6, dpi = 300)
```
**Variance Stabilization and PCA for Sample Exploration**

A variance-stabilizing transformation (VST) is applied to the filtered DESeq2 dataset to normalize variance across expression levels, making the data suitable for dimensionality reduction methods such as PCA. The PCA plot illustrates how samples cluster based on their expression profiles and provides an overview of major sources of variation, such as cell line differences or hypoxia treatment.
vsd <- vst(dds_filtered, blind = TRUE)
```
plot_PCA = function (vsd.obj) {
  pcaData <- plotPCA(vsd.obj, intgroup = c("condition"), returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  ggplot(pcaData, aes(PC1, PC2, color = condition)) +
    geom_point(size = 3) +
    labs(x = paste0("PC1: ", percentVar[1], "% variance"),
         y = paste0("PC2: ", percentVar[2], "% variance"),
         title = "PCA Plot colored by condition") +
    ggrepel::geom_text_repel(aes(label = name), color = "black")
}

png(filename = "pcab.png", width = 2000, height = 2000, res = 300)
plot_PCA(vsd)
dev.off()
```

**Running the DESeq2 Model and Extracting Normalized Counts**

Once low-expression genes have been removed, the DESeq2 model is run on the filtered dataset. This step performs normalization, dispersion estimation, and model fitting based on the specified design. Normalized counts are extracted and saved for downstream visualization or exploratory analysis.
```
dds <- DESeq(dds_filtered)
dds

normalized_counts <- counts(dds, normalized = TRUE)
normalized_counts_df <- as.data.frame(normalized_counts)
write.csv(normalized_counts_df, file = "normalized_counts.csv", row.names = TRUE)
```
**Variance Stabilization and PCA Visualization**

Variance-stabilizing transformation (VST) is applied to the dataset to reduce heteroscedasticity and make expression values comparable across different expression levels. PCA is then used to visualize sample relationships and major sources of variation across conditions and cell lines. The resulting PCA plot is exported as an image.
```
vsd <- vst(dds, blind = TRUE)

plot_PCA = function (vsd.obj) {
  pcaData <- plotPCA(vsd.obj, intgroup = c("condition"), returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  ggplot(pcaData, aes(PC1, PC2, color = condition)) +
    geom_point(size = 3) +
    labs(x = paste0("PC1: ", percentVar[1], "% variance"),
         y = paste0("PC2: ", percentVar[2], "% variance"),
         title = "PCA Plot colored by condition") +
    ggrepel::geom_text_repel(aes(label = name), color = "black")
}

png(filename = "pcaa.png", width = 2000, height = 2000, res = 300)
plot_PCA(vsd)
dev.off()
```
**Sample Distance Heatmap**

A heatmap of pairwise sample distances is generated to assess similarities and detect any potential outliers. This visualization complements PCA by showing hierarchical relationships among samples based on global expression patterns.

Variance-stabilizing transformation (VST) is applied to the dataset to reduce heteroscedasticity and make expression values comparable across different expression levels. PCA is then used to visualize sample relationships and major sources of variation across conditions and cell lines. The resulting PCA plot is exported as an image.
```
vsd <- vst(dds, blind = TRUE)

plotDists = function (vsd.obj) {
  sampleDists <- dist(t(assay(vsd.obj)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(vsd.obj$condition)
  colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(55)
  pheatmap::pheatmap(sampleDistMatrix,
                     clustering_distance_rows = sampleDists,
                     clustering_distance_cols = sampleDists,
                     col = colors,
                     fontsize_row = 4, fontsize_col = 4, fontsize_legend = 4, fontsize = 4)
}

png(filename = "sampleheatmap1.png", width = 1000, height = 900, res = 300)
plotDists(vsd)
dev.off()
```
**Heatmap of Top Variable Genes**

To highlight genes contributing most strongly to differences between samples, the most variable genes across the dataset are identified and visualized using a heatmap. Gene names are added by matching Ensembl IDs to the annotation file. This visualization provides insight into biological patterns and clustering by condition or cell line.
```
variable_gene_heatmap <- function (vsd.obj, num_genes = 500, annotation, title = "") {
  brewer_palette <- "RdBu"
  ramp <- colorRampPalette(RColorBrewer::brewer.pal(11, brewer_palette))
  mr <- ramp(256)[256:1]

  stabilized_counts <- assay(vsd.obj)
  row_variances <- rowVars(stabilized_counts)
  top_variable_genes <- stabilized_counts[order(row_variances, decreasing = TRUE)[1:num_genes],]
  top_variable_genes <- top_variable_genes - rowMeans(top_variable_genes, na.rm = TRUE)

  gene_names <- annotation$Gene.name[match(rownames(top_variable_genes), annotation$Gene.stable.ID)]
  rownames(top_variable_genes) <- gene_names

  coldata <- as.data.frame(vsd.obj@colData)
  coldata$sizeFactor <- NULL

  pheatmap::pheatmap(top_variable_genes,
                     color = mr,
                     annotation_col = coldata,
                     fontsize_col = 8,
                     fontsize_row = 250/num_genes,
                     border_color = NA,
                     main = title)
}

png(filename = "variable_gene_heatmap.png", width = 1000, height = 1000, res = 300)
variable_gene_heatmap(vsd, num_genes = 40, annotation = annotation)
dev.off()

dim(dds)
```

#!/bin/bash

# Directory containing FASTQ files
FASTQ_DIR=~/rna_seq_fastq_files/fastq_output
GENOME_INDEX=~/rna_seq_fastq_files/grch38/genome
LOGFILE="alignment_log.txt"

# Clear or create logfile
> $LOGFILE

# List of FASTQ files
FILES=(
    "LNCAP_Hypoxia_S1.fastq.gz"
    "LNCAP_Hypoxia_S2.fastq.gz"
    "LNCAP_Normoxia_S1.fastq.gz"
    "LNCAP_Normoxia_S2.fastq.gz"
    "PC3_Hypoxia_S1.fastq.gz"
    "PC3_Hypoxia_S2.fastq.gz"
    "PC3_Normoxia_S1.fastq.gz"
    "PC3_Normoxia_S2.fastq.gz"
)

# Number of threads for HISAT2 and Samtools
THREADS=2

# Loop through each file
for f in "${FILES[@]}"; do
    SAMPLE_NAME=$(basename "$f" .fastq.gz)
    echo "Processing $SAMPLE_NAME at $(date)" | tee -a $LOGFILE
    START_TIME=$(date +%s)

    # Run HISAT2 alignment and Samtools sorting/indexing
    hisat2 -p $THREADS -q -x $GENOME_INDEX -U $FASTQ_DIR/$f 2>> $LOGFILE | \
    samtools sort -@ $THREADS -o ${SAMPLE_NAME}.bam
    samtools index ${SAMPLE_NAME}.bam

    END_TIME=$(date +%s)
    ELAPSED=$((END_TIME - START_TIME))
    echo "Finished $SAMPLE_NAME in $ELAPSED seconds at $(date)" | tee -a $LOGFILE
    echo "--------------------------------------" | tee -a $LOGFILE

done

echo "All files processed successfully at $(date)" | tee -a $LOGFILE

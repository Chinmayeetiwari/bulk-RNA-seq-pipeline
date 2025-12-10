#!/bin/bash

FASTQ_DIR=~/rna_seq_fastq_files/fastq_output
GENOME_INDEX=~/rna_seq_fastq_files/grch38/genome
LOGFILE=~/rna_seq_fastq_files/PC3_Normoxia_S2_rerun.log

# Clear or create logfile
> $LOGFILE

SAMPLE="PC3_Normoxia_S2.fastq.gz"
SAMPLE_NAME="PC3_Normoxia_S2"

echo "====================================" | tee -a $LOGFILE
echo "Starting alignment for $SAMPLE_NAME at $(date)" | tee -a $LOGFILE
START_TIME=$(date +%s)

# Run HISAT2 + Samtools Sort
hisat2 -p 2 -q -x $GENOME_INDEX -U $FASTQ_DIR/$SAMPLE \
    2>> $LOGFILE \
| samtools sort -@ 2 -m 500M -o ${SAMPLE_NAME}.bam \
    2>> $LOGFILE

# Check if BAM was created
if [[ ! -f ${SAMPLE_NAME}.bam ]]; then
    echo "ERROR: BAM file was not created. Sorting or alignment failed." | tee -a $LOGFILE
    exit 1
fi

# Index BAM
samtools index ${SAMPLE_NAME}.bam 2>> $LOGFILE

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

echo "Finished $SAMPLE_NAME successfully at $(date)" | tee -a $LOGFILE
echo "Total time: $ELAPSED seconds" | tee -a $LOGFILE
echo "====================================" | tee -a $LOGFILE


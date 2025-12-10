#!/bin/bash

# Go to your aligned reads folder
ANNOTATION=~/rna_seq_fastq_files/Homo_sapiens.GRCh38.114.gtf
QUANTS_DIR=~/rna_seq_fastq_files/quants2

# Create output folder if not present
mkdir -p $QUANTS_DIR

# Loop over all BAM files
for bam in *.bam; do
    start=$(date +%s)  # start time

    echo "Processing $bam ..."
    featureCounts -s 2 -a $ANNOTATION \
        -o $QUANTS_DIR/${bam%.bam}_featurecounts.txt \
        "$bam"

    end=$(date +%s)  # end time
    runtime=$(( (end - start) / 60 ))  # in minutes

    echo "Completed $bam in $runtime minutes."
    echo "------------------------------------"
done

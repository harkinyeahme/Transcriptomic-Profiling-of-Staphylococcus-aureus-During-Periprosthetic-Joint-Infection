#!/bin/bash
# Create directory for raw reads quality control with fastqc
mkdir -p qc_raw
# Loop over each paired end 
for R1 in raw_reads/*_1.fastq.gz; do
SAMPLE=$(basename "$R1" _1.fastq.gz)
R2=raw_reads/${SAMPLE}_2.fastq.gz
fastqc "$R1" "$R2" -o qc_raw
done
# Generate consolidated MultiQC reports for all raw reads
multiqc qc_raw/ -o multiqc_raw

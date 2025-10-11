#!/bin/bash
# Create directory for trimmed data
mkdir -p trimmed_reads
# Loop over each paired end data
for R1 in raw_reads/*_1.fastq.gz; do
SAMPLE=$(basename "$R1" _1.fastq.gz)
R2=raw_reads/${SAMPLE}_2.fastq.gz
# Perform quality trimming and adapter trimming using fastp
fastp \
	-i "$R1" \
	-I "$R2" \
	-o trimmed_reads/${SAMPLE}_1_trimmed_fastq.gz \
	-O trimmed_reads/${SAMPLE}_2_trimmed_fastq.gz \
	-h trimmed_reads/${SAMPLE}_trimmed.html \
	-j trimmed_reads/${SAMPLE}_trimmed.json
done
# Create directories for the trimmed data quality control 
mkdir -p qc_trimmed multiqc_trimmed
# Loop through all the trimmed paired end data
for R1 in trimmed_reads/*1_trimmed_fastq.gz; do
SAMPLE=$(basename "$R1" _1_trimmed_fastq.gz)
R2=trimmed_reads/${SAMPLE}_2_trimmed_fastq.gz

# Perform Quality control on the trimmed paired end data
fastqc "$R1" "$R2" -o qc_trimmed
done
# Generate consolidated MultiQC reports for all trimmed data
multiqc qc_trimmed/ -o multiqc_trimmed

#!/bin/bash
# ============================================================
# RNA-Seq Mapping Pipeline using STAR and Samtools
# Organism: Staphylococcus aureus subsp. aureus NCTC8325
# Input: Paired-end trimmed FASTQ files
# Output: Sorted and indexed BAM files for each sample
# ============================================================

# --- Step 1: Create a folder for the reference genome
mkdir -p reference && cd reference/

# --- Step 2: Download reference genome (FASTA)
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.fna.gz

# --- Step 3: Rename for simplicity and unzip the file
mv GCF_000013425.1_ASM1342v1_genomic.fna.gz S_aureus.fna.gz
gunzip S_aureus.fna.gz

# --- Step 4: Download corresponding GFF annotation
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.gff.gz

# --- Step 5: Rename and unzip the GFF annotation file
mv GCF_000013425.1_ASM1342v1_genomic.gff.gz S_aureus.gff.gz
gunzip S_aureus.gff.gz

# --- Step 6: Go back to main project directory
cd ..

# --- Step 7: Create STAR genome index directory
mkdir -p genome

# --- Step 8: Generate STAR genome index
STAR --runThreadN 8 \
--runMode genomeGenerate \
--genomeDir genome \
--genomeFastaFiles reference/S_aureus.fna \
--limitGenomeGenerateRAM 80000000000

# --- Step 10: Create directory for output BAM files
mkdir -p mapped_reads

# --- Step 11: Loop through all paired-end FASTQ files
for R1 in trimmed_reads/*_1_trimmed_fastq.gz; do
SAMPLE=$(basename "$R1" _1_trimmed_fastq.gz)
R2=trimmed_reads/${SAMPLE}_2_trimmed_fastq.gz
echo " Processing sample: $SAMPLE"

# --- Step 12: Align reads using STAR
STAR --runThreadN 8 \
--genomeDir genome \
--readFilesIn "$R1" "$R2" \
--readFilesCommand zcat \
--outFileNamePrefix mapped_reads/${SAMPLE}_ \
--outSAMtype BAM SortedByCoordinate

# --- Step 13: Rename and index BAM file
mv "mapped_reads/${SAMPLE}_Aligned.sortedByCoord.out.bam" "mapped_reads/${SAMPLE}_sorted.bam"
samtools index "mapped_reads/${SAMPLE}_sorted.bam"

echo " Completed: $SAMPLE"
done

# --- Step 14: Prepare files for IGV visualization
mkdir -p IGV
cp mapped_reads/*sorted.bam IGV/
cp mapped_reads/*sorted.bam.bai IGV/

# ============================================================
# OUTPUTS:
#   mapped_reads/<sample>_sorted.bam      → Sorted BAM file
#   mapped_reads/<sample>_sorted.bam.bai  → BAM index
# These can be used in IGV or featureCounts.
# ============================================================


#!/bin/bash
# Create a directory for the counts
mkdir -p counts
# Perform feature counts on the sorted.bam file using the reference annotation file
featureCounts -p -B -C -O -t gene -g locus_tag -a reference/S_aureus.gff -o counts/count.txt IGV/*_sorted.bam

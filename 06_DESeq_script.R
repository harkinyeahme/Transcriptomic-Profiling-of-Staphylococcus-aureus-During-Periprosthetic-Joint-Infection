#############################################
# Differential Expression Analysis (DESeq2)
# Project: Staphylococcus aureus (Acute vs Chronic PJI)
# Description:
# This script performs RNA-seq differential analysis using DESeq2.
# It includes data cleaning, normalization, visualization and exports results safely.
#############################################

# --- 1. Install Required Packages ---
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("DESeq2")) BiocManager::install("DESeq2")
if (!require("pheatmap")) install.packages("pheatmap")
if (!require("ggplot2")) install.packages("ggplot2")

# --- 2. Load Libraries ---
library(DESeq2)
library(pheatmap)
library(ggplot2)

# --- 3. Set Working Directory (Change this path to your own folder) ---
setwd("C:/Users/Olanrewaju/Downloads/")   # <--- Update this to your local path

# --- 4. Load Data ---
# Read featureCounts output
counts <- read.delim("count.txt", header = TRUE, comment.char = "#", sep = "\t")
# Read metadata file
s_a_meta  <- read.csv("metadata.csv", header = TRUE)

# --- 5. Remove non-gene columns such as genomic coordinates.
raw_counts <- counts[, -(1:6)]
rownames(raw_counts) <- counts$Geneid  # Use gene IDs as row names

# --- 6. Clean Sample Names ---
# Removes unnecessary prefixes and suffixes added by aligners or IGV.
colnames(raw_counts) <- gsub("^IGV\\.", "", colnames(raw_counts))
colnames(raw_counts) <- gsub("_GSM.*$", "", colnames(raw_counts))
colnames(raw_counts) <- gsub("_Staphylococcus_aureus_RNA-Seq_sorted.bam", "", colnames(raw_counts))
#Preview
head(raw_counts)
# --- 7.  Load and clean metadata
# -----------------------------
# --- 8.  Metadata should have two columns:
#   1. sample - matching the sample names in count data
#   2. state  - biological condition (e.g., acute or chronic)
s_a_meta <- read.csv("metadata.csv", header = TRUE, stringsAsFactors = FALSE)

# --- 9.  Ensure the correct column names
colnames(s_a_meta) <- c("sample", "state")

# --- 10. Keep only samples present in both metadata and count matrix
s_a_meta <- s_a_meta[s_a_meta$sample %in% colnames(raw_counts), ]

# --- 11. Set sample names as row names (required for DESeq2)
rownames(s_a_meta) <- s_a_meta$sample

# --- 12.  Convert condition (state) to factor and set reference level
# Reference level determines the baseline for log2FoldChange interpretation
s_a_meta$state <- factor(s_a_meta$state, levels = c("chronic", "acute"))

# -----------------------------
# --- 13. Verify alignment of metadata and counts
# -----------------------------
# This ensures that sample order matches between metadata and counts.
stopifnot(all(colnames(raw_counts) == rownames(s_a_meta)))

# -----------------------------
# --- 14. Create DESeq2 dataset object
# -----------------------------
dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                              colData = s_a_meta,
# Preview                              design = ~ state)
dds
dds$sample
dds$state

# -----------------------------
# --- 15. Filter out lowly expressed genes
# -----------------------------
# Removes genes with total counts <= 10 across all samples
dds <- dds[rowSums(counts(dds)) > 10, ]


# -----------------------------
# --- 16. Run DESeq2 analysis
# -----------------------------
# Performs normalization, dispersion estimation, and differential testing
dds <- DESeq(dds)


# --- 17. Extract Results (Acute vs Chronic) ---
final_res <- results(dds)
head(final_res)

# --- 18. Distribution of p-values
plot(density(x = na.omit(final_res$pvalue)))


# --- 19. Visualizing the differentially expressed genes
plot(x = final_res$log2FoldChange,
     y = -log10(final_res$padj),
     cex= 0.5,
     pch = 19,
     col = 'grey',
     ylim = c(0,15),
     ylab = 'Adjusted P-Value',
     xlab = 'Log2 Fold Change')
abline(v = c(-2, 2), h = -log10(0.05), lwd = 0.5, lty = 2)

# --- 20. Subsetting upregulated genes for the visualization
upregulated <- subset(final_res, padj < 0.05 & log2FoldChange > 2)
points(upregulated$log2FoldChange,
       y = -log10(upregulated$padj),
       cex = 0.35,
       pch = 19,
       col = 'salmon')
# --- 21. Subsetting downregulated genes for the visualization
downregulated <- subset(final_res, padj < 0.05 & log2FoldChange < -2)
points(downregulated$log2FoldChange,
       y = -log10(downregulated$padj),
       cex = 0.35,
       pch = 19,
       col = 'lightblue')
mtext('DESeq Volcano Plot')

# --- 22. Heatmap of DEGs ---
if (nrow(upregulated) > 0 | nrow(downregulated) > 0) {
  degs <- rbind(
    raw_counts[rownames(upregulated), ],
    raw_counts[rownames(downregulated), ]
  )
# Combine both sets for visualization
  degs <- rbind(raw_counts[rownames(upregulated), ],
                raw_counts[rownames(downregulated), ])
  
  
  pheatmap(degs,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           show_rownames = FALSE,
           scale = "row",
           show_colnames = TRUE,
           annotation_col = s_a_meta,
           main = "Differentially Expressed Genes")

}
rownames(upregulated)
rownames(downregulated)

# --- 23. Export Results Safely
if (!dir.exists("results")) dir.create("results")
write.csv(upregulated, 'results/upregulated.csv')
write.csv(downregulated, 'results/downregulated.csv')
write.csv(raw_counts, 'results/raw_counts.csv')

# --- 24. Summary Counts

up_count <- nrow(upregulated)
down_count <- nrow(downregulated)
cat("Upregulated genes:", up_count, "\n")
cat("Downregulated genes:", down_count, "\n")

#Functional Enrichment Analysis
# Visit https://bioinformatics.sdstate.edu/go/

#############################################
# END OF SCRIPT
#############################################


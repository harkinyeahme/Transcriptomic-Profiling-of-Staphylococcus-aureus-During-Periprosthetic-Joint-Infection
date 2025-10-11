Transcriptomic Profiling of Staphylococcus aureus During Acute vs Chronic Phases of Periprosthetic Joint Infection (PJI)
# 🧬 Transcriptomic Profiling of *Staphylococcus aureus* During Chronic and Acute Infection

This repository contains the workflow, scripts, and analysis outputs for the **transcriptomic profiling of *Staphylococcus aureus*** during chronic and acute infection conditions. The aim of this study was to identify differentially expressed genes (DEGs) that may be associated with persistence or virulence during infection.

---

## 📊 Project Overview

**Project Title:** Transcriptomic Profiling of *Staphylococcus aureus* During Chronic and Acute Infection  
**BioProject Accession:** [PRJNA867318](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA867318)  
**SRA Accessions:** SRR20959676, SRR20959677, SRR20959678, SRR20959679, SRR20959680, SRR20959681  
**Organism:** *Staphylococcus aureus* subsp. *aureus* NCTC 8325  
**Data Source:** NCBI Sequence Read Archive (SRA)  
**Objective:** To compare the gene expression profiles of *S. aureus* during chronic and acute infections and identify significant transcriptomic signatures linked to infection states.

---

## 🧪 Methodology

### 1. **Data Retrieval**
Raw RNA-Seq datasets were downloaded from the NCBI SRA database using the listed accession numbers. Each dataset represents biological replicates of *S. aureus* infection models under acute and chronic conditions.

### 2. **Quality Control and Preprocessing**
Tools used:
- **FastQC** – Initial read quality assessment  
- **Fastp** – Adapter trimming, low-quality base filtering, and read quality enhancement  
- **MultiQC** – Aggregated summary report of QC results  

### 3. **Reference Genome and Alignment**
- Reference genome: *Staphylococcus aureus* NCTC 8325  
- Reads were aligned to the reference genome using **STAR** aligner.  
- Indexing was performed using the reference FASTA and GFF annotation files.

### 4. **Read Quantification**
- Gene-level read counts were generated using **FeatureCounts**.  
- Output: `results/counts.txt` containing raw counts for each gene across samples.

### 5. **Differential Expression Analysis**
- Conducted in **R** using the **DESeq2** package.  
- Genes were considered significantly differentially expressed if:  
  - **Adjusted p-value (FDR) < 0.05**  
  - **|log₂ Fold Change| ≥ 2**  
- Normalization and variance stabilization were performed prior to visualization.

### 6. **Functional Enrichment**
- Differentially expressed genes were analyzed for Gene Ontology (GO) and pathway enrichment using **ShinyGO**.

---

## 🧠 Results Summary

After quality filtering, alignment, and differential expression testing:

- Only **one gene**, **SAOUHSC_02191**, was found to be **significantly upregulated** between the chronic and acute infection conditions.  
- No genes were found to be significantly downregulated under the defined threshold (**|log₂FC| ≥ 2, adj. p < 0.05**).  
- Functional enrichment analysis using **ShinyGO** showed **no associated pathways** for SAOUHSC_02191 within current annotation databases.

### 🔍 Interpretation
- The expression profile revealed limited differential activity, possibly due to the **restricted sample size and limited dataset availability** in the SRA project.  
- The single upregulated gene, **SAOUHSC_02191**, encodes a **hypothetical protein** with no known biological function or pathway association.  
- This suggests potential involvement in unexplored mechanisms of bacterial adaptation or persistence.  
- Additional data or experimental validation (e.g., qPCR) may be required to confirm its biological relevance.

---

##  Tools and Environment

| Tool / Package | Version | Purpose |
|----------------|----------|----------|
| FastQC | v0.12+ | Raw read quality assessment |
| Fastp | v0.23+ | Adapter trimming and quality filtering |
| MultiQC | v1.15+ | Aggregated QC visualization |
| STAR | v2.7+ | Read alignment to reference genome |
| FeatureCounts | v2.0+ | Gene quantification |
| R / DESeq2 | R v4.3+, DESeq2 v1.42+ | Differential expression analysis |
| ShinyGO | Online | GO and pathway enrichment |

---

## 🧾 Key Output Files

| File | Description |
|------|--------------|
| `results/upregulated.csv` | Upregulated genes |
| `results/raw_counts.csv` | raw featurecounts of the genes |
| `multiqc_report.html` | Combined QC report (FastQC + Fastp) |
| `Aligned.sortedByCoord.out.bam` | STAR-aligned reads |
| `count.txt` | Gene count matrix (FeatureCounts output) |
| `DESeq2_results.csv` | Differential expression results |

---

## 🧬 Biological Insight

| Gene ID | Regulation | Description | Pathway Involvement |
|----------|-------------|--------------|----------------------|
| SAOUHSC_02191 | Upregulated | Hypothetical protein | None detected |

---

## ⚙️ Scripts Used

All analytical scripts are available in the `/scripts` directory:

- `01_download_data.sh` - Downloads the raw reads
- `02_qc_raw.sh` – Performs FastQC, and MultiQC for raw reads
- `03_trim.sh` - Performs trimming of raw data and generates post trim qc report
- `04_mapping_script..sh` – Builds genome index and aligns reads with STAR  
- `05_featureCounts.sh` – Quantifies aligned reads per gene  
- `06_DESeq2_script.R` – Runs DESeq2 and generates DEG results  

---

##  References

1. Love, M. I., Huber, W., & Anders, S. (2014). *Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.* Genome Biology, 15(12), 550.  
2. Dobin, A., et al. (2013). *STAR: ultrafast universal RNA-seq aligner.* Bioinformatics, 29(1), 15–21.  
3. Chen, S., et al. (2018). *fastp: an ultra-fast all-in-one FASTQ preprocessor.* Bioinformatics, 34(17), i884–i890.  

---


📧 [akinyemiolanrewaju@gmail.com](mailto:akinyemiolanrewaju@gmail.com)  
🌐 [LinkedIn](https://www.linkedin.com/in/akinyemi-olanrewaju)  

---

## 📜

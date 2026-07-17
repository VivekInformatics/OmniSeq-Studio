# 🧬 OmniSeq Studio

OmniSeq Studio is a robust, end-to-end R Shiny application designed for automated RNA-seq transcriptomics analysis. This tool bridges the gap between raw sequencing data and biological interpretation, featuring a hardened analytical engine designed to handle diverse datasets from small test samples to large-scale experiments.

## 🚀 Key Features

* **Hardened Analysis Pipeline:** Incorporates statistical fallbacks for dispersion estimation, ensuring successful execution even with minimal replicates or low variance.
* **Automated Data Sanitization:** Automatically removes constant/zero-variance genes and handles missing values, preventing common pipeline crashes during normalization and PCA.
* **Comprehensive Diagnostics:** Includes 18 integrated visualizations for end-to-end quality control, differential expression analysis, and functional enrichment.
* **Flexible Normalization:** Automatically handles standard DESeq2 workflows, including automated variance-stabilized transformation (VST) and log2-fallback normalization.
* **Interactive Exploration:** Built-in DataTables for gene list filtering and interactive `plotly` charts for sample-level diagnostics.
  
<img width="1917" height="955" alt="Screenshot 2026-07-17 221323" src="https://github.com/user-attachments/assets/5382f7c5-e9e3-4e00-bc4a-6969426ad830" />


## 🛠️ Technical Workflow

OmniSeq Studio follows a standardized bioinformatics workflow:
1.  **Data Ingestion & Cleaning:** Validation of counts and metadata, automatic sample matching, and variance-based gene filtering.
2.  **Statistical Modeling:** Negative Binomial modeling via `DESeq2` with automated fallback logic for dispersion estimation.
3.  **Visualization:** Automated generation of QC (PCA, Correlation, Density), DEG (MA, Volcano), and Enrichment (GO/KEGG) plots.

## 📋 Installation

To run this application, ensure you have R (v4.0+) installed, then run the following in your R console:

```r
# Install CRAN dependencies
install.packages(c("shiny", "ggplot2", "plotly", "DT", "pheatmap", "enrichR"))

# Install Bioconductor dependencies
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2"))

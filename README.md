# Omni Seq 🔬

> **Complete Transcriptomics Analysis Suite**

Omni Seq is a powerful, open-source R Shiny application built for comprehensive, end-to-end RNA-seq and transcriptomics data analysis. Designed with a sleek, modern dark-mode dashboard interface, it allows researchers to process data, perform differential expression analysis, uncover co-expression networks using WGCNA, and run functional enrichment analysis without writing code.

---

## ✨ Key Features

* **📁 Multi-Format Data Input:** Upload raw count matrices and sample metadata in CSV, TSV, TXT, or Excel (`.xlsx`) formats with a built-in 100 MB upload limit.


* **⚙️ Robust Data Preprocessing & Quality Control:**
* Automated data quality reporting and sample distribution boxplots.


* Flexible gene filtering based on low counts and zero-variance removal.


* Multiple normalization options including **TMM (edgeR)**, **RLE (DESeq2)**, **CPM**, and **Log2CPM**.


* Interactive PCA plots and sample correlation heatmaps (pre- and post-normalization).




* **🎯 Differential Gene Expression (DEG) Analysis:**
* Powered by `edgeR` GLM pipeline.


* Interactive **Volcano plots** and **MA plots** via Plotly.


* Top DEGs hierarchical clustering heatmaps and downloadable result tables.




* **🕸️ Weighted Gene Co-Expression Network Analysis (WGCNA):**
* Automated soft-threshold power selection and scale-free topology fit plots.


* Dynamic tree cutting for module detection and close-module merging.


* **Module-Trait relationship heatmaps** and module eigengene correlations.


* Identification and extraction of **hub genes** per co-expression module.




* **🧬 Functional Enrichment Analysis:**
* Automated Gene Ontology (GO) Biological Process and KEGG pathway enrichment using `clusterProfiler`.


* Interactive dot plots and downloadable enrichment tables for both human and mouse datasets.




* **📦 Complete Bundle Export:** Export individual result tables or download a complete compressed (`.zip`) bundle containing all processed files and tabular outputs.



---

## 🛠️ Installation & Dependencies

Omni Seq requires **R** (version 4.0 or higher) along with several CRAN and Bioconductor packages. The application script automatically checks for missing dependencies and installs them upon startup, but you can also install them manually:

```r
# Install CRAN packages
pkgs <- c("shiny", "shinydashboard", "DT", "plotly", "ggplot2", 
          "pheatmap", "dplyr", "tidyr", "readxl", "stringr", "reshape2",
          "WGCNA", "flashClust", "zip")

for(pkg in pkgs){
  if(!require(pkg, character.only = TRUE)) install.packages(pkg)
}

# Install Bioconductor packages
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

bioc_pkgs <- c("edgeR", "limma", "clusterProfiler", "org.Hs.eg.db", 
               "org.Mm.eg.db", "enrichplot", "WGCNA")

for(pkg in bioc_pkgs){
  if(!require(pkg, character.only = TRUE)) BiocManager::install(pkg)
}
```[cite: 1]

---

## 🚀 Running the Application

1. Clone this repository or download the main script file (`app.R`).
2. Open the file in **RStudio**.
3. Run the application using `shiny::runApp()` or click **Run App** in RStudio:

```r
shiny::shinyApp(ui = ui, server = server)
```[cite: 1]

---

## 📊 Workflow Guide

1. **Data Upload Tab:** Upload your raw count matrix (genes in rows, samples in columns) and sample metadata file[cite: 1]. Select your condition column, control group, and treatment group[cite: 1].
2. **Data Prep Tab:** Review the quality metrics report, inspect sample boxplots and PCA plots, configure your filtering/normalization thresholds, and click **Apply Preprocessing**[cite: 1].
3. **DEG Analysis Tab:** Click **Run DEG Analysis** to calculate fold changes and p-values. Explore interactive volcano/MA plots and download your differential expression tables[cite: 1].
4. **WGCNA Tab:** Set your soft threshold power, minimum module size, and gene subset limits, then run network analysis to explore modules, trait correlations, and hub genes[cite: 1].
5. **Enrichment Tab:** Map significant genes to Entrez IDs and generate GO and KEGG functional enrichment dot plots[cite: 1].
6. **Downloads Tab:** Download the comprehensive results bundle containing all data matrices, module assignments, and pathway tables in a single `.zip` file[cite: 1].

```

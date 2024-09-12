# RNA-Seq Transcriptome Profiling Analysis Project

## Project Overview

In this project, I performed a comprehensive RNA-Seq transcriptome profiling analysis to investigate the effects of Dexamethasone treatment on various cell lines. The primary aim was to identify differentially expressed genes and gain insights into the gene expression changes induced by this glucocorticoid.

## Technical Details

### Software and Tools Used

- **R Programming Language**: The primary language used for data analysis and visualization.
- **Bioconductor Packages**: Essential for RNA-Seq data analysis.
- **DESeq2**: Used for differential expression analysis of count data.
- **clusterProfiler**: For Gene Set Enrichment Analysis (GSEA) and functional annotation.
- **enrichplot**: To visualize enrichment results.
- **ComplexHeatmap**: For creating heatmaps with hierarchical clustering.

### Data Processing

#### Data Import

- RNA-Seq count data was read from CSV files using `read.csv()`.
- Metadata containing sample information was also imported to match with the count data.

#### Data Transformation and Filtering

- Applied variance-stabilizing transformation (VST) to normalize the data using `vst()` function from DESeq2.
- Filtered out low-count reads to retain only those with sufficient expression levels.

### Differential Expression Analysis

- **DESeq2**: Used to perform differential expression analysis with the `DESeq()` function.
- **Result Extraction**: Results were extracted with an adjusted p-value cutoff of 0.05. Significant results were annotated with gene symbols using `org.Hs.eg.db` package.

### Visualization

- **Volcano Plot**: Generated to visualize the relationship between log fold changes and -log10 p-values of the genes. Significant genes were highlighted with fold change and adjusted p-value thresholds.
- **PCA Plot**: Principal Component Analysis (PCA) was performed to visualize sample clustering based on gene expression.
- **Heatmap**: Created to visualize the expression patterns of the top 30 differentially expressed genes across samples.
- **MA Plot**: Visualized the log fold changes against mean expression levels to identify genes with strong differential expression.

### Gene Set Enrichment Analysis (GSEA)

- Performed to identify biological processes and pathways significantly enriched in the gene expression data.
- Visualized results using dot plots, network plots, ridge plots, and GSEA plots.

### Additional Analysis

- **Pairwise Term Similarity**: Analyzed term similarity to explore relationships between different gene sets.
- **Publication Metrics Plot**: Examined trends and metrics related to top terms over a range of years.

## Summary

This project demonstrates a comprehensive approach to RNA-Seq data analysis, from raw count data processing and differential expression analysis to advanced functional enrichment and visualization techniques. By integrating various R packages and statistical methods, the analysis provided valuable insights into the impact of Dexamethasone on gene expression and identified key genes and pathways associated with the treatment.

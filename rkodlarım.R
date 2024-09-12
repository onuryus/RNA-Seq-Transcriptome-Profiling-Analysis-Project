# Installation of necessary libraries
# Installation is done only once, so we check if the packages are already installed.
# install.packages("ggplot2")

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install()

# Installing additional necessary Bioconductor packages for RNA-Seq analysis
# BiocManager::install(c("GenomicFeatures"))
# BiocManager::install("DESeq2")


# Load the count data (gene expression counts for each sample)
cnt <- read.csv("counts.csv")
str(cnt) # Check the structure of the count data


# Load the metadata (sample information such as treatment and condition)
met <- read.csv("metadata2.csv")
str(met) # Check the structure of the metadata



# Ensuring that the column names in the count data match the row names in the metadata
all(colnames(cnt) %in% rownames(met))



# Checking the order of row names and column names to ensure they align correctly
all(colnames(cnt) == rownames(met))



# Loading the DESeq2 library, which is used for differential gene expression analysis
library(DESeq2)


# Creating a DESeq2 dataset object from the count matrix and metadata
# 'design = ~dexamethasone' indicates that we want to analyze the effect of Dexamethasone treatment
dds <- DESeqDataSetFromMatrix(countData = cnt, 
                              colData = met,
                              design = ~dexamethasone)


dds # Display the DESeqDataSet object


# Optional step: Removing genes with low counts (filtering out genes with fewer than 10 reads across all samples)
# This step is done because genes with very low counts are often not expressed in the tissue of interest
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds  # Display the filtered dataset


# Setting the reference level for Dexamethasone treatment
# The 'untreated' group is set as the reference to compare against the treated samples
dds$dexamethasone <- relevel(dds$dexamethasone, ref = "untreated")

# Performing differential expression analysis
deg <- DESeq

# Extracting the results of the analysis
res <- results(deg)

# Writing the results to a CSV file for further examination
write.csv(res,"test_udemy.csv")


# Summarizing the results (overview of significant genes, p-values, etc.)
summary(res)


# Filtering results with an adjusted p-value (alpha) threshold of 0.05
res0.05 <-results(deg, alpha = 0.05)

# Displaying a summary of the filtered results (genes with an adjusted p-value < 0.05)
summary(res0.05)


# Installing the human genome annotation package for mapping gene IDs to symbols (only run once)
# BiocManager::install("org.Hs.eg.db")

# Loading the org.Hs.eg.db package, which contains human genome annotation data
library(org.Hs.eg.db)


# Converting the results to a data frame for easier manipulation
res0.05.df <- as.data.frame(res0.05)
str(res0.05.df)# Checking the structure of the data frame


# Mapping ENSEMBL gene IDs to gene symbols using the human genome database
# 'ENSEMBL' represents the keytype of the IDs in the results, and we map them to 'SYMBOL'
res0.05.df$Symbol <- mapIds(org.Hs.eg.db, rownames(res0.05.df), keytype = "ENSEMBL", column = "SYMBOL" )


# Checking the structure of the data frame after adding gene symbols
str(res0.05.df)


# Writing the final results with gene symbols to a CSV file for further analysis or reporting
write.csv(res0.05.df, "final_test_udemy.csv")


# Building a PCA plot to visualize sample clustering based on gene expression
# 'vst' stands for variance stabilizing transformation, which normalizes the data
vsd <- vst(deg, blind= FALSE)

# Plotting the PCA using the 'dexamethasone' group to distinguish between treatments
plotPCA(vsd, intgroup= "dexamethasone")

# Displaying the size factors used for normalization, which adjusts for differences in sequencing depth
sizeFactors(deg)


# Estimating and visualizing the dispersion of gene expression
# Dispersion measures how much the expression of genes varies between samples
plotDispEsts(deg)


# Loading necessary libraries for further data manipulation and visualization
# dplyr for data manipulation and ggplot2 for plotting
# install.packages("dplyr")
# install.packages("ggplot2")
library(dplyr)
library(ggplot2)



# Building an MA plot, which visualizes log fold changes (y-axis) vs. mean expression (x-axis)
# It helps to identify genes with strong differential expression
plotMA(res0.05)


# Identifying the top 30 genes with the most significant p-adjusted values (smallest padj)
best_genes <- res0.05.df %>%
  arrange(padj)%>%   # Sorting the results by adjusted p-value
  head(30)           # Selecting the top 30 genes

# Displaying the top 30 genes
best_genes

# Writing the best 30 genes to a CSV file for further exploration
write.csv(best_genes,"best_genes.csv")


# Generation of a Volcano plot
# A volcano plot visualizes log fold change (x-axis) against the -log10 adjusted p-value (y-axis)
# It helps to identify significantly differentially expressed genes
vol <- res0.05.df %>%
  filter(!is.na(padj))       # Remove rows with NA values in the padj column



# Creating the volcano plot
# Points are colored based on whether they meet significance (padj < 0.05) and log2 fold change thresholds
ggplot(vol, aes(x=log2FoldChange, y =-log10(padj),color = padj < 0.05 & abs(log2FoldChange)>1)) + 
  geom_point() +  # Scatter plot of log fold change vs -log10 adjusted p-value
  geom_text(data = best_genes, aes(label = Symbol), hjust = -0.2, vjust=0.5)         # Label top genes


# Building a heatmap for the top 30 differentially expressed genes
# The ComplexHeatmap package is used for creating complex heatmaps with hierarchical clustering
#install.packages("rjson")
#BiocManager::install("ComplexHeatmap")

library(ComplexHeatmap)


# Selecting the top 30 genes with the lowest adjusted p-values (most significant)
top_genes <- res0.05.df %>%
  arrange(padj)%>%
  head(30)


# Extracting the normalized counts for the top 30 genes
mat <- counts(deg, normalized=T)[rownames(top_genes),]
head(mat, 5)  # Displaying the first 5 rows of the matrix


# Scaling the gene expression values for each gene (z-score transformation)
# This step helps in visualizing relative expression across samples
mat.z <- t(apply(mat,1,scale))
head(mat.z, 5) # Displaying the first 5 rows of the scaled matrix



# Renaming the columns with sample names from the metadata

colnames(mat.z) <- rownames(met)
head(mat.z, 5)



# Creating the heatmap
# Rows (genes) and columns (samples) are clustered based on expression patterns
Heatmap(mat.z, cluster_rows = T, cluster_columns= T, column_labels = colnames(mat.z), row_labels = top_genes$Symbol )




# Install necessary packages for gene set enrichment analysis and visualization (run only once)


#BiocManager::install("clusterProfiler", version = "3.19")
#BiocManager::install("pathview")
#BiocManager::install("enrichplot")

# Load required libraries for gene set enrichment analysis and visualization

library(clusterProfiler)
library(enrichplot)


# Prepare the gene list for enrichment analysis
# 'best_genes$log2FoldChange' provides the log fold changes for the top genes
# Gene names are assigned as names of the list
original_gene_list <- best_genes$log2FoldChange
names(original_gene_list) <- rownames(best_genes)



# Sorting the gene list in descending order of fold changes
original_gene_list = sort(original_gene_list, decreasing = TRUE)


# Performing Gene Set Enrichment Analysis (GSEA)
# gseGO function performs enrichment analysis for gene ontology terms
gse <- gseGO(geneList = original_gene_list, 
             ont = "ALL",                  # Analyze all gene ontology categories
             keyType = "ENSEMBL",          # Gene IDs are in ENSEMBL format
             minGSSize = 3,                # Minimum gene set size
             maxGSSize = 800,              # Maximum gene set size
             pvalueCutoff = 0.05,          # P-value cutoff for significance
             verbose = TRUE,               # Print progress messages
             OrgDb = org.Hs.eg.db,         # Annotation database for human genes
             pAdjustMethod = "none")       # No adjustment for p-values



# Load the DOSE package for additional enrichment analysis tools
require(DOSE)


# Creating a dot plot of the GSEA results
# Shows the top 10 categories, splitting the plot by the sign of enrichment
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)


# Performing pairwise term similarity analysis
# This helps in visualizing the relationships between different gene sets
gsee <- pairwise_termsim(gse)
emapplot(gsee, showCategory = 10)



# Creating a network plot of the enrichment results
# The size of each node represents the p-value or gene number
cnetplot(gse, categorySize="pvalue", foldChange=original_gene_list, showCategory = 3)


# Creating a ridge plot to visualize the distribution of enrichment scores
# Customizing the y-axis label size
ridgeplot(gse) + 
  labs(x = "enrichment distribution") +
  theme(axis.text.y = element_text(size = 7))  # Y eksenindeki etiketlerin boyutunu küçült



# Creating a GSEA plot for a specific gene set
# 'geneSetID' specifies which gene set to plot, and 'by' determines the plot type
gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 1)


# Plotting publication metrics for the top terms over the years
# 'terms' specifies which terms to plot, and '2010:2018' specifies the range of years
terms <- gse$Description[1:3]
pmcplot(terms, 2010:2018, proportion=FALSE)





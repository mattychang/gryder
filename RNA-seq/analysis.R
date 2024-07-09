### Goal: generate heatmaps to visualize RNA-seq data
### Pre-requisites: TPM matrix, gene list

library(pheatmap)
library(dendsort)

setwd("/Volumes/rc/SOM_GENE_BEG33/RNA_seq/hg38/projects/RMS_IHK/Practice_MSC/ExpMatrices")
project.folder = "/Volumes/SOM_GENE_BEG33/RNA_seq/hg38/projects/RMS_IHK/Practice_MSC"
sample.set = "RMS_IHK_RNA-seq_022924"



### 4. Heatmaps
### Note: First heatmap focuses on the relative expression of genes using log2-transformed and row-scaled values. It is useful for visualizing the relative expression levels of genes.
### Note: Second heatmap focuses on the clustering of genes and samples based on their expression profiles using raw values, row scaling, and hierarchical clustering. It is useful for identifying co-expressed genes and similar samples.

## Read the TPM matrix
EXP.coding.matrix = read.table(file.choose(), header = T)

## Define the gene list
gene.list = read.table(file = "/Volumes/rc/SOM_GENE_BEG33/RNA_seq/hg38/projects/Genesets/IO_Custom_Genesets/A485_and_dCBP1_down.txt", sep="\t", header=F)

## Set the cutoff for minimal expression
cutoff.expression.min = 10

## Generating the filtered matrix
EXP.expressed.matrix = EXP.coding.matrix                                                                      # create a copy of the expression matrix for filtering
EXP.expressed.matrix$maxTPM = apply(EXP.expressed.matrix[, 2:(ncol(EXP.expressed.matrix) - 1)], 1, FUN=max)   # calculate the maximum TPM value for each gene across all samples

## Filtering process
EXP.expressed.matrix = subset(EXP.expressed.matrix, EXP.expressed.matrix$maxTPM > cutoff.expression.min)      # subset the matrix to only include genes with a maximum TPM value above the threshold
EXP.expressed.matrix.TFs = subset(EXP.expressed.matrix, EXP.expressed.matrix$gene_id %in% gene.list$V1)       # subset the matrix to only include genes in the gene.list

## Generate the heatmap (log2 transformation)
## Note: drop the first (gene_id) and last (maxTPM) columns
pheatmap(log2(as.matrix(EXP.expressed.matrix.TFs[, 2:(ncol(EXP.expressed.matrix.TFs) - 1)]) + 1), 
         scale = 'row',
         main = paste(sample.set, " log2 TPM heatmap", sep = ""))                                              # each row (gene) is standardized (mean = 0, std = 1)

## Dendrogram
## Perform hierarchical clustering on samples and plot the dendrogram
## Sort for better visualization and replot the dendrogram
mat_cluster_cols <- hclust(dist(t(EXP.expressed.matrix[, 2:(ncol(EXP.expressed.matrix) - 1)])))
sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
mat_cluster_cols <- sort_hclust(mat_cluster_cols)
plot(mat_cluster_cols, 
     main = "Sorted Dendrogram", 
     xlab = "", 
     ylab = "", 
     sub = "", 
     axes = FALSE)
    
## Perform hierarchial clustering on genes and plot heatmap
mat_cluster_rows <- sort_hclust(hclust(dist(EXP.expressed.matrix[, 2:(ncol(EXP.expressed.matrix) - 1)])))
pheatmap(EXP.expressed.matrix.TFs[, 2:(ncol(EXP.expressed.matrix.TFs) - 1)], 
         scale = 'row', 
         cluster_cols = mat_cluster_cols, 
         cluster_rows = TRUE, 
         main = paste(sample.set, " TPM heatmap with hierarchial clustering", sep = ""))




######################### end of heatmap #########################

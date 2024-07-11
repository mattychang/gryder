################################
### RNA Matrix Analysis ########
## pairs with buildTPM_Matrix ##
## gryderlab.com 2021  #########
##################################################
###  build transcripts per million data matrix ###
###  then perform GSEA, heatmap comparisons,   ###
###  PCA plots, ranked scatter plots, and more ###
##################################################
###  by Berkley Gryder, gryderart@gmail.com    ###
###  edited by Diana 2022 for Mac              ###
##################################################

#build expression matrix from ChIP_seq/RNA_DATA/
#in buildTPM_Matrix.R file




### 4. Heatmap
### Goal: generate heatmaps to visualize RNA-seq data
### Pre-requisites: TPM matrix, gene list
### Note: First heatmap focuses on the relative expression of genes using log2-transformed and row-scaled values. It is useful for visualizing the relative expression levels of genes.
### Note: Second heatmap focuses on the clustering of genes and samples based on their expression profiles using raw values, row scaling, and hierarchical clustering. It is useful for identifying co-expressed genes and similar samples.

## Setup
library(plyr); library(dplyr)
library(pheatmap)
library(dendsort)

setwd("/Volumes/rc/SOM_GENE_BEG33/RNA_seq/hg38/projects/RMS_IHK/Practice_MSC/ExpMatrices")
project.folder = "/Volumes/SOM_GENE_BEG33/RNA_seq/hg38/projects/RMS_IHK/Practice_MSC"
sample.set = "RMS_IHK_RNA-seq_022924"

## Read the TPM matrix
EXP.coding.matrix = read.table(file.choose(), header = T)

## Drop the last 4 columns (IHK45)
EXP.coding.matrix <- EXP.coding.matrix %>% select(-((ncol(EXP.coding.matrix) - 3):ncol(EXP.coding.matrix)))

## Define the gene list
gene.list = read.table(file = "/Volumes/rc/SOM_GENE_BEG33/RNA_seq/hg38/projects/RMS_IHK/Practice_MSC/MSC_TF_list_custom.txt", sep="\t", header=F)

## Set the cutoff threshold for minimal expression
cutoff.expression.min = 10

## Generating the filtered matrix
EXP.expressed.matrix = EXP.coding.matrix                                                                      # create a copy of the expression matrix for filtering
EXP.expressed.matrix$maxTPM = apply(EXP.expressed.matrix[, 2:(ncol(EXP.expressed.matrix) - 1)], 1, FUN=max)   # calculate the maximum TPM value for each gene across all samples

## Filtering steps
EXP.expressed.matrix = subset(EXP.expressed.matrix, EXP.expressed.matrix$maxTPM > cutoff.expression.min)      # subset the matrix to only include genes with a maximum TPM value above the threshold
EXP.expressed.matrix.TFs = subset(EXP.expressed.matrix, EXP.expressed.matrix$gene_id %in% gene.list$V1)       # subset the matrix to only include genes in the gene.list
rownames(EXP.expressed.matrix.TFs) <- EXP.expressed.matrix.TFs$gene_id                                        # match rownames with the gene_id

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
    





######################### come back to later #########################
###matrix of scatters and Pearson Corr.
library(ggplot2)
library(GGally)
GGally::ggpairs(log2((EXP.expressed.matrix.TFs[, 2:(ncol(EXP.expressed.matrix.TFs) - 1)]) + 1))

###scatter
library(LSD)
######################### come back to later #########################






### 5. PCA
### Goal: Perform PCA on the log2-transformed expression matrix to visualize sample clustering.
### Prerequisites: TPM matrix, sample list with subtype/drug Tx and timepoint info
### Note: subtype and timepoint info is optional. Also, ensure that samples are labeled properly

## Setup
library(tidyverse)
library(ggrepel)

## Convert the expression matrix to a numeric matrix
EXP.pca = as.matrix(EXP.coding.matrix[, 2:ncol(EXP.coding.matrix)])                               # create PCA matrix
rownames(EXP.pca) = EXP.coding.matrix$gene_id                                                     # add rownames

## Apply a log2 transformation to the matrix (adding 1 to avoid log(0))
EXP.pca.log2 = log2(EXP.pca + 1)

## Get the list of sample names (column names of the PCA matrix)
sample.name.list = c(colnames(EXP.pca))

## OPTIONAL FOR GROUPING PURPOSES
## Get sample subtype and timepoint information
sample_class = read.table(file.choose(), header = F)
sample.subtypes <- sample_class$V2[match(sample.name.list, sample_class$V1)]
timepoints <- sample_class$V3[match(sample.name.list, sample_class$V1)]
timepoints <- as.factor(timepoints)                                                               # convert timepoints to a factor (interpreting as continuous, when it is actually discrete)


# Create a dataframe with sample names and set the row names to sample names
EXP.coldata = data.frame(sample.name.list, sample.subtypes, timepoints)
rownames(EXP.coldata) = sample.name.list

# Perform PCA on the log-transformed data
pca = EXP.pca.log2 %>% t %>% prcomp

# Convert PCA results to a data frame and add sample names
EXP.d = pca$x %>% as.data.frame
EXP.d$sample.name.list = EXP.coldata$sample.name.list

# Join the PCA data frame with the sample metadata
EXP.dmeta = join(EXP.d, EXP.coldata, by = "sample.name.list")
EXP.dmeta$sample.name.list <- gsub("_RNA_022924_CWRU", "", EXP.dmeta$sample.name.list)            # remove "_RNA_022924_CWRU" from sample.name.list
EXP.dmeta = na.omit(EXP.dmeta)                                                                    # omit missing values

# Calculate the proportion of variability explained by each principal component
pcv = round((pca$sdev)^2 / sum(pca$sdev^2) * 100, 2)
      
# Graph the PCA plot
plot.pca = ggplot(EXP.dmeta, aes(PC2, PC3, colour = sample.subtypes, shape = timepoints)) + 
    geom_point() +
    xlab(label = paste0("PC2 (", pcv[2], "%)")) +
    ylab(label = paste0("PC3 (", pcv[3], "%)")) +
    theme_bw() +
    geom_label_repel(aes(label = sample.name.list), show.legend = FALSE) +
    theme(axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15)) +
    labs(title = "PCA of IHK44 and comparison drugs",
         subtitle = "Grouped by drug treatment and timepoint")

# Print the PCA plot
print(plot.pca)


######################### end of PCA #########################




      
### 5. GSEA ranklist maker
### Goal: To create rank lists for Gene Set Enrichment Analysis (GSEA) by calculating log2 fold changes (log2FC) for comparisons and saving the results.
### Pre-requisites: TPM matrix
### Notes: Must specify the log2FC for the samples of interest by hand

## Load in TPM expression matrix
setwd("/Volumes/SOM_GENE_BEG33/RNA_seq/hg38/")
project.folder = "/Volumes/SOM_GENE_BEG33/RNA_seq/hg38/projects/IHK_RMS/Practice_MSC"
EXP.coding.matrix = read.table(file.choose(), header = T)

## Create a directory for the GSEA ranklists and create the log2FC matrix
dir.create(file.path(project.folder, "GSEA_ranklist"))
EXP.GSEA = EXP.coding.matrix

## Calculate log2FC for each sample and add as new column
EXP.GSEA$RH4_DMSO_2h = log2(EXP.GSEA$RH4_DMSO_2h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_NT_2h_RNA_022924_CWRU + 1)
EXP.GSEA$RH4_DMSO_6h = log2(EXP.GSEA$RH4_DMSO_6h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_NT_6h_RNA_022924_CWRU + 1)
EXP.GSEA$RH4_A485_100nM_2h = log2(EXP.GSEA$RH4_A485_100nM_2h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_NT_2h_RNA_022924_CWRU + 1)
EXP.GSEA$RH4_A485_1uM_2h = log2(EXP.GSEA$RH4_A485_1uM_2h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_NT_2h_RNA_022924_CWRU + 1)
EXP.GSEA$RH4_A485_100nM_6h = log2(EXP.GSEA$RH4_A485_100nM_6h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_NT_6h_RNA_022924_CWRU + 1)
EXP.GSEA$RH4_A485_1uM_6h = log2(EXP.GSEA$RH4_A485_1uM_6h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_NT_6h_RNA_022924_CWRU + 1)
EXP.GSEA$RH4_JQAD_100nM_2h = log2(EXP.GSEA$RH4_JQAD_100nM_2h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_NT_2h_RNA_022924_CWRU + 1)
EXP.GSEA$RH4_JQAD_1uM_2h = log2(EXP.GSEA$RH4_JQAD_1uM_2h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_NT_2h_RNA_022924_CWRU + 1)
EXP.GSEA$RH4_JQAD_100nM_6h = log2(EXP.GSEA$RH4_JQAD_100nM_6h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_NT_6h_RNA_022924_CWRU + 1)
EXP.GSEA$RH4_JQAD_1uM_6h = log2(EXP.GSEA$RH4_JQAD_1uM_6h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_NT_6h_RNA_022924_CWRU + 1)
EXP.GSEA$RH4_dCBP_100nM_2h = log2(EXP.GSEA$RH4_dCBP_100nM_2h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_NT_2h_RNA_022924_CWRU + 1)
EXP.GSEA$RH4_dCBP_1uM_2h = log2(EXP.GSEA$RH4_dCBP_1uM_2h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_NT_2h_RNA_022924_CWRU + 1)
EXP.GSEA$RH4_dCBP_100nM_6h = log2(EXP.GSEA$RH4_dCBP_100nM_6h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_NT_6h_RNA_022924_CWRU + 1)
EXP.GSEA$RH4_dCBP_1uM_6h = log2(EXP.GSEA$RH4_dCBP_1uM_6h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_NT_6h_RNA_022924_CWRU + 1)
EXP.GSEA$RH4_QL_100nM_2h = log2(EXP.GSEA$RH4_QL_100nM_2h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_NT_2h_RNA_022924_CWRU + 1)
EXP.GSEA$RH4_QL_1uM_2h = log2(EXP.GSEA$RH4_QL_1uM_2h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_NT_2h_RNA_022924_CWRU + 1)
EXP.GSEA$RH4_QL_100nM_6h = log2(EXP.GSEA$RH4_QL_100nM_6h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_NT_6h_RNA_022924_CWRU + 1)
EXP.GSEA$RH4_QL_1uM_6h = log2(EXP.GSEA$RH4_QL_1uM_6h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_NT_6h_RNA_022924_CWRU + 1)
EXP.GSEA$RH4_LS_100nM_2h = log2(EXP.GSEA$RH4_LS_100nM_2h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_NT_2h_RNA_022924_CWRU + 1)
EXP.GSEA$RH4_LS_1uM_2h = log2(EXP.GSEA$RH4_LS_1uM_2h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_NT_2h_RNA_022924_CWRU + 1)
EXP.GSEA$RH4_LS_100nM_6h = log2(EXP.GSEA$RH4_LS_100nM_6h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_NT_6h_RNA_022924_CWRU + 1)
EXP.GSEA$RH4_LS_1uM_6h = log2(EXP.GSEA$RH4_LS_1uM_6h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_NT_6h_RNA_022924_CWRU + 1)
EXP.GSEA$RH4_IHK44_100nM_2h = log2(EXP.GSEA$RH4_IHK44_100nM_2h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_NT_2h_RNA_022924_CWRU + 1)
EXP.GSEA$RH4_IHK44_1uM_2h = log2(EXP.GSEA$RH4_IHK44_1uM_2h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_NT_2h_RNA_022924_CWRU + 1)
EXP.GSEA$RH4_IHK44_100nM_6h = log2(EXP.GSEA$RH4_IHK44_100nM_6h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_NT_6h_RNA_022924_CWRU + 1)
EXP.GSEA$RH4_IHK44_1uM_6h = log2(EXP.GSEA$RH4_IHK44_1uM_6h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_NT_6h_RNA_022924_CWRU + 1)
EXP.GSEA$RH4_IHK45_100nM_2h = log2(EXP.GSEA$RH4_IHK45_100nM_2h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_NT_2h_RNA_022924_CWRU + 1)
EXP.GSEA$RH4_IHK45_1uM_2h = log2(EXP.GSEA$RH4_IHK45_1uM_2h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_NT_2h_RNA_022924_CWRU + 1)
EXP.GSEA$RH4_IHK45_100nM_6h = log2(EXP.GSEA$RH4_IHK45_100nM_6h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_NT_6h_RNA_022924_CWRU + 1)
EXP.GSEA$RH4_IHK45_1uM_6h = log2(EXP.GSEA$RH4_IHK45_1uM_6h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_NT_6h_RNA_022924_CWRU + 1)
# there should be 63 columns


## Round the log2FC to 5 sig figs
EXP.GSEA[,(ncol(EXP.coding.matrix) + 1) : ncol(EXP.GSEA)] = round(EXP.GSEA[,(ncol(EXP.coding.matrix) + 1) : ncol(EXP.GSEA)], digits = 5)

## Loop through the new columns, create rank lists, and save them to files
for (i in (ncol(EXP.coding.matrix)+1):ncol(EXP.GSEA)) {
  
    # initialize rank list with gene names
    Ranklist <- data.frame(EXP.GSEA[, 4])
    
    # Add the log2FC values to the rank list
    Ranklist$DeltaTPM <- EXP.GSEA[, i]
    
    # Sort the rank list by DeltaTPM in descending order
    Ranklist = Ranklist[rev(order(Ranklist$DeltaTPM)),]
    
    # Get the name of the current comparison
    SampleName = colnames(EXP.GSEA)[i]
    
    # Get the current date
    mytime <- format(Sys.time(), "%b_%d_%Y")
    
    # Define the file path
    myfile <- file.path(project.folder, 
                        "GSEA_ranklist", 
                        paste0(SampleName,"_",mytime,".rnk"))
    
    # Save the rank list to a file
    write.table(Ranklist, 
                file = myfile, 
                sep = "\t", 
                row.names = FALSE, 
                col.names = FALSE,
                quote = FALSE, 
                append = FALSE)
}












#load GSEA library
#library(GSEA)
#pass ranklists into GSEA
#pass .gmt file (gene sets specific to the project)
#write out the output into a folder on the Google Drive
#write out cool PDF figures

### 6. Heatmap the data
## PREP DATA   
    GeneSet1.name = "RMS_proteins"
    #GeneSet1 = read.table(file = "K:/projects/ChIP_seq/RNA_DATA/RNA_projects/Genesets/Qlucore_format/IAD_MYC_SLAMseq.genelist.txt", sep="\t", header=F)
    GeneSet1 = as.data.frame(c("MYCN", "MYOD1", "MYOG", "SOX8", "PAX3", "FOXO1"))  #custom cut
    EXP.GSEA.GeneSet1 = subset(EXP.GSEA, EXP.GSEA$GeneID %in% GeneSet1[,1])
    EXP.GSEA.GeneSet1.log2 = EXP.GSEA.GeneSet1;# EXP.GSEA.GeneSet1.log2[,c(5,6,7,8)] = log2(EXP.GSEA.GeneSet1.log2[,c(5,6,7,8)]+1) #manual column selection 
   
    # GeneSet2.name = "Bromodomain proteins"
    # GeneSet2 = as.data.frame(c("BRD2", "BRD3", "BRD4","BRDT"))
    # EXP.GSEA.GeneSet2 = subset (EXP.GSEA, EXP.GSEA$GeneID %in% GeneSet2[,1])
    #     EXP.GSEA.all.GeneSet2 = subset (EXP.GSEA.all, EXP.GSEA.all$GeneID %in% GeneSet2[,1])
    
    sample.set = "RMS_all_NT"
    
    dir.create(file.path(project.folder,sample.set))
    myfile <- file.path(project.folder,sample.set, paste0(GeneSet1.name,".GSEA_matrix.txt"))
    write.table(EXP.GSEA.GeneSet1, file = myfile, sep = "\t", row.names = F, col.names = T, quote = FALSE, append = FALSE)
    myfile.log2 <- file.path(project.folder,sample.set, paste0(GeneSet1.name,".GSEA_matrix.log2.txt"))
    write.table(EXP.GSEA.GeneSet1.log2, file = myfile.log2, sep = "\t", row.names = F, col.names = T, quote = FALSE, append = FALSE)
    
    allfile <- file.path(project.folder,sample.set, paste0(sample.set,".allgenes.GSEA_matrix.txt"))
    
    write.table(EXP.GSEA, file = allfile, sep = "\t", row.names = F, col.names = T, quote = FALSE, append = FALSE)

    EXP.GSEA.GeneSet1 = gsub(EXP.GSEA.GeneSet1, pattern = "KCCF14_", replacement = "RH41_")
    
### PLOT DATA
#plot log2 scale, not delta values
    EXP.GSEA.GeneSet1$minTPM <- NULL
    EXP.GSEA.GeneSet1$maxTPM <- NULL
    
    EXP.GSEA.GeneSet1.plot = EXP.GSEA.GeneSet1[,c(1:34)]
    # EXP.GSEA.GeneSet1.plot = EXP.GSEA.GeneSet1[,c(1:4, 34:45)] #manual column selection  
    matrix.GeneSet1.plot = as.matrix(EXP.GSEA.GeneSet1.plot[,5:(ncol(EXP.GSEA.GeneSet1.plot))])
    rownames(matrix.GeneSet1.plot) = EXP.GSEA.GeneSet1.plot$GeneID
    matrix.GeneSet1.plot.log2 = log2(matrix.GeneSet1.plot+1)
    matrix.log.bounds = 11
    matrix.GeneSet1.plot.log2[matrix.GeneSet1.plot.log2>matrix.log.bounds] <- matrix.log.bounds; #matrix.GeneSet1.plot.log2[matrix.GeneSet1.plot<(-matrix.log.bounds)] <- -matrix.log.bounds
    #library(pheatmap)
    breaklist = seq(-matrix.log.bounds, matrix.log.bounds, by = 0.1)
    pheatmap(matrix.GeneSet1.plot.log2,scale = "none", cluster_rows = F, cluster_cols = T,main=paste(GeneSet1.name," log2 TPM heatmap",sep=""))#, color = colorRampPalette(c("blue", "white", "orange")))

#plot log 2 fold change values
    library(pheatmap)
    EXP.GSEA.GeneSet1.plot.log2fc = EXP.GSEA.GeneSet1[c(36:ncol(EXP.GSEA.GeneSet1))]
      #EXP.GSEA.GeneSet1.plot.log2fc = EXP.GSEA.GeneSet1[,(ncol(EXP.coding.matrix)+1):ncol(EXP.GSEA)] #manual column selection  
    matrix.GeneSet1.plot.log2fc = as.matrix(EXP.GSEA.GeneSet1.plot.log2fc)
    rownames(matrix.GeneSet1.plot.log2fc) = EXP.GSEA.GeneSet1$GeneID
      #matrix.log.bounds = 11
      #matrix.GeneSet1.plot.log2[matrix.GeneSet1.plot.log2>matrix.log.bounds] <- matrix.log.bounds; #matrix.GeneSet1.plot.log2[matrix.GeneSet1.plot<(-matrix.log.bounds)] <- -matrix.log.bounds
    breaklist = seq(-0.9, 0.9, by = 0.1)
    pheatmap(matrix.GeneSet1.plot.log2fc,scale = "column", cluster_rows = T, cluster_cols = T,breaks = breaklist[50] ,main=paste(sample.set," ",GeneSet1.name," log2 TPM heatmap",sep=""))#, color = colorRampPalette(c("blue", "white", "orange")))

    EXP.GSEA.GeneSet2.plot.log2fc = EXP.GSEA.GeneSet2[c(36:ncol(EXP.GSEA.GeneSet2))]
    matrix.GeneSet2.plot.log2fc = as.matrix(EXP.GSEA.GeneSet2.plot.log2fc)
    rownames(matrix.GeneSet2.plot.log2fc) = EXP.GSEA.GeneSet2$GeneID
    breaklist = seq(-0.9, 0.9, by = 0.1)
    pheatmap(matrix.GeneSet2.plot.log2fc,scale = "column", cluster_rows = T, cluster_cols = T,breaks = breaklist[50] ,main=paste(sample.set," ",GeneSet2.name," log2 TPM heatmap",sep=""))#, color = colorRampPalette(c("blue", "white", "orange")))
        EXP.GSEA.all.GeneSet2.plot.log2fc = EXP.GSEA.all.GeneSet2[c(54:76)]
        matrix.GeneSet2.all.plot.log2fc = as.matrix(EXP.GSEA.all.GeneSet2.plot.log2fc)
        rownames(matrix.GeneSet2.all.plot.log2fc) = EXP.GSEA.all.GeneSet2$GeneID
        breaklist = seq(-0.9, 0.9, by = 0.1)
        pheatmap(matrix.GeneSet2.all.plot.log2fc,scale = "none", cluster_rows = T, cluster_cols = T,breaks = breaklist[50] ,main=paste(sample.set," ",GeneSet2.name," log2 TPM heatmap",sep="")) +
          theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) #, color = colorRampPalette(c("blue", "white", "orange")))
        
    

## Prettier heatmaps
### with improved clustering
#sort it out a plot again. from http://slowkow.com/notes/heatmap-tutorial/

mat_cluster_cols <- hclust(dist(t(matrix.GeneSet1.plot.log2fc)))
plot(mat_cluster_cols, main = "Unsorted Dendrogram", xlab = "", sub = "")

library(dendsort)
sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
mat_cluster_cols <- sort_hclust(mat_cluster_cols)
plot(mat_cluster_cols, main = "Sorted Dendrogram", xlab = "", sub = "")
mat_cluster_rows <- sort_hclust(hclust(dist(matrix.GeneSet1.plot.log2fc)))

pheatmap(matrix.GeneSet1.plot.log2fc,scale='column',cluster_cols=F,cluster_rows=T,  main=paste(sample.set," ",GeneSet1.name," log2 TPM heatmap",sep=""))

### 7. Bar plots of a gene
#library(ggplot2)
#library(tidyr)
GeneSet1 = as.data.frame(c("KLK3", "GAPDH"))  #custom cut
EXP.GSEA.GeneSet1 = subset(EXP.GSEA, EXP.GSEA$GeneID %in% GeneSet1[,1])
EXP.GSEA.GeneSet1 = EXP.GSEA.GeneSet1[,c(1:4, 34:45)]
EXP.GSEA.GeneSet1.long = gather(EXP.GSEA.GeneSet1, Sample, TPM,5:ncol(EXP.GSEA.GeneSet1))
#EXP.GSEA.GeneSet1.long$Sample = with(EXP.GSEA.GeneSet1.long, reorder(Sample, TPM, median))
ggplot(EXP.GSEA.GeneSet1.long, aes(x = Sample, y = TPM, fill =GeneID))+geom_bar(stat = 'identity')+ coord_cartesian(ylim=c(-3,1))+
  theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))

ggplot(EXP.GSEA.GeneSet1.long %>% mutate(Sample = factor(Sample, levels = c("L2FC_LAMP_BG15n", "L2FC_LNCaP_BG15a", "L2FC_LNCaP_BG15n", "L2FC_LNCaP_XIP_BG15a", "L2FC_LNCaP_XIP_BG15n", "L2FC_X22RV1_BG15a", "L2FC_X22RV1_BG15n", "L2FC_LuCaP170_BG15_BU", "L2FC_LuCaP167_BG15n", "L2FC_LuCaP173_BG15_BU","L2FC_PC3_BG15a","L2FC_PC3_BG15n"))),aes(x = Sample, y = TPM))+
geom_bar(stat = 'identity')+facet_wrap(~GeneID, nrow = 3)+
theme_bw() + theme(axis.text.x = element_text(angle = 45)) 

### 8. Boxplots for gene lists
library(ggplot2)
library(tidyverse)
GeneSet1.name = "AR_pathway_down"
GeneSet1 = as.data.frame(c("KLK3"))  #custom cut => "AR","KLK2","TMPRSS2","FKBP5","NKX3-1","HPGD"

EXP.GSEA.GeneSet1 = subset(EXP.GSEA, EXP.GSEA$GeneID %in% GeneSet1[,1])
  colnames(EXP.GSEA.GeneSet1)
  EXP.GSEA.GeneSet1 <- EXP.GSEA.GeneSet1[, c(1,2,3,4,5,7,9,6,8,10)]
  EXP.GSEA.GeneSet1.long = gather(EXP.GSEA.GeneSet1, Sample, TPM, 5:ncol(EXP.GSEA.GeneSet1), factor_key = T) %>%
  group_by(Sample)

EXP.GSEA.L2FC.GeneSet1 = subset(EXP.GSEA.L2FC, EXP.GSEA$GeneID %in% GeneSet1[,1])
  colnames(EXP.GSEA.L2FC.GeneSet1)
  #EXP.GSEA.L2FC.GeneSet1 <- EXP.GSEA.L2FC.GeneSet1[, c(1,2,3,4,5,7,9,6,8,10)]
  EXP.GSEA.L2FC.GeneSet1.long = gather(EXP.GSEA.L2FC.GeneSet1, Sample, TPM, 5:ncol(EXP.GSEA.L2FC.GeneSet1), factor_key = T) %>%
    group_by(Sample)

ggplot(EXP.GSEA.GeneSet1.long, aes(x = Sample, y = TPM))+
  geom_bar(stat = 'identity')+
  ylab("L2FC(TPM)")
  facet_wrap(~GeneID, nrow = 3)+ theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust= 1))

ggplot(EXP.GSEA.L2FC.GeneSet1.long, aes(x = Sample, y = TPM))+
  geom_histogram(stat = 'identity', color ="black")+facet_wrap(~GeneID, nrow = 3) +
  ylab("L2FC(TPM)")+
  theme_bw()+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust= 1))

ggplot(subset(EXP.GSEA.GeneSet1[,c(36:ncol(EXP.GSEA.GeneSet1))], EXP.GSEA.GeneSet1$GeneID %in% c("AR","KLK2","KLK3","TMPRSS2","FKBP5","NKX3-1","HPGD") ), aes(x=Sample, y = TPM ,fill=GeneSet1)) + geom_violin() + geom_boxplot(outlier.shape = NA, width = 0.5) + theme_bw()
  # annotate("text",x=1.5,y=-1.7,label=paste("Welch t-test, )) +annotate("segment", x = 1.1, xend = 1.9, y = 1.8, yend = 1.8)
  #annotate("text",x=1.5,y=-1.7,label=paste("Welch t-test, p = ",ttestpVal)) +annotate("segment", x = 1.1, xend = 1.9, y = 1.8, yend = 1.8)


GeneSet2.name = "HOUSEKEEP"
GeneSet2 = read.table(file = "Q:/sxg1131/pRNA/RNA_projects/Genesets/Qlucore_format/HouseKeep_1819_min1.txt", sep="\t", header=F)
EXP.GSEA.GeneSet2 = subset(EXP.GSEA, EXP.GSEA$GeneID %in% GeneSet2[,1])
EXP.GSEA.GeneSet2$GeneSet = GeneSet2.name 

Geneset3.name = "AR"
Geneset3 = read.table(file = "Q:/sxg1131/pRNA/RNA_projects/Genesets/Qlucore_format/", sep="\t", header=F)
EXP.GSEA.GeneSet3 = subset(EXP.GSEA, EXP.GSEA$GeneID %in% GeneSet3[,1])
EXP.GSEA.GeneSet3$GeneSet = GeneSet3.name
EXP.Genesets = rbind(EXP.GSEA.GeneSet1,EXP.GSEA.GeneSet2, EXP.GSEA)

ttestpVal = round(t.test(EXP.GSEA.GeneSet1$L2FC_BG15n_rep1, EXP.GSEA.GeneSet2$L2FC_BG15n_rep1, paired = F)$p.value, 4)

ggplot(EXP.Genesets,aes(y=L2FC_BG15n_rep1,x=GeneSet,fill=GeneSet))+geom_violin()+geom_boxplot(outlier.shape = NA, width = 0.5)+theme_bw()+
  # annotate("text",x=1.5,y=-1.7,label=paste("Welch t-test, )) +annotate("segment", x = 1.1, xend = 1.9, y = 1.8, yend = 1.8)
  annotate("text",x=1.5,y=-1.7,label=paste("Welch t-test, p = ",ttestpVal)) +annotate("segment", x = 1.1, xend = 1.9, y = 1.8, yend = 1.8)

#########################################################################################################
##  Collect GSEA results.  Works with a limited set of gene sets all in the same GSEA output folder.   ##
#########################################################################################################
##  Berkley Gryder, 2018-2019 (gryderart@gmail.com): https://github.com/GryderArt/VisualizeRNAseq      ##
#########################################################################################################

### 1. Find all the folder names

setwd("Q:/RNA_seq/hg19/projects/GSEA_ranklist")
GSEA.folder = "ARIA_LNCaP"
samples <- list.dirs(path=GSEA.folder, full.names=F, recursive=F)  #DONT USE "." in GSEA Analysis Names, will cause incorrect parsing
#samples = samples[grep("1645135", samples, invert = F)] 
samples = samples[grep("LNCaPXIP", samples, invert = F)]
samples = samples[grep("6991|9187", samples, invert = T)]
### 2. Edit string from folder names to get report names (xls file name is the same as the sheet name, minus a few digits)
sample.df = read.table(text = samples, sep=".")
sample.df$V3 = as.character(sample.df$V3)
colnames(sample.df)=c("Model_Tx","GSEA_Style","Digits")
sample.df$neg.xls.path = paste(GSEA.folder,"/",paste(sample.df$Model_Tx,sample.df$GSEA_Style,sample.df$Digits,sep="."),"/gsea_report_for_na_neg_",sample.df$Digits,".tsv",sep="")
sample.df$pos.xls.path = paste(GSEA.folder,"/",paste(sample.df$Model_Tx,sample.df$GSEA_Style,sample.df$Digits,sep="."),"/gsea_report_for_na_pos_",sample.df$Digits,".tsv",sep="")

#subset
#sample.df = sample.df[grep("weight", sample.df$Model_Tx), ]

### 3. Import GSEA data 
#initiate dataframe
allsamples = t(data.frame(sample = c("NAME","ES","NES","NOM.p.val","FDR.q.val","FWER.p.val","Model_Tx")))
colnames(allsamples) = c("NAME","ES","NES","NOM.p.val","FDR.q.val","FWER.p.val","Model_Tx")
allsamples=allsamples[-1,]; df.empty = allsamples

lapply(sample.df$Model_Tx, function(x) {
  #read in fake xls file
  temp.df <- subset(sample.df,sample.df$Model_Tx %in% x)
  
  if(file.exists(temp.df$neg.xls.path) == 'TRUE'){
    df.neg = read.table(temp.df$neg.xls.path,sep="\t",header=T)
  }else{df.neg = df.empty}
  
  if(file.exists(temp.df$pos.xls.path) == 'TRUE'){
    df.pos = read.table(temp.df$pos.xls.path,sep="\t",header=T)
  }else{df.pos = df.empty}
  
  df = rbind(df.pos,df.neg)
  #make it smaller, add column for identification
  df = df[,c("NAME","ES","NES","NOM.p.val","FDR.q.val","FWER.p.val")]
  df$Model_Tx = x
  allsamples <<- rbind(allsamples,df)
})

#simplify names
allsamples$Model_Tx = gsub(allsamples$Model_Tx, pattern = "HALLMARK_", replacement = "")
allsamples$Model_Tx = gsub(allsamples$Model_Tx, pattern = "HALLMARKS_", replacement = "")
allsamples$Model_Tx = gsub(allsamples$Model_Tx, pattern = "L2FC_", replacement = "")
#allsamples$Model_Tx = gsub(allsamples$Model_Tx, pattern = "_01062022", replacement = "")
allsamples$NAME = gsub(allsamples$NAME, pattern = "HALLMARK_", replacement = "")

#Account for missing data
#remove non-numeric values
#allsamples = allsamples[-grep("HALLMARK_COAGULATION",allsamples$NAME),]
#reindex
row.names(allsamples) <- 1:nrow(allsamples)
#as.numeric(allsamples$NOM.p.val)+0.00001


#calculate some statistics
allsamples$NOM.p.val[is.na(allsamples$NOM.p.val)] <- 0

allsamples$NOM.p.val = as.numeric(allsamples$NOM.p.val) + 0.00001
allsamples$log10.NOM.p.val= log10(allsamples$NOM.p.val)
allsamples$NES = as.numeric(allsamples$NES)

#reorder by Enrichment Score (NES)
NAME.ranks = aggregate(NES ~ NAME, allsamples, mean);colnames(NAME.ranks)=c("NAME","rankmetric")
library(plyr)
allsamples <- join(allsamples, NAME.ranks, by = "NAME")
allsamples$NAME <- factor(allsamples$NAME, levels=allsamples[order(unique(allsamples$rankmetric,decreasing=F)),]$NAME)  
#allsamples$NAME = gsub(allsamples$NAME, pattern = "HALLMARK_", replacement = "")
 allsamples <- allsamples[-c(10)]
ES.min = min(abs(allsamples$ES))
NES.min = min(abs(allsamples$NES))

#allsamples$NAME <- factor(allsamples$NAME, levels = c("HALLMARK_APOPTOSIS","GRYDER_HOUSEKEEPING","TFS_NO_EPIMACHINES","GRYDER_RH4_SE_GENES","GRYDER_RH4_TOP_SE_TFS")) #order genesets

write.table(allsamples, paste(GSEA.folder,"/",GSEA.folder,".GSEA.allsamples.summary.txt",sep=""),col.names = T,row.names = F,quote = F,sep="\t")

### 4. Plots
# a. Bubble chart plots
library(tidyr)
library(ggplot2)

ggplot(allsamples, aes(x = allsamples$Model_Tx, y = allsamples$NAME, size = -log10.NOM.p.val,colour = NES)) + 
  geom_point() +scale_colour_gradient2(low="red",mid="white", high="green4") +
  ylab("GSEA Hallmark Pathways") + xlab("Sample")+
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(allsamples %>% mutate(Model_Tx = factor(Model_Tx, levels = c("LAMP_BG15n_", "LNCaP_BG15a_", "LNCaP_BG15n_", "LNCaP_XIP_BG15a_", "LNCaP_XIP_BG15n_", "X22Rv1_BG15a_", "X22Rv1_BG15n_", "LuCaP170_BG15n_", "LuCaP167_BG15n_"))), aes(x = Model_Tx, y = NAME, size = -log10.NOM.p.val,colour = NES)) + 
  geom_point() +scale_colour_gradient2(low="red",mid="white", high="green4") +
  ylab("GSEA Hallmark Pathways") + xlab("Sample")+
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
#levels = c("LNCaPXIP_BG15n_10d1", "LNCaPXIP_BG15n_20d1", "LNCaPXIP_BG15n_30d1", "LNCaPXIP_BG15n_10d2", "LNCaPXIP_BG15n_20d2", "LNCaPXIP_BG15n_30d2"))), aes(x = Model_Tx, y = NAME, size = -log10.NOM.p.val,colour = NES)) + 
  

#ggplot(allsamples, aes(x = Model_Tx, y = NAME, size = -log10.NOM.p.val,colour = NES)) + geom_point() +scale_colour_gradient2(low="red",mid="white", high="green4") +
# theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

#import select list of gene sets
#gene.sets=read.table("shSNAI2/shSNAI2_consistentHALLMARKS.txt")
#ggplot(subset(allsamples,(allsamples$NAME %in% gene.sets$V1)), aes(x = Model_Tx, y = NAME, size = -log10.NOM.p.val,colour,colour = NES)) + geom_point() +scale_colour_gradient2(low="grey20",mid="white", high="orange") +
#  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#Heatmap of enrichments
library(pheatmap)
library(tidyr)
library(tibble)
library(viridis)
?viridis
allsamples.wide   = allsamples[,c("NAME","NES","Model_Tx")]
allsamples.wide   = spread(allsamples.wide, NAME, NES)
allsamples.matrix = as.matrix(allsamples.wide[,2:ncol(allsamples.wide)])
dimnames(allsamples.matrix) = list( c(allsamples.wide$Model_Tx), c(allsamples$NAME))
rownames(allsamples.matrix) = allsamples.wide$Model_Tx
mycolors = colorRampPalette(c("firebrick", "red", "pink", "white","dodgerblue","blue"))
pheatmap(allsamples.matrix,scale = "none", cluster_cols= T, cluster_rows = F,  show_rownames = TRUE, show_colnames = TRUE, color = mycolors(50),  main="GSEA Hallmark Pathways", border_color = NA) 

library(dendsort)
mat_cluster_cols <- hclust(dist(t(allsamples.matrix)))
plot(mat_cluster_cols, main = "Unsorted Dendrogram", xlab = "", sub = "")
sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
mat_cluster_cols <- sort_hclust(mat_cluster_cols)
plot(mat_cluster_cols, main = "Sorted Dendrogram", xlab = "", sub = "")
mat_cluster_rows <- sort_hclust(hclust(dist(allsamples.matrix)))

pheatmap(allsamples.matrix,scale='none',cluster_cols=mat_cluster_cols,cluster_rows=F,  main="Sorted GSEA Hallmarks",color = mycolors(50),border_color = NA)
# b. GSEA ranklist plot for 1 gene set across many samples
#extract out RANK and RUNNING ES data for a given gene set
Geneset = "HALLMARK_ANDROGEN_RESPONSE"
Geneset.allsamples = data.frame(RANK.IN.GENE.LIST=numeric(), RUNNING.ES=numeric(), SYMBOL=character(), Model_Tx=character())

##HEATMAP A SUBSET OF PATHWAYS##
miniGSEA<- subset(allsamples, allsamples$NAME == c("ANDROGEN_RESPONSE"))
miniGSEA<- subset(allsamples, allsamples$NAME == c("ANDROGEN_RESPONSE") | allsamples$NAME ==c("MYC_TARGETS_V1") | allsamples$NAME ==c("MYC_TARGETS_V2")| allsamples$NAME ==c("G2M_CHECKPOINT")| allsamples$NAME ==c("E2F_TARGETS"))
miniGSEA.wide   = miniGSEA[,c("NAME","NES","Model_Tx")]
miniGSEA.wide   = spread(miniGSEA.wide, NAME, NES)
miniGSEA.matrix = as.matrix(miniGSEA.wide[,2:ncol(miniGSEA.wide)])
dimnames(miniGSEA.matrix) = list( c(miniGSEA.wide$Model_Tx), c(miniGSEA$NAME))
rownames(miniGSEA.matrix) = miniGSEA.wide$Model_Tx
mycolors = colorRampPalette(c("firebrick", "red", "pink", "white","dodgerblue","blue"))
pheatmap(miniGSEA.matrix,scale = "none", cluster_cols= F, cluster_rows = F,  show_rownames = TRUE, show_colnames = TRUE, color = mycolors(50),  main="GSEA Hallmark Pathways", border_color = NA) 


#project.folder = paste("projects/GSEA_ranklist")
setwd("Q:/RNA_seq/hg19/projects/GSEA_ranklist")
GSEA.folder = "ARIA_all"

lapply(sample.df$Model_Tx, function(x) {
  #read in fake xls file
  ESplot.sample <- subset(sample.df,sample.df$Model_Tx %in% x)
  
  ESplot.sample$geneset.path = paste(GSEA.folder,"/", paste(ESplot.sample$Model_Tx, ESplot.sample$GSEA_Style, ESplot.sample$Digits, sep="."),"/", Geneset,".tsv",sep="")
  
  Geneset.df = read.table(ESplot.sample$geneset.path,sep="\t",header=T)
  
  #make it smaller, add column for identification
  Geneset.df = Geneset.df[,c("RANK.IN.GENE.LIST","RUNNING.ES","SYMBOL")]
  
  Geneset.df$Model_Tx = x
  
  Geneset.allsamples <<- rbind(Geneset.allsamples,Geneset.df)
})

downgenesMYC <- file.path(project.folder, sample.set, paste(sample.set,".MYCdown.GSEA_matrix.txt"))
write.table(Geneset.allsamples, file = downgenesMYC, sep = "\t", row.names = F, col.names = T, quote = FALSE, append = FALSE)
downMYC = read.table(file = file.choose(), header = T) #sample.set,".MYCdown.GSEA_matrix.txt"
#X22Rv1_BG15a_downMYC = read.table(file = downgenesMYC, sep="\t", header=T, nrows=59)
#X22Rv1_BG15n_downMYC = read.table(file = downgenesMYC, sep="\t",  header=T, skip = 59, nrows=57)
#colnames(X22Rv1_BG15n_downMYC) <- c("RANK.IN.GENE.LIST", "RUNNING.ES", "SYMBOL", "Model_Tx")
#downgenesMYC_wide <- reshape(downgenesMYC, idvar = "Model_Tx", timevar = "SYMBOL", direction = "wide")#pivot_wider(downMYC,names_from = "Model_Tx", values_from = "SYMBOL" )


#Geneset.allsamples = na.omit(Geneset.allsamples)
#Geneset.allsamples = subset(Geneset.allsamples, Geneset.allsamples$Model_Tx %in% c("MKL2_LQ","MKL2_Mocet","MKL2_Pano"))
#Geneset.allsamples.pick = subset(Geneset.allsamples, Geneset.allsamples$SYMBOL %in% c("MYC"))
library(grid); library(gridExtra); library(ggplot2)
Model_Txcolors = c("tomato","red","blue","darkcyan","red4","dodgerblue","green","goldenrod","purple","magenta") # http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf

p1 = ggplot(Geneset.allsamples,aes(x = RANK.IN.GENE.LIST, y = RUNNING.ES, group = Model_Tx, color = Model_Tx)) + geom_line(size=1.2) + geom_hline(yintercept=0,lty='dashed')+ scale_color_manual(values=Model_Txcolors)+
  theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position = c(0.25, 0.25))+
  ylab("Enrichment Score") + xlab(paste("Genes Ranked by log2 fold change")) +ggtitle(paste("Geneset:", Geneset,sep=" "))#+  
#geom_point(data=Geneset.allsamples.pick, color = "black")+geom_text(data=Geneset.allsamples.pick, label = Geneset.allsamples.pick$SYMBOL)

p2 = ggplot(Geneset.allsamples,group = Model_Tx) + stat_density(aes(x=RANK.IN.GENE.LIST,y=0.5,fill=..density..),geom="tile",position="identity")+facet_wrap(~Model_Tx, ncol=1)+ scale_fill_gradient(low="white",high="darkgrey")+
  theme(axis.line=element_blank(),
        axis.title.x=element_blank(), legend.position="none",panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

p3 = p2+ geom_linerange(aes(x=RANK.IN.GENE.LIST,ymin=0,ymax=1,color = Model_Tx))+facet_wrap(~Model_Tx, ncol=1)+theme_bw() + scale_color_manual(values=Model_Txcolors)+
  theme(axis.line=element_blank(),
        axis.title.x=element_blank(), legend.position="none",panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

grid.arrange(arrangeGrob(p1,ncol=1))
grid.arrange(arrangeGrob(p1,p3,ncol=1))


#
Geneset.edge = subset(Geneset.allsamples, Geneset.allsamples$RANK.IN.GENE.LIST > 7000)
Geneset.edge = Geneset.edge[-grep("ENZA",Geneset.edge$Model_Tx),]
Geneset.edge = unique(Geneset.edge$SYMBOL)
#leading edge genes
LeadEdge.all = data.frame(SYMBOL=character(),RANK.IN.GENE.LIST=numeric(),RANK.METRIC.SCORE=numeric(),Model_Tx=character())
lapply(sample.df$Model_Tx, function(x) {
  #read in fake xls file
  ESplot.sample <- subset(sample.df,sample.df$Model_Tx %in% x)
  
  ESplot.sample$geneset.path = paste(GSEA.folder,"/",paste(ESplot.sample$Model_Tx,ESplot.sample$GSEA_Style,ESplot.sample$Digits,sep="."),"/",Geneset,".xls",sep="")
  
  Geneset.df = read.table(ESplot.sample$geneset.path,sep="\t",header=T)
  
  #make it smaller, add column for identification
  Geneset.df = Geneset.df[,c("SYMBOL","RANK.IN.GENE.LIST","RANK.METRIC.SCORE")]
  
  Geneset.df$Model_Tx = x
  
  LeadEdge.all <<- rbind(LeadEdge.all,Geneset.df)
})

LeadEdge.rep1 = subset(LeadEdge.all, LeadEdge.all$Model_Tx %in% "RMS_DLp300_lo")
LeadEdge.rep2 = subset(LeadEdge.all, LeadEdge.all$Model_Tx %in% "RMS_DLp300_hi")
LeadEdge.comb = merge(LeadEdge.rep1, LeadEdge.rep2, by = "SYMBOL")
LeadEdge.pick = subset(LeadEdge.comb, LeadEdge.comb$SYMBOL %in% c("MYC","NOLC1","PNP"))
ggplot(LeadEdge.comb, aes(x=RANK.METRIC.SCORE.x, y=RANK.METRIC.SCORE.y))+geom_point()+
  geom_point(data=LeadEdge.pick, color = "red")+geom_text(data=LeadEdge.pick, label = LeadEdge.pick$SYMBOL)
ggplot(LeadEdge.comb, aes(x=RANK.IN.GENE.LIST.x, y=RANK.IN.GENE.LIST.y))+
  geom_point(data=LeadEdge.pick, color = "red")+geom_text(data=LeadEdge.pick, label = LeadEdge.pick$SYMBOL)+
  theme_bw()+coord_fixed(ratio=1)+xlim(1, 19000)+ylim(1, 19000)

LeadEdge.pick = subset(LeadEdge.all, LeadEdge.all$SYMBOL %in% c("MYC","NOLC1","PNP"))
ggplot(LeadEdge.pick, aes(x=Model_Tx, y=RANK.IN.GENE.LIST))+
  geom_point(data=LeadEdge.pick, color = "red")+geom_text(data=LeadEdge.pick, label = LeadEdge.pick$SYMBOL)+
  theme_bw()+ylim(1, 19000)+facet_wrap(~SYMBOL)



## Setup: define libraries, set project folder, read TPM matrix and gene list (with optional sample subsetting)
library(plyr)
library(dplyr)
library(pheatmap)
library(dendsort)

project.folder = "/Users/matthewchang/Documents/R/"

setwd("/Users/matthewchang/Documents/R/ExpMatrices/")
EXP.coding.matrix = read.table("IHK_samples.coding.norm.matrix.txt", header = T)                                # read in TPM matrix

setwd("/Users/matthewchang/Documents/R/GeneSets/")
gene.list = read.table(file.choose(), sep = "\t", header = F)                                                   # read in gene list

## subset samples
EXP.coding.matrix = EXP.coding.matrix %>% select(-((ncol(EXP.coding.matrix) - 3):ncol(EXP.coding.matrix)))      # dropping the last 4 columns (IHK45)





### 4. Heatmap
### Goal: generate heatmap to visualize RNA-seq data across all samples
### Pre-requisites: TPM matrix, gene list

## Filtering based off of minimum expression and genes
cutoff.expression.min = 10                                                                                      # set the cutoff threshold for minimal expression

EXP.expressed.matrix = EXP.coding.matrix                                                                        # create a copy of the expression matrix for filtering
EXP.expressed.matrix$maxTPM = apply(EXP.expressed.matrix[, 2:(ncol(EXP.expressed.matrix) - 1)], 1, FUN = max)   # calculate the maximum TPM value for each gene across all samples
EXP.expressed.matrix = subset(EXP.expressed.matrix, EXP.expressed.matrix$maxTPM > cutoff.expression.min)        # subset the matrix to only include genes with a maximum TPM value above the threshold
EXP.expressed.matrix.TFs = subset(EXP.expressed.matrix, EXP.expressed.matrix$gene_id %in% gene.list$V1)         # subset the matrix to only include genes in the gene.list
rownames(EXP.expressed.matrix.TFs) = EXP.expressed.matrix.TFs$gene_id                                           # rownames set to gene_id
EXP.expressed.matrix.TFs = EXP.expressed.matrix.TFs[, 2:(ncol(EXP.expressed.matrix.TFs) - 1)]                   # dropped first (gene_id) and last (maxTPM)
EXP.expressed.matrix.TFs = as.matrix(t(EXP.expressed.matrix.TFs))                                               # transposed matrix
rownames(EXP.expressed.matrix.TFs) = gsub("RH4_", "", rownames(EXP.expressed.matrix.TFs))                       # removed "RH4_"
rownames(EXP.expressed.matrix.TFs) = gsub("_RNA_022924_CWRU", "", rownames(EXP.expressed.matrix.TFs))           # removed "_RNA_022924_CWRU"

## OPTIONAL: separate samples by group
annotations = data.frame(
  SampleGroup = c("Control/selective degraders", "Control/selective degraders", "Control/selective degraders", "Control/selective degraders", 
                  "Dual inhibitors/degraders", "Dual inhibitors/degraders", "Dual inhibitors/degraders", "Dual inhibitors/degraders", 
                  "Control/selective degraders", "Control/selective degraders", "Control/selective degraders", "Control/selective degraders",
                  "Dual inhibitors/degraders", "Dual inhibitors/degraders", "Dual inhibitors/degraders", "Dual inhibitors/degraders",
                  "Control/selective degraders", "Control/selective degraders", "Control/selective degraders", "Control/selective degraders", 
                  "Control/selective degraders", "Control/selective degraders", "Control/selective degraders", "Control/selective degraders", 
                  "Dual inhibitors/degraders", "Dual inhibitors/degraders", "Dual inhibitors/degraders", "Dual inhibitors/degraders")
)
rownames(annotations) = rownames(EXP.expressed.matrix.TFs)
annotation_colors = list(
  SampleGroup = c("Control/selective degraders" = "darkgoldenrod1", 
                  "Dual inhibitors/degraders" = "tomato")
)

# Generate the custom color scale from white to blue
color_scale = colorRampPalette(c("gray99", "plum2", "mediumvioletred"))(100)

pheatmap(log2(EXP.expressed.matrix.TFs + 1), 
         cluster_rows = T, 
         cluster_cols = T, 
         scale = 'none', 
         show_rownames = T, 
         show_colnames = T, 
         angle_col = 315, 
         annotation_row = annotations,
         main = paste("Log2(TPM + 1)"),
         color = color_scale,
         annotation_colors = annotation_colors)













## heatmap of log2FC vs DMSO
setwd("/Users/matthewchang/Documents/R/GSEA_log2FC/p300_samples/")
EXP.log2FC.heat = read.table("EXP.log2FC.txt", header = T, row.names = 1)                                       # read in log2FC values

colnames(EXP.log2FC.heat) = gsub("^RH4_", "", colnames(EXP.log2FC.heat))                                        # removed "RH4_"
EXP.log2FC.heat = EXP.log2FC.heat[rownames(EXP.log2FC) %in% gene.list$V1, ]                                     # filter genes from gene.list
EXP.log2FC.heat = as.matrix(t(EXP.log2FC.heat))                                                                 # transposed matrix

## OPTIONAL: separate samples by group
annotations = data.frame(
  SampleGroup = c("Dual inhibitors/degraders", "Dual inhibitors/degraders", "Dual inhibitors/degraders", "Dual inhibitors/degraders", 
                  "Selective degraders", "Selective degraders", "Selective degraders", "Selective degraders",
                  "Dual inhibitors/degraders", "Dual inhibitors/degraders", "Dual inhibitors/degraders", "Dual inhibitors/degraders",
                  "Selective degraders", "Selective degraders", "Selective degraders", "Selective degraders", 
                  "Selective degraders", "Selective degraders", "Selective degraders", "Selective degraders", 
                  "Dual inhibitors/degraders", "Dual inhibitors/degraders", "Dual inhibitors/degraders", "Dual inhibitors/degraders")
)
rownames(annotations) = rownames(EXP.log2FC.heat)
annotation_colors = list(
  SampleGroup = c("Selective degraders" = "darkgoldenrod1", 
                  "Dual inhibitors/degraders" = "tomato")
)

color_scale = colorRampPalette(c("cornflowerblue", "cadetblue3", "gray99"))(100)

pheatmap(EXP.log2FC.heat, 
         cluster_rows = T, 
         cluster_cols = T, 
         scale = 'none', 
         show_rownames = T, 
         show_colnames = T, 
         angle_col = 315, 
         annotation_row = annotations,
         main = paste("Log2FC"),
         color = color_scale,
         annotation_colors = annotation_colors)



######################### end of heatmap #########################








### 5. PCA
### Goal: Perform PCA to visualize sample clustering.
### Prerequisites: TPM matrix, sample list with subtype/drug Tx and timepoint info

library(tidyverse)
library(ggrepel)

# create PCA matrix
EXP.pca = as.matrix(EXP.coding.matrix[, 2:ncol(EXP.coding.matrix)])
rownames(EXP.pca) = EXP.coding.matrix$gene_id                                                                   # rownames set to gene_id
EXP.pca.log2 = log2(EXP.pca + 1)                                                                                # apply log2 transformation
sample.name.list = c(colnames(EXP.pca))                                                                         # extract sample names

## OPTIONAL: color and shape grouping for drug tx and timepoint
setwd("/Users/matthewchang/Documents/R/SampleList/")
sample_class = read.table(file.choose(), header = F)
drugs = sample_class$V2[match(sample.name.list, sample_class$V1)]                                               # acquire drug tx info
timepoints = sample_class$V3[match(sample.name.list, sample_class$V1)]                                          # acquire timepoint info
timepoints = as.factor(timepoints)                                                                              # convert timepoints to a factor
EXP.coldata = data.frame(sample.name.list, drugs, timepoints)                                                   # create df with sample name, drug tx, and timepoint info
rownames(EXP.coldata) = sample.name.list                                                                        # rownames set to sample name

## Perform PCA
pca = EXP.pca.log2 %>% t %>% prcomp                                                                             # pca function
EXP.pca.df = pca$x %>% as.data.frame                                                                            # convert to df
EXP.pca.df$sample.name.list = EXP.coldata$sample.name.list                                                      # add sample names

# Combine the PCA data frame with the sample metadata
EXP.pca.df.meta = join(EXP.pca.df, EXP.coldata, by = "sample.name.list")                                        # add in EXP.coldata from sample.name.list
EXP.pca.df.meta$sample.name.list = gsub("_RNA_022924_CWRU", "", EXP.pca.df.meta$sample.name.list)               # removed "_RNA_022924_CWRU"
EXP.pca.df.meta$sample.name.list = gsub("RH4_", "", EXP.pca.df.meta$sample.name.list)                           # removed "RH4_"
EXP.pca.df.meta = na.omit(EXP.pca.df.meta)                                                                      # omit missing values

# Calculate the proportion of variability explained by each principal component
pcv = round((pca$sdev)^2 / sum(pca$sdev^2) * 100, 2)
      
# Graph the PCA plot
plot.pca = ggplot(EXP.pca.df.meta, aes(PC1, PC2, colour = drugs, shape = timepoints)) + 
    geom_point() +
    xlab(label = paste0("PC1 (", pcv[1], "%)")) +
    ylab(label = paste0("PC2 (", pcv[2], "%)")) +
    theme_bw() +
    geom_label_repel(aes(label = sample.name.list), show.legend = FALSE) +
    theme(axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15)) +
    labs(title = "PCA",
         subtitle = "Grouped by drug treatment and timepoint")
print(plot.pca)

######################### end of PCA #########################










### 5. make GSEA ranklists
### Goal: To create rank lists for Gene Set Enrichment Analysis (GSEA) by calculating log2 fold changes (log2FC) for comparisons.
### Pre-requisites: TPM matrix
### Notes: You must specify the log2FC for the samples of interest by hand

## Load in TPM expression matrix
setwd("/Users/matthewchang/Documents/R/ExpMatrices/")
project.folder = "/Users/matthewchang/Documents/R/"
sample.set = "p300_samples"
EXP.coding.matrix = read.table("IHK_samples.coding.norm.matrix.txt", header = T)

## Create a directory for the GSEA ranklists and create the log2FC matrix
dir.create(file.path(project.folder, "GSEA_ranklist", sample.set), recursive = T)
EXP.GSEA = EXP.coding.matrix

## Calculate log2FC for each sample and add as new column
# EXP.GSEA$RH4_DMSO_2h = log2(EXP.GSEA$RH4_DMSO_2h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_DMSO_2h_RNA_022924_CWRU + 1)
# EXP.GSEA$RH4_DMSO_6h = log2(EXP.GSEA$RH4_DMSO_6h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_DMSO_6h_RNA_022924_CWRU + 1)
EXP.GSEA$RH4_A485_100nM_2h = log2(EXP.GSEA$RH4_A485_100nM_2h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_DMSO_2h_RNA_022924_CWRU + 1)
EXP.GSEA$RH4_A485_1uM_2h = log2(EXP.GSEA$RH4_A485_1uM_2h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_DMSO_2h_RNA_022924_CWRU + 1)
EXP.GSEA$RH4_A485_100nM_6h = log2(EXP.GSEA$RH4_A485_100nM_6h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_DMSO_6h_RNA_022924_CWRU + 1)
EXP.GSEA$RH4_A485_1uM_6h = log2(EXP.GSEA$RH4_A485_1uM_6h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_DMSO_6h_RNA_022924_CWRU + 1)
EXP.GSEA$RH4_JQAD_100nM_2h = log2(EXP.GSEA$RH4_JQAD_100nM_2h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_DMSO_2h_RNA_022924_CWRU + 1)
EXP.GSEA$RH4_JQAD_1uM_2h = log2(EXP.GSEA$RH4_JQAD_1uM_2h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_DMSO_2h_RNA_022924_CWRU + 1)
EXP.GSEA$RH4_JQAD_100nM_6h = log2(EXP.GSEA$RH4_JQAD_100nM_6h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_DMSO_6h_RNA_022924_CWRU + 1)
EXP.GSEA$RH4_JQAD_1uM_6h = log2(EXP.GSEA$RH4_JQAD_1uM_6h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_DMSO_6h_RNA_022924_CWRU + 1)
EXP.GSEA$RH4_dCBP_100nM_2h = log2(EXP.GSEA$RH4_dCBP_100nM_2h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_DMSO_2h_RNA_022924_CWRU + 1)
EXP.GSEA$RH4_dCBP_1uM_2h = log2(EXP.GSEA$RH4_dCBP_1uM_2h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_DMSO_2h_RNA_022924_CWRU + 1)
EXP.GSEA$RH4_dCBP_100nM_6h = log2(EXP.GSEA$RH4_dCBP_100nM_6h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_DMSO_6h_RNA_022924_CWRU + 1)
EXP.GSEA$RH4_dCBP_1uM_6h = log2(EXP.GSEA$RH4_dCBP_1uM_6h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_DMSO_6h_RNA_022924_CWRU + 1)
EXP.GSEA$RH4_QL_100nM_2h = log2(EXP.GSEA$RH4_QL_100nM_2h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_DMSO_2h_RNA_022924_CWRU + 1)
EXP.GSEA$RH4_QL_1uM_2h = log2(EXP.GSEA$RH4_QL_1uM_2h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_DMSO_2h_RNA_022924_CWRU + 1)
EXP.GSEA$RH4_QL_100nM_6h = log2(EXP.GSEA$RH4_QL_100nM_6h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_DMSO_6h_RNA_022924_CWRU + 1)
EXP.GSEA$RH4_QL_1uM_6h = log2(EXP.GSEA$RH4_QL_1uM_6h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_DMSO_6h_RNA_022924_CWRU + 1)
EXP.GSEA$RH4_LS_100nM_2h = log2(EXP.GSEA$RH4_LS_100nM_2h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_DMSO_2h_RNA_022924_CWRU + 1)
EXP.GSEA$RH4_LS_1uM_2h = log2(EXP.GSEA$RH4_LS_1uM_2h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_DMSO_2h_RNA_022924_CWRU + 1)
EXP.GSEA$RH4_LS_100nM_6h = log2(EXP.GSEA$RH4_LS_100nM_6h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_DMSO_6h_RNA_022924_CWRU + 1)
EXP.GSEA$RH4_LS_1uM_6h = log2(EXP.GSEA$RH4_LS_1uM_6h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_DMSO_6h_RNA_022924_CWRU + 1)
EXP.GSEA$RH4_IHK44_100nM_2h = log2(EXP.GSEA$RH4_IHK44_100nM_2h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_DMSO_2h_RNA_022924_CWRU + 1)
EXP.GSEA$RH4_IHK44_1uM_2h = log2(EXP.GSEA$RH4_IHK44_1uM_2h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_DMSO_2h_RNA_022924_CWRU + 1)
EXP.GSEA$RH4_IHK44_100nM_6h = log2(EXP.GSEA$RH4_IHK44_100nM_6h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_DMSO_6h_RNA_022924_CWRU + 1)
EXP.GSEA$RH4_IHK44_1uM_6h = log2(EXP.GSEA$RH4_IHK44_1uM_6h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_DMSO_6h_RNA_022924_CWRU + 1)
EXP.GSEA$RH4_IHK45_100nM_2h = log2(EXP.GSEA$RH4_IHK45_100nM_2h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_DMSO_2h_RNA_022924_CWRU + 1)
EXP.GSEA$RH4_IHK45_1uM_2h = log2(EXP.GSEA$RH4_IHK45_1uM_2h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_DMSO_2h_RNA_022924_CWRU + 1)
EXP.GSEA$RH4_IHK45_100nM_6h = log2(EXP.GSEA$RH4_IHK45_100nM_6h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_DMSO_6h_RNA_022924_CWRU + 1)
EXP.GSEA$RH4_IHK45_1uM_6h = log2(EXP.GSEA$RH4_IHK45_1uM_6h_RNA_022924_CWRU + 1) - log2(EXP.GSEA$RH4_DMSO_6h_RNA_022924_CWRU + 1)

## Round the log2FC to 5 sig figs
EXP.GSEA[,(ncol(EXP.coding.matrix) + 1) : ncol(EXP.GSEA)] = round(EXP.GSEA[,(ncol(EXP.coding.matrix) + 1) : ncol(EXP.GSEA)], digits = 5)

## Loop through the new columns, create rank lists, and save them to files
for (i in (ncol(EXP.coding.matrix) + 1):ncol(EXP.GSEA)) {
  Ranklist <- data.frame(EXP.GSEA[, 1])                                                                         # initialize rank list with gene names
  Ranklist$DeltaTPM <- EXP.GSEA[, i]                                                                            # add the log2FC values to the rank list
  Ranklist = Ranklist[rev(order(Ranklist$DeltaTPM)), ]                                                          # sort the rank list by DeltaTPM in descending order
  SampleName = colnames(EXP.GSEA)[i]                                                                            # get the name of the current comparison
  mytime <- format(Sys.time(), "%b_%d_%Y")                                                                      # get the current date
  myfile <- file.path(project.folder, "GSEA_ranklist", sample.set, paste0(SampleName,"_",mytime,".rnk"))        # define file path
  write.table(Ranklist, file = myfile, sep = "\t", row.names = F, col.names = F, quote = F, append = F)         # save the ranklist
}

## In GSEA, use the following datasets and settings
#  Gene sets database: GRYDERLAB_RMS_gene_sets.gmt
#  Collapse/Remap to gene symbols: No_collapse
#  Chip platoform: Human_Gene_Symbol_with_Remapping_MSigDBv2023.2.Hs.chip
#  Max size: 5000
#  Min size: 10

# OPTIONAL: export log2FC df to project.folder
exported.columns = EXP.GSEA[, 30:ncol(EXP.GSEA)]
rownames(exported.columns) = EXP.coding.matrix$gene_id
export.file = file.path(project.folder, "GSEA_log2FC", sample.set, "EXP.log2FC.txt")
dir.create(file.path(project.folder, "GSEA_log2FC", sample.set), recursive = T)
write.table(exported.columns, file = export.file, sep = "\t", row.names = T, col.names = T, quote = F)


######################### end of making GSEA ranklists #########################

##################################################
###  build transcripts per million data matrix ###
###  then perform GSEA, heatmap comparisons,   ###
###  PCA plots, ranked scatter plots, and more ###
##################################################
###  by Berkley Gryder, gryderart@gmail.com    ###
###  edited by Diana 2022 for Mac              ###
###  edited by Bhava 2023 for ###4             ###
###  edited by Matt 2024                      ###
##################################################

### Goal: build TPM expression matrix from RNA-seq data

### Pre-requisites: .fastq files have been run through the RNA-seq pipeline

### 1. Set directories, define samples to work on.
getwd()
setwd("/Volumes/rc/SOM_GENE_BEG33/RNA_seq/hg38/projects/IHK_RMS/Practice_MSC")

## Read sample list from a file chosen from the client
## Note: this file should be a .txt file with the first entry labeled "Sample"
sample.list.all = read.table(file.choose(), header = T, sep = "\t")                           # acquiring the sample list
sample.list = sample.list.all$Sample
project.folder = "/Volumes/rc/SOM_GENE_BEG33/RNA_seq/hg38/projects/RMS_IHK/Practice_MSC"      # define the project folder path
sample.set = "IHK_samples"                                                                    # define the sample set name




### 2. Make Protein coding TPM file, normalizing based on expected counts

setwd("/Volumes/rc/SOM_GENE_BEG33/RNA_seq/hg38/")

## Check to see if TPM.txt exist for each sample in the sample list
## Remove samples that do not have the necessary files
file.exists = file.exists(paste0("DATA/", sample.list,"/", sample.list, ".gene.TPM.txt", sep="")) 
sample.list = sample.list[file.exists]

## Loop through each sample in the sample list
lapply(sample.list, function(x) {
  
    ## Load more stringent protein coding list (from Diana)
    coding <- read.table("ref/RSEM_hg38/HGNC_protein-coding_gene_19229_2022.txt", sep="\t", header=T)
    
    ## Load RSEM (genes.results) output file for the current sample with expected count, TPM, FPKM
    EXP <- read.table(paste("DATA/",x,"/",x,".genes.results",sep=""), sep="\t", header=T)
    
    ## Remove non-coding RNA entries
    EXP$coding = EXP$gene_id %in% coding$symbol
    EXP.coding <<- subset(EXP, EXP$coding %in% c("TRUE"))
    
    ## Sum all coding expected counts
    expected_sum = sum(EXP.coding$expected_count)
    
    ## Normalize expected counts to counts per million
    EXP.coding$count_norm = EXP.coding$expected_count / expected_sum * 1000000
    
    ## Write the normalized data to a file
    write.table(EXP.coding[,1:9], 
                file=paste("DATA/",x,"/",x,".gene.coding.norm.txt",sep=""), 
                sep="\t", 
                row.names=FALSE, 
                col.names=TRUE, 
                quote=FALSE)
})






### 3. Build matrix


## Initiate matrix with gene IDs
EXP.coding.matrix = EXP.coding["gene_id"]

## Loop through each sample, again, to extract TPM (normalized counts)
lapply(sample.list, function(x) {
  
    ## Read the normalized counts column
    EXP.sample <- read.table(paste("DATA/",x,"/",x,".gene.coding.norm.txt",sep=""), sep="\t", header=T)
  
    ## Extract the normalized counts column
    EXP.sample = as.data.frame(EXP.sample[,9])
    removable.string = "Sample_"
    sample.name = gsub(removable.string,"",x)
    colnames(EXP.sample) = c(sample.name)
    
    ## Combine the current sample data into the matrix
    EXP.coding.matrix <<- cbind(EXP.coding.matrix, EXP.sample)
})

colnames(EXP.coding.matrix) <- gsub("-", "", colnames(EXP.coding.matrix))                       # remove hyphen

## Write the final expression matrix to a file
## Note: ExpMatrices directory must exist in project.folder before executing

## Note: Can create directory
dir.create(file.path(project.folder, "ExpMatrices"))

write.table(EXP.coding.matrix, 
            file=paste(project.folder,"/", "ExpMatrices/",sample.set,".coding.norm.matrix.txt",sep=""), 
            sep="\t", 
            row.names=FALSE, 
            col.names=TRUE, 
            quote=FALSE)


################### end of building the TPM matrix ###################















### 4. Make some comparison plots

### load in data, set parameters
setwd("/Volumes/rc/SOM_GENE_BEG33/RNA_seq/hg38/projects/RMS_IHK/Practice_MSC/ExpMatrices")
project.folder = "/Volumes/rc/SOM_GENE_BEG33/RNA_seq/hg38/projects/RMS_IHK/Practice_MSC" 

#### 4.1. heatmaps
cutoff.expression.min = 10
EXP.expressed.matrix = EXP.coding.matrix
#EXP.expressed.matrix_2 = EXP.coding.matrix_2
EXP.expressed.matrix$maxTPM = apply(EXP.expressed.matrix[,c(2:ncol(EXP.expressed.matrix))], 1, FUN=max)
EXP.expressed.matrix = subset(EXP.expressed.matrix, EXP.expressed.matrix$maxTPM > cutoff.expression.min)
AllHumanTFs = read.table(file = "/Volumes/rc/SOM_GENE_BEG33/RNA_seq/hg38/projects/Genesets/List_format/TranscriptionFactorsGryder.txt", sep="\t", header=F)
EXP.expressed.matrix.TFs = subset(EXP.expressed.matrix, EXP.expressed.matrix$gene_id %in% AllHumanTFs$V1)
library(pheatmap)
pheatmap(log2(as.matrix(EXP.expressed.matrix[,2:(ncol(EXP.expressed.matrix.TFs)-1)])+1),scale='row')


##Add the subtype information and give a different color to each subtype
library(RColorBrewer)
library(viridis)
sample.info = sample.list.all
sample.diagnosis = sample.info$Diagnosis[match(as.character(colnames(EXP.expressed.matrix)), as.character(sample.info$Sample))]
sample.diagnosis = sample.diagnosis[-c(1, 349)]
sample.diagnosis[is.na(sample.diagnosis)] <- "Neuroblastoma"

diagnosis = unique(sample.diagnosis)

group_colors <- brewer.pal(length(diagnosis), "Set1")
names(group_colors) <- diagnosis
mat_colors <- list(Group = group_colors)

## All Coding ##
#Dendrogram sort
# install.packages("dendsort")
library(dendsort)
mat_cluster_cols <- hclust(dist(t(EXP.expressed.matrix[,c(2:(ncol(EXP.expressed.matrix)-1))]))) ##include only the columns with numeric values
plot(mat_cluster_cols, main = "Unsorted Dendrogram", xlab = "", sub = "")
sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
mat_cluster_cols <- sort_hclust(mat_cluster_cols)
plot(mat_cluster_cols, main = "Sorted Dendrogram", xlab = "", sub = "")

# library(dendextend)
# colored_dend <- color_branches(mat_cluster_cols, k = length(sample.diagnosis), col = group_colors[diagnosis])
# colored_dend <- color_labels(colored_dend, col = group_colors[diagnosis])
# plot(colored_dend, main = "Sorted Dendrogram", xlab = "", sub = "")
# 
# mat_cluster_rows <- sort_hclust(hclust(dist(EXP.expressed.matrix[,(2:110)])))

pheatmap(EXP.expressed.matrix[,2:(ncol(EXP.expressed.matrix.TFs)-1)],scale='row',cluster_cols=mat_cluster_cols, cluster_rows=T, main=paste(sample.set," "," EXP matrix heatmap_all coding genes",sep=""))

##Adding cancer diagnosis information
sample.diagnosis.df <- as.data.frame(sample.diagnosis)
rownames(sample.diagnosis.df) <- colnames(EXP.expressed.matrix)[2:(ncol(EXP.expressed.matrix)-1)]

pheatmap(EXP.expressed.matrix[,2:(ncol(EXP.expressed.matrix.TFs)-1)],scale='row',
         cluster_cols=mat_cluster_cols, 
         cluster_rows=T, 
         annotation_col = sample.diagnosis.df, 
         main=paste(sample.set," "," EXP matrix heatmap_all coding genes",sep=""))
pheatmap(EXP.expressed.matrix.TFs[,2:(ncol(EXP.expressed.matrix.TFs)-1)],scale='row',
         cluster_cols=mat_cluster_cols, 
         cluster_rows=T, 
         annotation_col = sample.diagnosis.df, 
         main=paste(sample.set," "," EXP matrix heatmap_TFs",sep=""))

#### 4.2. PCA plot.  Metadata needs to be curated to meet your particular sample set
EXP.pca = as.matrix(EXP.coding.matrix[,c(2:ncol(EXP.coding.matrix))])
rownames(EXP.pca) = EXP.coding.matrix$geneID
EXP.pca.log2 = log2(EXP.pca+1)
sample.name.list = c(colnames(EXP.pca)) 
EXP.coldata = data.frame(sample.name.list);rownames(EXP.coldata) = sample.name.list

# Generate PCA Data & Proportion of variability
library(tidyverse) 
library(ggrepel)   
pca        <- EXP.pca.log2 %>% t %>% prcomp
EXP.d      <- pca$x %>% as.data.frame; EXP.d$sample.name.list = EXP.coldata$sample.name.list
library(plyr)
EXP.dmeta  <- join(EXP.d, EXP.coldata, by = "sample.name.list")
pcv        <- round((pca$sdev)^2 / sum(pca$sdev^2)*100, 2)

plot.pca   <- ggplot(EXP.dmeta, aes(PC2,PC3, colour = sample.diagnosis)) + #symbol = sample.list # can add shape = time.list
  geom_point() +
  xlab(label=paste0("PC2 (", pcv[2], "%)")) +
  ylab(label=paste0("PC3 (", pcv[3], "%)")) +
  theme_bw() +
  geom_label_repel(aes(label = sample.name.list), show.legend = F) +
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15)) +
  labs(title    = "PCA",
       subtitle = "pan-cancer")


print(plot.pca)

##### 4.2.1 pca for TFs
EXP.pca.TF = as.matrix(EXP.expressed.matrix.TFs[,2:(ncol(EXP.expressed.matrix.TFs)-1)])
rownames(EXP.pca.TF) = EXP.expressed.matrix.TFs$geneID
EXP.pca.log2.TF = log2(EXP.pca.TF+1)
sample.name.list.TF = c(colnames(EXP.pca.TF)) 
EXP.coldata.TF = data.frame(sample.name.list.TF);rownames(EXP.coldata.TF) = sample.name.list.TF

# Generate PCA Data & Proportion of variability
pca.TF        <- EXP.pca.log2.TF %>% t %>% prcomp
EXP.d.TF      <- pca.TF$x %>% as.data.frame; EXP.d.TF$sample.name.list.TF = EXP.coldata.TF$sample.name.list.TF
EXP.dmeta.TF  <- join(EXP.d.TF, EXP.coldata.TF, by = "sample.name.list.TF")
pcv.TF        <- round((pca.TF$sdev)^2 / sum(pca.TF$sdev^2)*100, 2)

#sample.subtypes <- sample.subtypes[-c(1, 111)] #removing the NA vectors

# Graph the PCA plot
plot.pca.TF   <- ggplot(EXP.dmeta.TF, aes(PC1,PC2, colour = sample.diagnosis)) + #symbol = sample.list # can add shape = time.list
  geom_point() +
  xlab(label=paste0("PC1 (", pcv[1], "%)")) +
  ylab(label=paste0("PC2 (", pcv[2], "%)")) +
  theme_bw() +
  geom_label_repel(aes(label = sample.name.list.TF), show.legend = F) +
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15)) +
  labs(title    = "PCA",
       subtitle = "Pan-cancer")

print(plot.pca.TF)








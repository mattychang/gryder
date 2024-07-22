## Setup: define libraries, set project folder, read TPM matrix and gene list (with optional sample subsetting)
library(plyr)
library(dplyr)
library(pheatmap)
library(dendsort)

project.folder = "/Users/matthewchang/Documents/R/"

setwd("/Users/matthewchang/Documents/R/ExpMatrices/")
EXP.coding.matrix = read.table("IHK_samples.coding.norm.matrix.txt", header = T)                                                       # read TPM matrix

setwd("/Users/matthewchang/Documents/R/GeneSets/")
gene.list = read.table(file.choose(), sep = "\t", header = F)                                                   # read gene list

## OPTIONAL: subset samples
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
rownames(annotations) = rownames(EXP.expressed.matrix.TFs)                                                      # match rownames in annotations with EXP.expressed.matrix.TFs

## Generate the heatmap (log2 transformation)
## OPTIONAL: z-score / scale = 'row' or 'none'
pheatmap(log2(EXP.expressed.matrix.TFs + 1), 
         cluster_rows = T, 
         cluster_cols = T, 
         scale = 'none', 
         show_rownames = T, 
         show_colnames = T, 
         angle_col = 315, 
         annotation_row = annotations,
         main = paste("Log2(TPM + 1)"))

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
setwd("/Volumes/rc/SOM_GENE_BEG33/RNA_seq/hg38/")
project.folder = "/Volumes/rc/SOM_GENE_BEG33/RNA_seq/hg38/projects/RMS_IHK/Practice_MSC"
sample.set = "p300_samples"
EXP.coding.matrix = read.table(file.choose(), header = T)

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
  myfile <- file.path(project.folder, "GSEA_ranklist", sample.set, paste0(SampleName,"_",mytime,".rnk"))                    # define file path
  write.table(Ranklist, file = myfile, sep = "\t", row.names = F, col.names = F, quote = F, append = F)         # save the ranklist
}

## In GSEA, use the following datasets and settings
#  Gene sets database: GRYDERLAB_RMS_gene_sets.gmt
#  Collapse/Remap to gene symbols: No_collapse
#  Chip platoform: Human_Gene_Symbol_with_Remapping_MSigDBv2023.2.Hs.chip
#  Max size: 5000
#  Min size: 10

######################### end of making GSEA ranklists #########################












# boxplot/violin plot by log2FC by condition separated by multiple genesets
# note: all coding genes did not show me anything good



# Load necessary libraries
library(tidyverse)

# Define the sample names of interest
sample_names_of_interest <- c("DMSO_6h", "A485_1uM_6h", "dCBP_1uM_6h", "IHK44_1uM_6h", "JQAD_1uM_6h", "LS_1uM_6h", "QL_1uM_6h")

# Set working directory and read the gene set of interest
setwd("/Users/matthewchang/Documents/R/GeneSets/")
RH4_CR_TFs <- read.table("GRYDER_RH4_CR_TFs.genelist.txt", header = FALSE, stringsAsFactors = FALSE)

# Create log2FC df
EXP.log2FC <- EXP.GSEA[, 30:ncol(EXP.GSEA)]
colnames(EXP.log2FC) <- gsub("^RH4_", "", colnames(EXP.log2FC))
rownames(EXP.log2FC) <- EXP.coding.matrix$gene_id

# Transpose the DataFrame
EXP.filtered_t <- as.data.frame(t(EXP.log2FC))
EXP.filtered_t <- rownames_to_column(EXP.filtered_t, var = "sample_name")

# Gather the data to long format
EXP_long <- EXP.filtered_t %>%
  gather(key = "gene", value = "log2FC", -sample_name)

# Filter for the specified sample names
EXP_filtered_samples <- EXP_long %>%
  filter(sample_name %in% sample_names_of_interest)

# Filter for genes in the specified gene set
EXP_filtered_samples <- EXP_filtered_samples %>%
  filter(gene %in% RH4_CR_TFs$V1)

# Set the order of sample names
EXP_filtered_samples$sample_name <- factor(EXP_filtered_samples$sample_name, levels = sample_names_of_interest)

# Create the violin plot with overlaid boxplot for the specified samples, without showing outliers
ggplot(EXP_filtered_samples, aes(x = sample_name, y = log2FC, fill = sample_name)) +
  geom_violin(alpha = 0.5) + 
  geom_boxplot(width = 0.1, color = "black", alpha = 0.7, outlier.shape = NA) +  # Boxplot without outliers
  labs(x = "Sample Name", 
       y = "Log2FC", 
       title = "Boxplot of Log2FC values for specified samples (RH4 CR TFs)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Calculate summary statistics for each specified sample
summary_stats <- EXP_filtered_samples %>%
  group_by(sample_name) %>%
  summarize(
    count = n(),
    mean = mean(log2FC),
    median = median(log2FC),
    sd = sd(log2FC),
    min = min(log2FC),
    max = max(log2FC),
    q1 = quantile(log2FC, 0.25),
    q3 = quantile(log2FC, 0.75),
    iqr = IQR(log2FC)
  )

# Print the summary statistics
print(summary_stats)











# rankplot (single sample across all coding genes), prints the list of the 20 most downregulated genes by log2FC
# can change to include multiple sameples (overlap)
# can change to include a specific gene set

# Define the sample name of interest
sample_name_of_interest <- "QL_100nM_2h"

# Set working directory and read the gene set of interest
setwd("/Users/matthewchang/Documents/R/GeneSets/")
RH4_CR_TFs <- read.table("GRYDER_RH4_CR_TFs.genelist.txt", header = FALSE, stringsAsFactors = FALSE)

# Filter for the specified sample name
EXP_single_sample <- EXP_long %>%
  filter(sample_name == sample_name_of_interest)

# Rank the Log2FC values for the single sample
EXP_ranked <- EXP_single_sample %>%
  arrange(log2FC) %>%
  mutate(rank = row_number(),
         is_target_gene = gene %in% RH4_CR_TFs$V1)

# Create the rank plot for the specified sample with highlighted target genes
ggplot(EXP_ranked, aes(x = rank, y = log2FC)) +
  geom_point(aes(color = is_target_gene), alpha = 0.5) +
  labs(x = "Rank", 
       y = "Log2FC", 
       title = paste("Rank Plot of Log2FC Values for", sample_name_of_interest, "across All Coding Genes")) +
  scale_color_manual(values = c("grey", "red"), name = "Target Gene", labels = c("Non-target", "Target")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Print the ranks and Log2FC values of the target genes
target_gene_ranks <- EXP_ranked %>%
  filter(is_target_gene) %>%
  select(gene, rank, log2FC)

print(target_gene_ranks)




# rankplot of log2FC for proteomics IHK44_100nM data

# Load data
data <- read.csv("/Users/matthewchang/Documents/R/240313_DMSOvs100nM_IHK44_stats_table.csv")

# Set working directory and read the gene set of interest
setwd("/Users/matthewchang/Documents/R/GeneSets/")
RH4_CR_TFs <- read.table("GRYDER_RH4_CR_TFs.genelist.txt", header = FALSE, stringsAsFactors = FALSE)

# Remove rows with missing values in log2_fold_change
data <- data %>%
  filter(!is.na(log2_fold_change))

# Remove "_HUMAN" suffix from Protein.Names
data <- data %>%
  mutate(Protein.Names = str_remove(Protein.Names, "_HUMAN"))

# Create a rank column based on log2_fold_change
data <- data %>%
  mutate(Rank = rank(log2_fold_change, ties.method = "first"),
         is_target_gene = Protein.Names %in% RH4_CR_TFs$V1)

# Plot the rank plot with highlighted target genes
ggplot(data, aes(x = Rank, y = log2_fold_change)) +
  geom_point(aes(color = is_target_gene), alpha = 0.5) +
  labs(title = "Rank Plot of Proteomics Data",
       x = "Rank",
       y = "Log2 Fold Change") +
  scale_color_manual(values = c("grey", "red"), name = "Target Gene", labels = c("Non-target", "Target")) +
  theme_minimal()

# Print the ranks and log2_fold_change values of the target genes
target_gene_ranks <- data %>%
  filter(is_target_gene) %>%
  select(Protein.Names, Rank, log2_fold_change)

print(target_gene_ranks)

# Display the 20 most downregulated genes
most_downregulated <- data %>%
  arrange(log2_fold_change) %>%
  head(20) %>%
  select(Protein.Names, log2_fold_change, Rank)

print(most_downregulated)
































# bar chart of multiple samples on a single gene, emphasizing effect of time based on dosage
# can change to show effect of concentration based on time
# can change to show TPM instead of log2(TPM + 1)

# Load necessary libraries
library(tidyverse)

# Create TPM df
EXP.TPM <- EXP.coding.matrix
colnames(EXP.TPM) <- gsub("^RH4_", "", colnames(EXP.TPM))
colnames(EXP.TPM) <- gsub("_RNA_022924_CWRU", "", colnames(EXP.TPM))
rownames(EXP.TPM) <- EXP.coding.matrix$gene_id

# Define the sample names of interest
sample_names_of_interest <- c("DMSO_2h", "A485_1uM_2h", "IHK44_1uM_2h", "dCBP_1uM_2h", "JQAD_1uM_2h", "LS_1uM_2h", "QL_1uM_2h",
                              "DMSO_6h", "A485_1uM_6h", "IHK44_1uM_6h", "dCBP_1uM_6h", "JQAD_1uM_6h", "LS_1uM_6h", "QL_1uM_6h")

# Define the desired order of samples
desired_order <- sample_names_of_interest

# Define the gene of interest
gene_of_interest <- "SOX8"

# Transpose the DataFrame
EXP.filtered_t <- as.data.frame(t(EXP.TPM))
EXP.filtered_t <- rownames_to_column(EXP.filtered_t, var = "sample_name")

# Gather the data to long format
EXP_long <- EXP.filtered_t %>%
  gather(key = "gene", value = "TPM", -sample_name)

# Ensure TPM values are numeric
EXP_long <- EXP_long %>%
  mutate(TPM = as.numeric(TPM))

# Filter for the specified sample names and the gene of interest
EXP_filtered_samples <- EXP_long %>%
  filter(sample_name %in% sample_names_of_interest & gene == gene_of_interest)

# Extract time point information
EXP_filtered_samples <- EXP_filtered_samples %>%
  mutate(time_point = ifelse(grepl("_2h$", sample_name), "2h", "6h"))

# Convert sample_name to factor and set the desired order
EXP_filtered_samples <- EXP_filtered_samples %>%
  mutate(sample_name = factor(sample_name, levels = desired_order))

# Calculate log2(TPM + 1)
EXP_filtered_samples <- EXP_filtered_samples %>%
  mutate(log2_TPM = log2(TPM + 1))

# Define color scheme
color_scheme <- c(
  "DMSO_2h" = "grey", "DMSO_6h" = "grey",
  "A485_1uM_2h" = "red", "A485_1uM_6h" = "red",
  "IHK44_1uM_2h" = "gold", "IHK44_1uM_6h" = "gold",
  "dCBP_1uM_2h" = "skyblue", "dCBP_1uM_6h" = "skyblue",
  "JQAD_1uM_2h" = "darkgreen", "JQAD_1uM_6h" = "darkgreen",
  "LS_1uM_2h" = "darkorchid1", "LS_1uM_6h" = "darkorchid1",
  "QL_1uM_2h" = "lemonchiffon3", "QL_1uM_6h" = "lemonchiffon3"
)

# Create the bar plot comparing the 2h and 6h changes
ggplot(EXP_filtered_samples, aes(x = sample_name, y = log2_TPM, fill = sample_name)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = color_scheme) +  # Apply color scheme
  facet_wrap(~ time_point, scales = "free_x") +  # Facet by time point
  labs(x = "Sample Name", 
       y = "log2(TPM + 1)", 
       title = paste("log2(TPM + 1) for", gene_of_interest, "across Samples")) +
  theme_minimal() +
  theme(axis.text.x = element_blank())


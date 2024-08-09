### Updated: 08.09.2024 (Matt)

### 1. More heatmaps for RNA-seq 
# note: looking by gene separated by timepoint, conc and drug tx

### Functions
# Extracts drug, dosage, and timepoint information from sample name
extract_info = function(name) {
  parts = strsplit(name, "_")[[1]]     # parse through sample name by "_"
  drug = parts[2]
  dosage = ""
  timepoint = parts[3]
  
  # check if timepoint is mislabeled
  if (grepl("uM|nM", parts[3])) {
    dosage = parts[3]
    timepoint = parts[4]
  }
  
  return(list(drug = drug, dosage = dosage, timepoint = timepoint))
}

### Setup: load libraries, get TPM matrix and gene list (and subset samples optionally)
library(pheatmap)
library(ggplot2)
library(reshape2)
library(gridExtra)

# EDIT HERE
gene = "PAX3"
desired.drug.order = c("JQAD", "LS", "QL", "dCBP", "A485", "IHK44")

setwd("/Volumes/rc/SOM_GENE_BEG33/RNA_seq/hg38/projects/RMS_IHK/Practice_MSC/ExpMatrices/")
EXP.coding.matrix = read.table("IHK_samples.coding.norm.matrix.txt", header = T)

# Filtering for gene
EXP.filtered.matrix = EXP.coding.matrix[EXP.coding.matrix$gene_id == gene, ]
rownames(EXP.filtered.matrix) = EXP.filtered.matrix$gene_id
EXP.filtered.matrix = EXP.filtered.matrix[, -1]

# Disinclude DMSO and NT (first 4 columns) from drug heatmap
dmso.nt.heatmap = EXP.filtered.matrix[, 1:4]
rownames(dmso.nt.heatmap) = EXP.filtered.matrix$gene_id
EXP.filtered.matrix = EXP.filtered.matrix[, -(1:4)]

# Log2 transformation: log2(TPM + 1)
dmso.nt.heatmap_log2 = log2(dmso.nt.heatmap + 1)
EXP.filtered.matrix_log2 = log2(EXP.filtered.matrix + 1)

# Create lists to store the sample's drug, dosage, and timepoint information
drugs = character(length(colnames(EXP.filtered.matrix_log2)))
dosages = character(length(colnames(EXP.filtered.matrix_log2)))
timepoints = character(length(colnames(EXP.filtered.matrix_log2)))
for (j in seq_along(colnames(EXP.filtered.matrix_log2))) {
  info = extract_info(colnames(EXP.filtered.matrix_log2)[j])
  drugs[j] = info$drug
  dosages[j] = info$dosage
  timepoints[j] = info$timepoint
}
unique.drugs = unique(drugs)
unique.timepoints = unique(timepoints)
unique.dosages = unique(dosages)
unique.drugs = intersect(desired.drug.order, unique.drugs)

heatmap.data = matrix(NA, nrow = length(unique.dosages) * length(unique.timepoints), ncol = length(unique.drugs))
rownames(heatmap.data) = paste(rep(unique.dosages, each = length(unique.timepoints)), unique.timepoints, sep = "_")
colnames(heatmap.data) = unique.drugs
for (j in seq_along(colnames(EXP.filtered.matrix_log2))) {
  drug = drugs[j]
  dosage = dosages[j]
  timepoint = timepoints[j]
  
  if (paste(dosage, timepoint, sep = "_") %in% rownames(heatmap.data) && drug %in% colnames(heatmap.data)) {
    heatmap.data[paste(dosage, timepoint, sep = "_"), drug] = EXP.filtered.matrix_log2[1, j]
  }
}
heatmap.data = as.data.frame(heatmap.data)
heatmap.data <- apply(heatmap.data, 2, function(x) as.numeric(x))
rownames(heatmap.data) = paste(rep(unique.dosages, each = length(unique.timepoints)), unique.timepoints, sep = "_")
heatmap.data = heatmap.data[, colSums(is.na(heatmap.data)) == 0]

colnames(dmso.nt.heatmap_log2) = c("DMSO_2h", "DMSO_6h", "NT_2h", "NT_6h")

combined_data_log2 = cbind(dmso.nt.heatmap_log2, EXP.filtered.matrix_log2)
local.min = min(combined_data_log2, na.rm = T)
local.max = max(combined_data_log2, na.rm = T)

# Plot
p1 = pheatmap(t(dmso.nt.heatmap_log2), 
              cluster_rows = F, 
              cluster_cols = F, 
              scale = "none", 
              show_rownames = T, 
              show_colnames = T,
              breaks = seq(local.min, local.max, length.out = 101),
              main = paste(gene),
              legend = F)
p2 = pheatmap(heatmap.data, 
              cluster_rows = F, 
              cluster_cols = F, 
              scale = "none", 
              show_rownames = T, 
              show_colnames = T,
              angle_col = 315,
              breaks = seq(local.min, local.max, length.out = 101), 
              main = paste(gene))
grid.arrange(p1$gtable, p2$gtable, ncol = 2, widths = c(0.75, 2.5))

# figure out color scale (show controls do not reflect in color change...)
################ end of custom heatmaps ################












### 2. box/violin plot by log2FC by condition separated by multiple genesets
# note: all coding genes did not show me anything good

library(tidyverse)

setwd("/Volumes/rc/SOM_GENE_BEG33/RNA_seq/hg38/projects/RMS_IHK/Practice_MSC/ExpMatrices")
EXP.log2FC = read.table("EXP.log2FC.txt", header = T, row.names = 1)
colnames(EXP.log2FC) = gsub("^RH4_", "", colnames(EXP.log2FC))

# EDIT HERE
samples =        c("A485_1uM_6h", 
                   "IHK44_1uM_6h", 
                   "dCBP_1uM_6h")
samples.colors = c("A485_1uM_6h" = "red", 
                   "IHK44_1uM_6h" = "darkorange", 
                   "dCBP_1uM_6h" = "green3")

# EDIT HERE
setwd("/Volumes/rc/SOM_GENE_BEG33/RNA_seq/hg38/projects/RMS_IHK/Practice_MSC/GeneSets/")
housekeeping_genes = read.table("house_keeping_genes.txt", header = F, stringsAsFactors = F) %>%
  mutate(gene_set = "Housekeeping genes")
P3F_target_genes = read.table("GRYDER_PAX3FOXO1_ENHANCERS_KO_DOWN.txt", header = F, stringsAsFactors = F) %>%
  mutate(gene_set = "P3F targets genes")
RH4_CR_TFs = read.table("GRYDER_RH4_CR_TFs.genelist.txt", header = F, stringsAsFactors = F) %>%
  mutate(gene_set = "RH4 CR TFs")
all_genesets = bind_rows(housekeeping_genes, P3F_target_genes, RH4_CR_TFs)

# Subset for genes and samples
EXP.log2FC.filter = EXP.log2FC[rownames(EXP.log2FC) %in% all_genesets$V1, ]                                     # subset genes
EXP.log2FC.filter = EXP.log2FC.filter[, colnames(EXP.log2FC.filter) %in% samples]                               # subset samples
EXP.log2FC.filter = EXP.log2FC.filter[, samples]                                                                # set the order of the samples

# Transform for plotting
EXP.log2FC.plot = as.data.frame(EXP.log2FC.filter)                                                              # copy over filtered df
EXP.log2FC.plot = rownames_to_column(EXP.log2FC.plot, var = "gene")                                             # transfer rownames (genes) to a column
EXP.log2FC.plot = pivot_longer(EXP.log2FC.plot, cols = -gene, names_to = "sample_name", values_to = "log2FC")   # transform to long format (gene, sample_name, log2FC)
EXP.log2FC.plot = left_join(EXP.log2FC.plot, all_genesets, by = c("gene" = "V1"))
EXP.log2FC.plot$sample_name = factor(EXP.log2FC.plot$sample_name, levels = samples)

# Plot
ggplot(EXP.log2FC.plot, aes(x = sample_name, y = log2FC, fill = sample_name)) +
  geom_violin(alpha = 0.5, scale = "width", width = 0.5) + 
  geom_boxplot(width = 0.1, color = "black", alpha = 0.7, outlier.shape = NA) +
  facet_wrap(~ gene_set, scales = "fixed") +
  scale_fill_manual(values = samples.colors) +
  labs(x = "Sample Name", 
       y = "Log2FC", 
       title = "Boxplot of Log2FC for 1uM 6h samples") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


################ end of box/violin plots ################ 











### 3. bar plots for a single gene

setwd("/Volumes/rc/SOM_GENE_BEG33/RNA_seq/hg38/projects/RMS_IHK/Practice_MSC/ExpMatrices/")
EXP.coding.matrix = read.table("IHK_samples.coding.norm.matrix.txt", header = T)
cols.to.keep = 1:(ncol(EXP.coding.matrix) - 4)
EXP.coding.matrix = EXP.coding.matrix[, cols.to.keep]

# EDIT HERE
gene = "MYOD1"
samples =        c("NT_6h",
                   "DMSO_6h",
                   "JQAD_100nM_6h",
                   "LS_100nM_6h",
                   "QL_100nM_6h",
                   "dCBP_100nM_6h",
                   "A485_100nM_6h", 
                   "IHK44_100nM_6h")
samples.colors = c("NT_6h" = "gray80",
                   "DMSO_6h" = "gray30",
                   "JQAD_100nM_6h" = "yellow",
                   "LS_100nM_6h" = "lightblue",
                   "QL_100nM_6h" = "violet",
                   "dCBP_100nM_6h" = "green3",
                   "A485_100nM_6h" = "red", 
                   "IHK44_100nM_6h" = "darkorange")

# Filter for gene and subset samples
EXP.single.gene = EXP.coding.matrix[EXP.coding.matrix$gene_id == gene, ]
colnames(EXP.single.gene) = gsub("^RH4_", "", colnames(EXP.single.gene))
colnames(EXP.single.gene) = gsub("_RNA_022924_CWRU", "", colnames(EXP.single.gene))
EXP.single.gene = EXP.single.gene[, c("gene_id", samples)]
EXP.single.gene = pivot_longer(EXP.single.gene, cols = -gene_id, names_to = "sample_name", values_to = "expression")

# Plot
ggplot(EXP.single.gene, aes(x = factor(sample_name, levels = samples), y = expression, fill = sample_name)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) + 
  scale_fill_manual(values = samples.colors) +
  labs(x = "sample", 
       y = "expression level", 
       title = paste("bar plot of gene expression for", gene)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

################ end of bar plots ################ 












### 4. rankplot (single sample across all coding genes), prints the list of the 20 most downregulated genes by log2FC
# can change to include multiple sameples (overlap)
# can change to include a specific gene set

library(tidyverse)

# Define sample
sample = "QL_100nM_2h"

# Read gene set to locate rank
setwd("/Volumes/rc/SOM_GENE_BEG33/RNA_seq/hg38/projects/RMS_IHK/Practice_MSC/GeneSets/")
geneset = read.table("GRYDER_RH4_CR_TFs.genelist.txt", header = F, stringsAsFactors = F)
geneset.name = "GRYDER_RH4_CR_TFs.genelist.txt"

# Create log2FC df
setwd("/Volumes/rc/SOM_GENE_BEG33/RNA_seq/hg38/projects/RMS_IHK/Practice_MSC/GSEA_log2FC/p300_samples/")
EXP.log2FC = read.table("EXP.log2FC.txt", header = T, row.names = 1)                                            # read in log2FC values
colnames(EXP.log2FC) = gsub("^RH4_", "", colnames(EXP.log2FC))                                                  # remove "RH4_" from column names

# Filtering
EXP.log2FC.sample = EXP.log2FC[, sample, drop = F]                                                              # subset for the sample

# Ranking
EXP.log2FC.df <- data.frame(gene = rownames(EXP.log2FC.sample), log2FC = as.numeric(EXP.log2FC.sample[, 1]))    # create new df with new columns (gene, log2FC)
EXP.log2FC.df <- EXP.log2FC.df[order(EXP.log2FC.df$log2FC), ]                                                   # orders the rows by log2FC values in ascending order
EXP.log2FC.df$rank <- seq_len(nrow(EXP.log2FC.df))                                                              # ranks the rows in ascending order in a new column (rank)
EXP.log2FC.df$target <- EXP.log2FC.df$gene %in% geneset$V1                                                      # add a new column (target) and checks if gene is in geneset

# Plot
ggplot(EXP.log2FC.df, aes(x = rank, y = log2FC)) +
  geom_point(aes(color = target), alpha = 0.5) +
  labs(x = "Rank", 
       y = "Log2FC", 
       title = paste("Rank plot for", sample, "across all coding genes"),
       subtitle = paste("Target genes are for", geneset.name)) +
  scale_color_manual(values = c("grey", "red"), name = "Target Gene", labels = c("Non-target", "Target")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

################ end of rank plots ################

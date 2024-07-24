### Author: Matthew Chang
### More heatmaps for RNA-seq (looking by gene separated by timepoint, conc and drug tx)

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

### Setup: import libraries, set working directories, get TPM matrix, and gene list (and subset samples optionally)
library(pheatmap)
library(ggplot2)
library(reshape2)
library(gridExtra)

setwd("/Users/matthewchang/Documents/R/ExpMatrices/")
EXP.coding.matrix = read.table("IHK_samples.coding.norm.matrix.txt", header = T)
EXP.select.matrix = EXP.coding.matrix    # copy over expression matrix

setwd("/Users/matthewchang/Documents/R/GeneSets/")
gene.list = read.table("MSC_TF_list.txt", header = F, stringsAsFactors = F)
genes = gene.list$V1

### OPTIONAL: Subset samples
samples.to.keep <- c("RH4_DMSO_2h_RNA_022924_CWRU", "RH4_DMSO_6h_RNA_022924_CWRU",
                     "RH4_IHK44_1uM_6h_RNA_022924_CWRU", "RH4_IHK44_100nM_6h_RNA_022924_CWRU",
                     "RH4_IHK44_1uM_2h_RNA_022924_CWRU", "RH4_IHK44_100nM_2h_RNA_022924_CWRU")
sample.indices = which(colnames(EXP.coding.matrix) %in% c("gene_id", samples.to.keep))
EXP.select.matrix = EXP.coding.matrix[, sample.indices]

### Define the desired drug order
desired.drug.order = c("A485", "IHK44", "dCBP", "JQAD", "QL", "LS") # Modify this vector to your desired order

### Loop through each gene and generate a heatmap of the log2(TPM + 1) expression values
for (gene in genes) {
  # Filter TPM matrix to only include the gene of interest
  EXP.filtered.matrix = EXP.select.matrix[EXP.select.matrix$gene_id == gene, ]
  
  # Safety check: if gene not found, then skip to the next gene
  if (nrow(EXP.filtered.matrix) == 0) {
    message("Gene not found: ", gene)
    next
  }
  
  # Rename row to gene id and remove the gene id column
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
  
  # Fill out lists with sample information, calling 'extract_info' function defined above
  for (j in seq_along(colnames(EXP.filtered.matrix_log2))) {
    info = extract_info(colnames(EXP.filtered.matrix_log2)[j])
    drugs[j] = info$drug
    dosages[j] = info$dosage
    timepoints[j] = info$timepoint
  }
  
  # Safety check: remove repeats if any
  unique.drugs = unique(drugs)
  unique.timepoints = unique(timepoints)
  unique.dosages = unique(dosages)
  
  # Ensure the heatmap follows the desired drug order
  unique.drugs = intersect(desired.drug.order, unique.drugs)
  
  # Initialize a matrix to store the heatmap data with NA values, and rename the rows and columns
  heatmap.data = matrix(NA, nrow = length(unique.dosages) * length(unique.timepoints), ncol = length(unique.drugs))
  rownames(heatmap.data) = paste(rep(unique.dosages, each = length(unique.timepoints)), unique.timepoints, sep = "_")
  colnames(heatmap.data) = unique.drugs

  
  # Fill the heatmap data
  for (j in seq_along(colnames(EXP.filtered.matrix_log2))) {
    drug = drugs[j]
    dosage = dosages[j]
    timepoint = timepoints[j]
    
    if (paste(dosage, timepoint, sep = "_") %in% rownames(heatmap.data) && drug %in% colnames(heatmap.data)) {
      heatmap.data[paste(dosage, timepoint, sep = "_"), drug] = EXP.filtered.matrix_log2[1, j]
    }
  }
  
  # Transform matrix to df for heatmap
  heatmap.data = as.data.frame(heatmap.data)
  
  # Convert logical values to numeric before any operations that require numeric input, rename rownames
  heatmap.data <- apply(heatmap.data, 2, function(x) as.numeric(x))
  rownames(heatmap.data) = paste(rep(unique.dosages, each = length(unique.timepoints)), unique.timepoints, sep = "_")
  
  # Safety check: remove columns with any NA values
  heatmap.data = heatmap.data[, colSums(is.na(heatmap.data)) == 0]
  
  # Check if heatmap.data is empty after filtering NAs
  if (ncol(heatmap.data) == 0 || nrow(heatmap.data) == 0) {
    message("No valid data to plot for gene: ", gene)
    next
  }
  
  # Update column names to format <DMSO_2h>
  colnames(dmso.nt.heatmap_log2) <- c("DMSO_2h", "DMSO_6h", "NT_2h", "NT_6h")
  
  # Calculate the local min and max for the current gene
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
                angle_col = 45,
                breaks = seq(local.min, local.max, length.out = 101),
                main = paste(gene),
                legend = F)
  
  p2 = pheatmap(heatmap.data, 
                cluster_rows = F, 
                cluster_cols = F, 
                scale = "none", 
                show_rownames = T, 
                show_colnames = T,
                angle_col = 45,
                breaks = seq(local.min, local.max, length.out = 101),  # Fixed log2 color scale
                main = paste(gene))
  
  grid.arrange(p1$gtable, p2$gtable, ncol = 2)
}

################ end of custom heatmaps ################



# boxplot/violin plot by log2FC by condition separated by multiple genesets
# note: all coding genes did not show me anything good

library(tidyverse)

# Define samples to work on
samples = c("A485_1uM_2h", "dCBP_1uM_2h", "IHK44_1uM_2h", "JQAD_1uM_2h", "LS_1uM_2h", "QL_1uM_2h")

# Read gene set
setwd("/Users/matthewchang/Documents/R/GeneSets/")
geneset = read.table("GRYDER_RH4_CR_TFs.genelist.txt", header = F, stringsAsFactors = F)

# Create log2FC df
setwd("/Users/matthewchang/Documents/R/GSEA_log2FC/p300_samples/")
EXP.log2FC = read.table("EXP.log2FC.txt", header = T, row.names = 1)                                            # read in log2FC values
colnames(EXP.log2FC) = gsub("^RH4_", "", colnames(EXP.log2FC))                                                  # remove "RH4_" from column names

# Filtering steps
EXP.log2FC.filter = EXP.log2FC[rownames(EXP.log2FC) %in% geneset$V1, ]                                          # subset genes
EXP.log2FC.filter = EXP.log2FC.filter[, colnames(EXP.log2FC.filter) %in% samples]                               # subset samples
EXP.log2FC.filter = EXP.log2FC.filter[, samples]                                                                # set the order of the samples

# Transforming steps
EXP.log2FC.plot = as.data.frame(EXP.log2FC.filter)                                                              # copy over filtered df
EXP.log2FC.plot = rownames_to_column(EXP.log2FC.plot, var = "gene")                                             # transfer rownames (genes) to a column
EXP.log2FC.plot = pivot_longer(EXP.log2FC.plot, cols = -gene, names_to = "sample_name", values_to = "log2FC")   # transform to long format (gene, sample_name, log2FC)

# Plot
ggplot(EXP.log2FC.plot, aes(x = sample_name, y = log2FC, fill = sample_name)) +
  geom_violin(alpha = 0.5) + 
  geom_boxplot(width = 0.1, color = "black", alpha = 0.7, outlier.shape = NA) +  # Boxplot without outliers
  labs(x = "Sample Name", 
       y = "Log2FC", 
       title = "Boxplot of Log2FC values for specified samples (RH4 CR TFs)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


################ end of box / violin plots ################ 











# rankplot (single sample across all coding genes), prints the list of the 20 most downregulated genes by log2FC
# can change to include multiple sameples (overlap)
# can change to include a specific gene set

library(tidyverse)

# Define sample
sample = "QL_100nM_2h"

# Read gene set to locate rank
setwd("/Users/matthewchang/Documents/R/GeneSets/")
geneset = read.table("GRYDER_RH4_CR_TFs.genelist.txt", header = F, stringsAsFactors = F)
geneset.name = "GRYDER_RH4_CR_TFs.genelist.txt"

# Create log2FC df
setwd("/Users/matthewchang/Documents/R/GSEA_log2FC/p300_samples/")
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


######## normalize by gene, separated by drug treatment and timepoint, with a fixed color gradient ########

# Author: Matthew Chang
library(pheatmap)

### Setup: get TPM matrix and gene list
setwd("/Volumes/rc/SOM_GENE_BEG33/RNA_seq/hg38/projects/RMS_IHK/Practice_MSC/ExpMatrices")
EXP.coding.matrix = read.table(file.choose(), header = T)

setwd("/Volumes/rc/SOM_GENE_BEG33/RNA_seq/hg38/projects/Genesets/List_format/")
gene.list = read.table(file.choose(), header = F, stringsAsFactors = F)  # selected GRYDER_RH4_CR_TFs_CRISPRTop.genelist.txt
genes = gene.list$V1
# Create temporary array to hold TPM matrix of selected samples (if any)
EXP.select.matrix = EXP.coding.matrix

### OPTIONAL: Subset samples
# Define a vector of sample names to keep
samples.to.keep <- c("RH4_IHK44_1uM_6h_RNA_022924_CWRU", "RH4_IHK44_100nM_6h_RNA_022924_CWRU",
                     "RH4_IHK44_1uM_2h_RNA_022924_CWRU", "RH4_IHK44_100nM_2h_RNA_022924_CWRU")

# Ensure the gene_id column is included in the subset
sample.indices = which(colnames(EXP.coding.matrix) %in% c("gene_id", samples.to.keep))
EXP.select.matrix = EXP.coding.matrix[, sample.indices]

### Functions
# Extracts drug, dosage, and timepoint information from sample name
extract_info <- function(name) {
  parts <- strsplit(name, "_")[[1]]
  drug <- parts[2]
  dosage <- ""
  timepoint <- parts[3]
  
  # check if timepoint is mislabeled
  if (grepl("uM|nM", parts[3])) {
    dosage <- parts[3]
    timepoint <- parts[4]
  }
  
  return(list(drug = drug, dosage = dosage, timepoint = timepoint))
}




# Loop through each gene and generates a heatmap
for (gene in genes) {
  # Filter TPM matrix to only include the gene of interest
  EXP.filtered.matrix = EXP.select.matrix[EXP.select.matrix$gene_id == gene, ]
  
  # Safety check: if gene not found, then skip to the next gene
  if (nrow(EXP.filtered.matrix) == 0) {
    next
  }
  
  # Rename row to gene id and remove the gene id column
  rownames(EXP.filtered.matrix) = EXP.filtered.matrix$gene_id
  EXP.filtered.matrix = EXP.filtered.matrix[, -1]
  
  # Normalize TPM values by the maximum value in each row (gene)
  max.TPM = apply(EXP.filtered.matrix, 1, max)
  EXP.filtered.matrix_normalized = sweep(EXP.filtered.matrix, 1, max.TPM, FUN = "/")
  
  # Create lists for the sample's drug, dosage, and timepoint information
  drugs = character(length(colnames(EXP.filtered.matrix_normalized)))
  dosages = character(length(colnames(EXP.filtered.matrix_normalized)))
  timepoints = character(length(colnames(EXP.filtered.matrix_normalized)))
  
  # Fill out lists
  for (j in seq_along(colnames(EXP.filtered.matrix_normalized))) {
    info = extract_info(colnames(EXP.filtered.matrix_normalized)[j])
    drugs[j] = info$drug
    dosages[j] = info$dosage
    timepoints[j] = info$timepoint
  }
  
  # Safety check: remove redundancies
  unique.drugs = unique(drugs)
  unique.timepoints = unique(timepoints)
  unique.dosages = unique(dosages)
  
  # Initialize the matrix with NA values
  heatmap.data = matrix(NA, nrow = length(unique.timepoints), ncol = length(unique.drugs) * length(unique.dosages))
  rownames(heatmap.data) = unique.timepoints
  colnames(heatmap.data) = paste(rep(unique.drugs, each = length(unique.dosages)), unique.dosages, sep = "_")
  
  # Fill the heatmap data
  for (j in seq_along(colnames(EXP.filtered.matrix_normalized))) {
    drug = drugs[j]
    dosage = dosages[j]
    timepoint = timepoints[j]
    heatmap.data[timepoint, paste(drug, dosage, sep = "_")] = EXP.filtered.matrix_normalized[1, j]
  }
  
  # Transform matrix to df for heatmap
  heatmap.data = as.data.frame(heatmap.data)
  
  # Safety check: remove columns with any NA values
  heatmap.data = heatmap.data[, colSums(is.na(heatmap.data)) == 0]
  
  # Plot heatmap with a fixed color scale from 0 to 1
  pheatmap(heatmap.data, 
           cluster_rows = F, 
           cluster_cols = F, 
           scale = "none", 
           show_rownames = T, 
           show_colnames = T,
           angle_col = 45,
           breaks = seq(0, 1, length.out = 101),  # Fixed color scale from 0 to 1
           main = paste("Normalization Expression for", gene))
}




######## log2(TPM + 1) calculations ########

# Author: Matthew Chang
library(pheatmap)

### Setup: get TPM matrix and gene list
setwd("/Volumes/rc/SOM_GENE_BEG33/RNA_seq/hg38/projects/RMS_IHK/Practice_MSC/ExpMatrices")
EXP.coding.matrix = read.table(file.choose(), header = T)

setwd("/Volumes/rc/SOM_GENE_BEG33/RNA_seq/hg38/projects/Genesets/List_format/")
gene.list = read.table(file.choose(), header = F, stringsAsFactors = F)  # selected GRYDER_RH4_CR_TFs_CRISPRTop.genelist.txt
genes = gene.list$V1

### Functions
# Extracts drug, dosage, and timepoint information from sample name
extract_info <- function(name) {
  parts <- strsplit(name, "_")[[1]]
  drug <- parts[2]
  dosage <- ""
  timepoint <- parts[3]
  
  # check if timepoint is mislabeled
  if (grepl("uM|nM", parts[3])) {
    dosage <- parts[3]
    timepoint <- parts[4]
  }
  
  return(list(drug = drug, dosage = dosage, timepoint = timepoint))
}

# Calculate log2(TPM + 1) for the entire matrix to find global min and max
EXP.coding.matrix.log2 = log2(EXP.coding.matrix[, -1] + 1)
global.min = min(EXP.coding.matrix.log2, na.rm = T)
global.max = max(EXP.coding.matrix.log2, na.rm = T)

# Loop through each gene and generate a heatmap
for (gene in genes) {
  # Filter TPM matrix to only include the gene of interest
  EXP.filtered.matrix = EXP.coding.matrix[EXP.coding.matrix$gene_id == gene, ]
  
  # Safety check: if gene not found, then skip to the next gene
  if (nrow(EXP.filtered.matrix) == 0) {
    next
  }
  
  # Rename row to gene id and remove the gene id column
  rownames(EXP.filtered.matrix) = EXP.filtered.matrix$gene_id
  EXP.filtered.matrix = EXP.filtered.matrix[, -1]
  
  # Log2 transformation: log2(TPM + 1)
  EXP.filtered.matrix_log2 = log2(EXP.filtered.matrix + 1)
  
  # Create lists for the sample's drug, dosage, and timepoint information
  drugs = character(length(colnames(EXP.filtered.matrix_log2)))
  dosages = character(length(colnames(EXP.filtered.matrix_log2)))
  timepoints = character(length(colnames(EXP.filtered.matrix_log2)))
  
  # Fill out lists
  for (j in seq_along(colnames(EXP.filtered.matrix_log2))) {
    info = extract_info(colnames(EXP.filtered.matrix_log2)[j])
    drugs[j] = info$drug
    dosages[j] = info$dosage
    timepoints[j] = info$timepoint
  }
  
  # Safety check: remove redundancies
  unique.drugs = unique(drugs)
  unique.timepoints = unique(timepoints)
  unique.dosages = unique(dosages)
  
  # Initialize the matrix with NA values
  heatmap.data = matrix(NA, nrow = length(unique.timepoints), ncol = length(unique.drugs) * length(unique.dosages))
  rownames(heatmap.data) = unique.timepoints
  colnames(heatmap.data) = paste(rep(unique.drugs, each = length(unique.dosages)), unique.dosages, sep = "_")
  
  # Fill the heatmap data
  for (j in seq_along(colnames(EXP.filtered.matrix_log2))) {
    drug = drugs[j]
    dosage = dosages[j]
    timepoint = timepoints[j]
    heatmap.data[timepoint, paste(drug, dosage, sep = "_")] = EXP.filtered.matrix_log2[1, j]
  }
  
  # Transform matrix to df for heatmap
  heatmap.data = as.data.frame(heatmap.data)
  
  # Safety check: remove columns with any NA values
  heatmap.data = heatmap.data[, colSums(is.na(heatmap.data)) == 0]
  
  # Plot heatmap with a fixed log2 color scale based on global min and max
  pheatmap(heatmap.data, 
           cluster_rows = F, 
           cluster_cols = F, 
           scale = "none", 
           show_rownames = T, 
           show_colnames = T,
           angle_col = 45,
           breaks = seq(global.min, global.max, length.out = 101),  # Fixed log2 color scale
           main = paste("Log2 Normalized Expression for", gene))
}


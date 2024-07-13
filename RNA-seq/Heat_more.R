### Author: Matthew Chang
### More heatmaps and bar charts for RNA-seq gene expression (normalization, log2(TPM + 1), z-score)

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

setwd("/Users/matthewchang/Documents/R/")
EXP.coding.matrix = read.table(file.choose(), header = T)
EXP.select.matrix = EXP.coding.matrix    # copy over expression matrix

setwd("/Users/matthewchang/Documents/R/")
gene.list = read.table(file.choose(), header = F, stringsAsFactors = F)     # selected GRYDER_RH4_CR_TFs_CRISPRTop.genelist.txt
genes = gene.list$V1

### OPTIONAL: Subset samples
samples.to.keep <- c("RH4_IHK44_1uM_6h_RNA_022924_CWRU", "RH4_IHK44_100nM_6h_RNA_022924_CWRU",
                     "RH4_IHK44_1uM_2h_RNA_022924_CWRU", "RH4_IHK44_100nM_2h_RNA_022924_CWRU",
                     "RH4_dCBP_1uM_6h_RNA_022924_CWRU", "RH4_dCBP_100nM_6h_RNA_022924_CWRU",
                     "RH4_dCBP_1uM_2h_RNA_022924_CWRU", "RH4_dCBP_100nM_2h_RNA_022924_CWRU")
sample.indices = which(colnames(EXP.coding.matrix) %in% c("gene_id", samples.to.keep))
EXP.select.matrix = EXP.coding.matrix[, sample.indices]







######## normalize by gene, separated by drug treatment and timepoint, with a fixed color gradient across all samples ########
### Functions
# Creates heatmap
create_heatmap_normalized = function(data, gene) {
  pheatmap(data,
           cluster_rows = F, 
           cluster_cols = F, 
           scale = "none", 
           show_rownames = T, 
           show_colnames = T,
           angle_col = 45,
           breaks = seq(0, 1, length.out = 101), # Fixed normalized expression gradient
           main = paste("Normalization Expression for", gene))
}

# Creates bar chart
create_bar_chart_normalized <- function(data, gene) {
  data <- cbind(Timepoint = rownames(data), data)
  data_melt <- reshape2::melt(data, id.vars = "Timepoint")
  colnames(data_melt) <- c("Timepoint", "Condition", "Value")
  ggplot(data_melt, aes(x = Timepoint, y = Value, fill = Condition)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste("Normalized Expression Levels for", gene), y = "Normalized Value") +
    ylim(0, 1)  # Set y-axis limits from 0 to 1
}

### Loop through each gene and generates a heatmap of the normalized expression values
for (gene in genes) {
  # Filter TPM matrix to only include the gene of interest
  EXP.filtered.matrix = EXP.select.matrix[EXP.select.matrix$gene_id == gene, ]
  
  # Safety check: if gene not found, then skip to the next gene
  if (nrow(EXP.filtered.matrix) == 0) {
    next
  }
  
  # Rename rows to gene id and remove the gene id column
  rownames(EXP.filtered.matrix) = EXP.filtered.matrix$gene_id
  EXP.filtered.matrix = EXP.filtered.matrix[, -1]
  
  # Normalize TPM values by the maximum value in each row (gene)
  max.TPM = apply(EXP.filtered.matrix, 1, max)
  EXP.filtered.matrix_normalized = sweep(EXP.filtered.matrix, 1, max.TPM, FUN = "/")
  
  # Create lists to store the sample's drug, dosage, and timepoint information
  drugs = character(length(colnames(EXP.filtered.matrix_normalized)))
  dosages = character(length(colnames(EXP.filtered.matrix_normalized)))
  timepoints = character(length(colnames(EXP.filtered.matrix_normalized)))
  
  # Fill out lists with sample information, calling 'extract_info' function defined above
  for (j in seq_along(colnames(EXP.filtered.matrix_normalized))) {
    info = extract_info(colnames(EXP.filtered.matrix_normalized)[j])
    drugs[j] = info$drug
    dosages[j] = info$dosage
    timepoints[j] = info$timepoint
  }
  
  # Safety check: remove repeats if any
  unique.drugs = unique(drugs)
  unique.timepoints = unique(timepoints)
  unique.dosages = unique(dosages)
  
  # Initialize a matrix to store the heatmap data with NA values, and rename the rows and columns
  heatmap.data = matrix(NA, nrow = length(unique.timepoints), ncol = length(unique.drugs) * length(unique.dosages))
  rownames(heatmap.data) = unique.timepoints
  colnames(heatmap.data) = paste(rep(unique.drugs, each = length(unique.dosages)), unique.dosages, sep = "_")
  
  # Fill in the heatmap data matrix from the lists with the normalized expression value
  for (j in seq_along(colnames(EXP.filtered.matrix_normalized))) {
    drug = drugs[j]
    dosage = dosages[j]
    timepoint = timepoints[j]
    heatmap.data[timepoint, paste(drug, dosage, sep = "_")] = EXP.filtered.matrix_normalized[1, j]
  }
  
  # Transform matrix to df
  heatmap.data = as.data.frame(heatmap.data)
  
  # Safety check: remove columns with any NA values
  heatmap.data = heatmap.data[, colSums(is.na(heatmap.data)) == 0]
  
  # Plot
  # create_heatmap_normalized(heatmap.data, gene)
  print(create_bar_chart_normalized(heatmap.data, gene))
}









######## log2(TPM + 1) by gene, separated by drug treatment and timepoint, with a fixed color gradient across all samples ########
### Functions
# Creates heatmap
create_heatmap_log2 = function(data, gene) {
  pheatmap(heatmap.data, 
           cluster_rows = F, 
           cluster_cols = F, 
           scale = "none", 
           show_rownames = T, 
           show_colnames = T,
           angle_col = 45,
           breaks = seq(global.min, global.max, length.out = 101),  # Fixed log2 color scale
           main = paste("Log2(TPM + 1) for", gene))
}

# Creates bar chart
create_bar_chart_log2 <- function(data, gene) {
  data <- cbind(Timepoint = rownames(data), data)
  data_melt <- reshape2::melt(data, id.vars = "Timepoint")
  colnames(data_melt) <- c("Timepoint", "Condition", "Value")
  ggplot(data_melt, aes(x = Timepoint, y = Value, fill = Condition)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste("Log2(TPM + 1) Expression Levels for", gene), y = "Log2(TPM + 1)") +
    ylim(global.min, global.max)  # Adjust y-axis limits according to your data range
}

### Calculate log2(TPM + 1) for the entire matrix to find global min and max
EXP.select.matrix.log2 = log2(EXP.select.matrix[, -1] + 1)
global.min = min(EXP.select.matrix.log2, na.rm = T)
global.max = max(EXP.select.matrix.log2, na.rm = T)

### Loop through each gene and generate a heatmap of the log2(TPM + 1) expression values
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
  
  # Log2 transformation: log2(TPM + 1)
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
  
  # Initialize a matrix to store the heatmap data with NA values, and rename the rows and columns
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
  
  # Plot
  # create_heatmap_log2(heatmap.data, gene)
  print(create_bar_chart_log2(heatmap.data, gene))
}





######## z-score by gene, separated by drug treatment and timepoint, with a fixed color gradient across all samples ########
### Functions
# Transform expression data to z-score
z_score_transform <- function(x) {
  (x - mean(x, na.rm = T)) / sd(x, na.rm = T)
}

# Creates heatmap
create_heatmap_z <- function(data, gene) {
  pheatmap(data, 
           cluster_rows = F, 
           cluster_cols = F, 
           scale = "none", 
           show_rownames = T, 
           show_colnames = T,
           angle_col = 45,
           breaks = seq(-3, 3, length.out = 101),  # Z-score color scale
           main = paste("Z-score Transformation for", gene))
}

# Creates bar chart
create_bar_chart_z <- function(data, gene) {
  data <- cbind(Timepoint = rownames(data), data)
  data_melt <- reshape2::melt(data, id.vars = "Timepoint")
  colnames(data_melt) <- c("Timepoint", "Condition", "Value")
  ggplot(data_melt, aes(x = Timepoint, y = Value, fill = Condition)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste("Z-score Expression Levels for", gene), y = "Z-score") +
    ylim(-3, 3)  # Set y-axis limits from -3 to 3 for z-score
}

### Loop through each gene and generate a heatmap of the log2(TPM + 1) expression values
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
  
  # Log2 transformation: log2(TPM + 1)
  EXP.filtered.matrix.z = t(apply(EXP.filtered.matrix, 1, z_score_transform))
  
  # Create lists to store the sample's drug, dosage, and timepoint information
  drugs = character(length(colnames(EXP.filtered.matrix.z)))
  dosages = character(length(colnames(EXP.filtered.matrix.z)))
  timepoints = character(length(colnames(EXP.filtered.matrix.z)))
  
  # Fill out lists with sample information, calling 'extract_info' function defined above
  for (j in seq_along(colnames(EXP.filtered.matrix.z))) {
    info = extract_info(colnames(EXP.filtered.matrix.z)[j])
    drugs[j] = info$drug
    dosages[j] = info$dosage
    timepoints[j] = info$timepoint
  }
  
  # Safety check: remove repeats if any
  unique.drugs = unique(drugs)
  unique.timepoints = unique(timepoints)
  unique.dosages = unique(dosages)
  
  # Initialize a matrix to store the heatmap data with NA values, and rename the rows and columns
  heatmap.data = matrix(NA, nrow = length(unique.timepoints), ncol = length(unique.drugs) * length(unique.dosages))
  rownames(heatmap.data) = unique.timepoints
  colnames(heatmap.data) = paste(rep(unique.drugs, each = length(unique.dosages)), unique.dosages, sep = "_")
  
  # Fill the heatmap data
  for (j in seq_along(colnames(EXP.filtered.matrix.z))) {
    drug = drugs[j]
    dosage = dosages[j]
    timepoint = timepoints[j]
    heatmap.data[timepoint, paste(drug, dosage, sep = "_")] = EXP.filtered.matrix.z[1, j]
  }
  
  # Transform matrix to df for heatmap
  heatmap.data = as.data.frame(heatmap.data)
  
  # Safety check: remove columns with any NA values
  heatmap.data = heatmap.data[, colSums(is.na(heatmap.data)) == 0]
  
  # Plot
  # create_heatmap_z(heatmap.data, gene)
  print(create_bar_chart_z(heatmap.data, gene))
}

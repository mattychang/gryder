### Updated: 08.08.2024 (Matt)

### Goal: to build a TPM expression matrix from RNA-seq data
### Pre-requisites: raw expression files exist for each sample (.tdf), a list containing sample information (.txt)
### Notes: the sample list must contain 'Sample' as its first entry

### 1. Set directories, define samples to work on.
getwd()
setwd("/Volumes/rc/SOM_GENE_BEG33/RNA_seq/hg38/projects/RMS_IHK/Practice_MSC/SampleList")

## Read sample list from a file chosen from the client
## Note: this file should be a .txt file with the first entry labeled "Sample"
sample.list.all = read.table(file.choose(), header = T, sep = "\t")                                             # acquire the sample list
sample.list = sample.list.all$Sample
project.folder = "/Volumes/rc/SOM_GENE_BEG33/RNA_seq/hg38/projects/RMS_IHK/Practice_MSC"                        # define the project folder path
sample.set = "IHK_samples"                                                                                      # define the sample set name




### 2. Make coding TPM file, normalizing based on expected counts
setwd("/Volumes/rc/SOM_GENE_BEG33/RNA_seq/hg38/")
file.exists = file.exists(paste0("DATA/", sample.list,"/", sample.list, ".gene.TPM.txt", sep=""))               # check to see if TPM.txt exists for each sample
sample.list = sample.list[file.exists]                                                                          # remove samples that do not have TPM.txt

## Loop through each sample in the sample list
lapply(sample.list, function(x) {
  coding <- read.table("ref/RSEM_hg38/HGNC_protein-coding_gene_19229_2022.txt", sep="\t", header=T)             # load more stringent protein coding list (from Diana)
  EXP <- read.table(paste("DATA/",x,"/",x,".genes.results",sep=""), sep="\t", header=T)                         # load RSEM (genes.results) output file for the current sample with expected count, TPM, FPKM
  EXP$coding = EXP$gene_id %in% coding$symbol                                                                   # remove non-coding RNA entries
  EXP.coding <<- subset(EXP, EXP$coding %in% c("TRUE"))
  expected_sum = sum(EXP.coding$expected_count)                                                                 # sum all coding expected counts
  EXP.coding$count_norm = EXP.coding$expected_count / expected_sum * 1000000                                    # normalize expected counts to counts per million
  write.table(EXP.coding[,1:9],
              file=paste("DATA/",x,"/",x,".gene.coding.norm.txt",sep=""), 
              sep="\t", 
              row.names=FALSE, 
              col.names=TRUE, 
              quote=FALSE)                                                                                      # write the normalized data to a file
})                                                                                                              # note: outputs should all be NULL







### 3. Build TPM matrix
EXP.coding.matrix = EXP.coding["gene_id"]                                                                       # initiate matrix with gene IDs
lapply(sample.list, function(x) {                                                                               # loop through each sample, again, to extract TPM (normalized counts)
    EXP.sample = read.table(paste("DATA/",x,"/",x,".gene.coding.norm.txt",sep=""), sep="\t", header=T)         # read the normalized counts column
    EXP.sample = as.data.frame(EXP.sample[,9])                                                                  # extract the normalized counts column
    removable.string = "Sample_"
    sample.name = gsub(removable.string,"",x)
    colnames(EXP.sample) = c(sample.name)
    EXP.coding.matrix <<- cbind(EXP.coding.matrix, EXP.sample)                                                  # combine the current sample data into the matrix
})
colnames(EXP.coding.matrix) = gsub("-", "", colnames(EXP.coding.matrix))                                       # remove hyphen

## Write the final expression matrix to a file
dir.create(file.path(project.folder, "ExpMatrices"))
write.table(EXP.coding.matrix, 
            file=paste(project.folder,"/", "ExpMatrices/",sample.set,".coding.norm.matrix.txt",sep=""), 
            sep="\t", 
            row.names=FALSE, 
            col.names=TRUE, 
            quote=FALSE)


################### end of building the TPM matrix ###################

### Goal: build TPM expression matrix from RNA-seq data
### Pre-requisites: .tdf files have been made for each sample, a list containing sample information
### Notes: the sample list must contain 'Sample' as its first entry

### 1. Set directories, define samples to work on.
getwd()
setwd("/Volumes/rc/SOM_GENE_BEG33/RNA_seq/hg38/projects/IHK_RMS/Practice_MSC")

## Read sample list from a file chosen from the client
## Note: this file should be a .txt file with the first entry labeled "Sample"
sample.list.all = read.table(file.choose(), header = T, sep = "\t")                           # acquiring the sample list
sample.list = sample.list.all$Sample
project.folder = "/Volumes/rc/SOM_GENE_BEG33/RNA_seq/hg38/projects/RMS_IHK/Practice_MSC"      # define the project folder path
sample.set = "IHK_samples"                                                                    # define the sample set name




### 2. Make coding TPM file, normalizing based on expected counts
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






### 3. Build TPM matrix
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
dir.create(file.path(project.folder, "ExpMatrices"))
write.table(EXP.coding.matrix, 
            file=paste(project.folder,"/", "ExpMatrices/",sample.set,".coding.norm.matrix.txt",sep=""), 
            sep="\t", 
            row.names=FALSE, 
            col.names=TRUE, 
            quote=FALSE)


################### end of building the TPM matrix ###################

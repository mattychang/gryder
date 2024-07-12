#########################################################################################################
##  Collect GSEA results.  Works with a limited set of gene sets all in the same GSEA output folder.   ##
#########################################################################################################
##  Berkley Gryder, 2018-2019 (gryderart@gmail.com): https://github.com/GryderArt/VisualizeRNAseq      ##
#########################################################################################################

######################## start of gathering GSEA data ########################

### 1. Setup
setwd("/Volumes/rc/SOM_GENE_BEG33/RNA_seq/hg38/projects/RMS_IHK/Practice_MSC/")
GSEA.folder = "GSEA_results"

## Defining samples to work on
## Note: Do not use "." in GSEA Analysis names, will cause incorrect parsing
samples <- list.dirs(path = GSEA.folder, 
                     full.names = F, 
                     recursive = F)

## OPTIONAL: subselect samples
samples <- samples[grepl("IHK44_1uM_6h|IHK44_100nM_6h|A485_1uM_6h|A485_100nM_6h|dCBP_1uM_6h|dCBP_100nM_6h", 
                         samples)];

### 2. Edit string from folder names to get report names
## Create a df from the sample names
sample.df = read.table(text = samples, sep=".")                                                     # read samples into the df parsing through "."
sample.df$V3 = as.character(sample.df$V3)                                                           # convert third column to char
colnames(sample.df) = c("Drug","GSEA_Style","Digits")                                               # change column names

## Generate paths to negative and positive GSEA report files
sample.df$neg.xls.path = paste(GSEA.folder, 
                               "/", 
                               paste(sample.df$Drug, sample.df$GSEA_Style, sample.df$Digits, sep="."),
                               "/gsea_report_for_na_neg_",
                               sample.df$Digits,
                               ".tsv",
                               sep="")
sample.df$pos.xls.path = paste(GSEA.folder, 
                               "/", 
                               paste(sample.df$Drug, sample.df$GSEA_Style, sample.df$Digits, sep="."),
                               "/gsea_report_for_na_pos_",
                               sample.df$Digits,
                               ".tsv",
                               sep="")

### 3. Import GSEA data 
library(plyr)
column_names <- c("NAME", "ES", "NES", "NOM.p.val", "FDR.q.val", "FWER.p.val", "Drug")
allsamples <- data.frame(matrix(ncol = length(column_names), nrow = 0))
colnames(allsamples) <- column_names
df.empty <- allsamples


## Loop through each row in the 'Drug' column
lapply(sample.df$Drug, function(x) {
  
    # create temp df to subselect sample in loop
    temp.df <- subset(sample.df,sample.df$Drug %in% x)
    
    # check if neg.xls.path for the sample exists, and add to df.neg if true
    if(file.exists(temp.df$neg.xls.path) == 'TRUE') {
      df.neg = read.table(temp.df$neg.xls.path,sep="\t",header=T)
    }
    else {
      df.neg = df.empty
    }
    
    # check if pos.xls.path for the sample exists, and add to df.pos if true
    if(file.exists(temp.df$pos.xls.path) == 'TRUE') {
      df.pos = read.table(temp.df$pos.xls.path,sep="\t",header=T)
    }
    else {
      df.pos = df.empty
    }
    
    # combine df.neg and df.pos
    df = rbind(df.pos,df.neg)
    
    # select specific columns and add a 'Drug' column for identification
    df = df[,c("NAME","ES","NES","NOM.p.val","FDR.q.val","FWER.p.val")]
    df$Drug = x
    
    # append the modified dataframe 'df' to 'allsamples'
    allsamples <<- rbind(allsamples, df)
})

## OPTIONAL: simplify naming convention
allsamples$NAME = gsub(allsamples$NAME, pattern = "PAX3FOXO1", replacement = "P3F")
allsamples$NAME = gsub(allsamples$NAME, pattern = "ENHANCERS", replacement = "ENH")

## Reindex the rows of 'allsamples'
row.names(allsamples) <- 1:nrow(allsamples)

## Ensure 'NOM.p.val' is numeric and add a small value to avoid log(0)
as.numeric(allsamples$NOM.p.val) + 0.00001

# Calculate some statistics
allsamples$NOM.p.val[is.na(allsamples$NOM.p.val)] <- 0                                          # Replace NA values in the 'NOM.p.val' column with 0
allsamples$NOM.p.val = as.numeric(allsamples$NOM.p.val) + 0.00001                               # Convert the 'NOM.p.val' column to numeric type and add a small value (0.00001) to avoid log transformation issues with zero values
allsamples$log10.NOM.p.val= log10(allsamples$NOM.p.val)                                         # Calculate the base-10 logarithm of 'NOM.p.val' and store it in a new column 'log10.NOM.p.val'
allsamples$NES = as.numeric(allsamples$NES)                                                     # Ensure the 'NES' column is of numeric type
    
## Rorder by Enrichment Score (NES)
NAME.ranks = aggregate(NES ~ NAME, allsamples, mean)
colnames(NAME.ranks) = c("NAME", "rankmetric")
allsamples <- join(allsamples, NAME.ranks, by = "NAME")
allsamples$NAME <- factor(allsamples$NAME, levels = allsamples[order(unique(allsamples$rankmetric, decreasing=F)),]$NAME)               

allsamples <- na.omit(allsamples)

## Calculate minimum ES and NES
ES.min = min(abs(allsamples$ES))
NES.min = min(abs(allsamples$NES))
    
#allsamples$NAME <- factor(allsamples$NAME, levels = c("HALLMARK_APOPTOSIS","GRYDER_HOUSEKEEPING","TFS_NO_EPIMACHINES","GRYDER_RH4_SE_GENES","GRYDER_RH4_TOP_SE_TFS")) #order genesets

# OPTIONAL: save the 'allsamples' data frame to a text file
write.table(allsamples, 
            paste(GSEA.folder, ".GSEA.allsamples.summary.txt",sep=""),
            col.names = T,
            row.names = F,
            quote = F,
            sep="\t")








######################## start of making GSEA plots ########################

### Bubble chart and heatmap plots

## Setup
library(ggplot2)
library(pheatmap)
library(tidyr)
library(grid)
library(gridExtra)

## OPTIONAL: subset
GENE.plot.subset = subset(allsamples, 
                          allsamples$Drug %in% c("A485_1uM_6h",
                                                 "A485_100nM_6h",
                                                 "dCBP_1uM_6h",
                                                 "dCBP_100nM_6h",
                                                 "IHK44_1uM_6h",
                                                 "IHK44_100nM_6h")) 

## Bubble chart of ES values
ggplot(GENE.plot.subset, aes(x = Drug, y = NAME, size = abs(ES) - ES.min, colour = ES)) + 
  geom_point() +
  scale_colour_gradient2(low = "tomato",mid = "white", high = "lightskyblue") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

## Bubble chart of -log10(NOM.p.val) values
ggplot(GENE.plot.subset, aes(x = Drug, y = NAME, size = -log10.NOM.p.val, colour = ES)) + 
  geom_point() +
  scale_colour_gradient2(low = "blue",mid = "white", high = "orange") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
        
## Heatmap of enrichments
# Prepare data
heatmap.data = GENE.plot.subset[,c("NAME","NES","Drug")]
heatmap.data.wide = spread(heatmap.data, NAME, NES)

# Convert the data frame to a matrix
heatmap.matrix = as.matrix(heatmap.data.wide[,-1])
rownames(heatmap.matrix) = heatmap.data.wide$Drug

# Generate the heatmap with clustering
pheatmap(heatmap.matrix,
         scale = "row", 
         cluster_rows = T, 
         cluster_cols = T,
         main = paste(GSEA.folder," 
                    GSEA heatmap",
                    sep=""))




### GSEA ranklist

## Define the gene set of interest
Geneset = "GRYDER_RH4_TOP_CRTFS"                                                #HALLMARK_MYC_TARGETS_V2
Geneset.allsamples <- data.frame(RANK.IN.GENE.LIST = numeric(), 
                                 RUNNING.ES = numeric(), 
                                 SYMBOL = character(), 
                                 Drug = character())


## Loop through each drug in sample.df    
lapply(sample.df$Drug, function(x) {
    
    # Subset the sample.df data frame for the current drug
    ESplot.sample <- subset(sample.df, sample.df$Drug %in% x)
    
    # Create the path to the gene set file for the current drug
    ESplot.sample$geneset.path = paste(GSEA.folder, "/", paste(ESplot.sample$Drug,
                                                               ESplot.sample$GSEA_Style,
                                                               ESplot.sample$Digits,
                                                               sep="."),
                                       "/", Geneset, ".tsv", sep="")
    
    # Read the gene set data
    Geneset.df = read.table(ESplot.sample$geneset.path, sep = "\t", header = T)
   
    # Select relevant columns and add a column for drug identification
    Geneset.df = Geneset.df[,c("RANK.IN.GENE.LIST", "RUNNING.ES", "SYMBOL")]
    Geneset.df$Drug = x
    
    # Combine the data for all samples
    Geneset.allsamples <<- rbind(Geneset.allsamples, Geneset.df)
})

## OPTIONAL: write out Geneset.allsamples table
RH4_CRTFs_path <- file.path(project.folder, paste0("RH4_CRTFs.GSEA_matrix.txt"))
write.table(Geneset.allsamples, 
            file = RH4_CRTFs_path, 
            sep = "\t", 
            row.names = F, 
            col.names = T, 
            quote = FALSE, 
            append = FALSE)
        
# Filtering for specific gene to add to the ranklist
Geneset.allsamples = na.omit(Geneset.allsamples)                                                # remove NA values
Geneset.allsamples = subset(Geneset.allsamples, 
                            Geneset.allsamples$Drug %in% c("A485_1uM_6h",
                                                           "A485_100nM_6h",
                                                           "dCBP_1uM_6h",
                                                           "dCBP_100nM_6h",
                                                           "IHK44_1uM_6h",
                                                           "IHK44_100nM_6h"))                   # subselect drug/samples
# OPTIONAL: pick a gene of interest
Geneset.allsamples.pick = subset(Geneset.allsamples, 
                                 Geneset.allsamples$SYMBOL %in% c("SOX8"))                      # subselect gene

## Define color schematic for different drugs. There should be n numbercolors for n drugs
## http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
Drugcolors = c("firebrick", "firebrick1", "antiquewhite4","blue", "chartreuse4", "violet")

# Plot 1: Enrichment Score plot
p1 = ggplot(Geneset.allsamples, aes(x = RANK.IN.GENE.LIST, y = RUNNING.ES, group = Drug, color = Drug)) + 
  geom_line(size = 1.2) + 
  geom_hline(yintercept = 0, lty = 'dashed') +
  scale_color_manual(values=Drugcolors) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = c(0.25, 0.25)) +
  ylab("Enrichment Score") + 
  xlab(paste("Genes Ranked by log2 fold change")) +
  ggtitle(paste("Geneset:", Geneset,sep=" ")) +                                                 # optional
  geom_point(data = Geneset.allsamples.pick, color = "black") +                                 # optional
  geom_text(data = Geneset.allsamples.pick, label = Geneset.allsamples.pick$SYMBOL)             # optional

# Plot 2: Density plot
p2 = ggplot(Geneset.allsamples,group = Drug) + 
  stat_density(aes(x = RANK.IN.GENE.LIST, y = 0.5, fill = ..density..), geom = "tile", position="identity") + 
  facet_wrap(~Drug, ncol=1) + 
  scale_fill_gradient(low = "white", high = "darkgrey") +
  theme(axis.line = element_blank(), 
        axis.title.x = element_blank(), 
        legend.position = "none",
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background=element_blank()
        )

# Plot 3: Enhanced density plot with linerange
p3 = p2 + geom_linerange(aes(x = RANK.IN.GENE.LIST, ymin = 0, ymax = 1, color = Drug)) +
  facet_wrap(~Drug, ncol=1) +
  theme_bw() + 
  scale_color_manual(values=Drugcolors) +
  theme(axis.line = element_blank(),
        axis.title.x = element_blank(), 
        legend.position = "none",
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background=element_blank()
        )

# Arrange and display the plots
grid.arrange(arrangeGrob(p1,ncol=1))
grid.arrange(arrangeGrob(p1,p3,ncol=1))
        
         



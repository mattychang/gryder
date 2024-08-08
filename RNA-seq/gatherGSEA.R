### Updated: 08.08.2024 (Matt)

### 1. Setup: load libraries, set project.folder, define samples
library(plyr)
library(ggplot2)
library(pheatmap)
library(tidyr)
library(grid)
library(gridExtra)

setwd("/Volumes/rc/SOM_GENE_BEG33/RNA_seq/hg38/projects/RMS_IHK/Practice_MSC/GSEA_results")
GSEA.folder = "GSEA_results_against_DMSO"

project.folder = "/Volumes/rc/SOM_GENE_BEG33/RNA_seq/hg38/projects/RMS_IHK/Practice_MSC/GSEA_results"

# Defining samples to work on and define paths
samples = list.dirs(path = GSEA.folder, full.names = F, recursive = F)                                          # note: do NOT use "." in GSEA Analysis names, will cause incorrect parsing
sample.df = read.table(text = samples, sep=".")                                                                 # read samples into the df parsing through "."
sample.df$V3 = as.character(sample.df$V3)                                                                       # convert third column to char
colnames(sample.df) = c("Drug","GSEA_Style","Digits")                                                           # change column names
sample.df$neg.xls.path = paste(GSEA.folder, "/", paste(sample.df$Drug, sample.df$GSEA_Style, sample.df$Digits, sep="."), "/gsea_report_for_na_neg_", sample.df$Digits,".tsv", sep="")
sample.df$pos.xls.path = paste(GSEA.folder, "/", paste(sample.df$Drug, sample.df$GSEA_Style, sample.df$Digits, sep="."), "/gsea_report_for_na_pos_", sample.df$Digits, ".tsv", sep="")

### 2. Import GSEA data 
column.names = c("NAME", "ES", "NES", "NOM.p.val", "FDR.q.val", "FWER.p.val", "Drug")
allsamples = data.frame(matrix(ncol = length(column.names), nrow = 0))
colnames(allsamples) = column.names
df.empty = allsamples


# Loop through each row in the 'Drug' column
lapply(sample.df$Drug, function(x) {
  temp.df <- subset(sample.df,sample.df$Drug %in% x)                                                            # create temp df to subselect sample in loop
  if(file.exists(temp.df$neg.xls.path) == 'TRUE') {                                                             # check if neg.xls.path for the sample exists, and add to df.neg if true
    df.neg = read.table(temp.df$neg.xls.path,sep="\t",header=T)
  }
  else {
    df.neg = df.empty
  }
  if(file.exists(temp.df$pos.xls.path) == 'TRUE') {                                                             # check if pos.xls.path for the sample exists, and add to df.pos if true
    df.pos = read.table(temp.df$pos.xls.path,sep="\t",header=T)
  }
  else {
    df.pos = df.empty
  }
  df = rbind(df.pos,df.neg)                                                                                     # combine df.neg and df.pos
  df = df[,c("NAME","ES","NES","NOM.p.val","FDR.q.val","FWER.p.val")]                                           # select specific columns and add a 'Drug' column for identification
  df$Drug = x
  allsamples <<- rbind(allsamples, df)                                                                          # append the modified dataframe 'df' to 'allsamples'
})

# OPTIONAL: simplify naming convention
allsamples$NAME = gsub(allsamples$NAME, pattern = "PAX3FOXO1", replacement = "P3F")
allsamples$NAME = gsub(allsamples$NAME, pattern = "ENHANCERS", replacement = "ENH")

row.names(allsamples) <- 1:nrow(allsamples)                                                                     # reindex the rows of 'allsamples'
as.numeric(allsamples$NOM.p.val) + 0.00001                                                                      # ensure the 'NOM.p.val' column is numeric and add a small value to avoid log(0)
allsamples$NOM.p.val[is.na(allsamples$NOM.p.val)] <- 0                                                          # replace NA values in the 'NOM.p.val' column with 0
allsamples$NOM.p.val = as.numeric(allsamples$NOM.p.val) + 0.00001                                               # convert the 'NOM.p.val' column to numeric type and add a small value (0.00001) to avoid log transformation issues with zero values
allsamples$log10.NOM.p.val= log10(allsamples$NOM.p.val)                                                         # calculate the base-10 logarithm of 'NOM.p.val' and store it in a new column 'log10.NOM.p.val'
allsamples$NES = as.numeric(allsamples$NES)                                                                     # ensure the 'NES' column is of numeric type
    
# Rorder by Enrichment Score (NES)
NAME.ranks = aggregate(NES ~ NAME, allsamples, mean)
colnames(NAME.ranks) = c("NAME", "rankmetric")
allsamples = join(allsamples, NAME.ranks, by = "NAME")
allsamples$NAME = factor(allsamples$NAME, levels = allsamples[order(unique(allsamples$rankmetric, decreasing=F)),]$NAME)               

allsamples <- na.omit(allsamples)

# Calculate minimum ES and NES
ES.min = min(abs(allsamples$ES))
NES.min = min(abs(allsamples$NES))


# OPTIONAL: save the 'allsamples' data frame to a .txt file
write.table(allsamples, 
            paste(GSEA.folder, ".GSEA.allsamples.summary.txt",sep=""),
            col.names = T,
            row.names = F,
            quote = F,
            sep="\t")







######################## end of inputting GSEA info ########################
######################## start of making GSEA plots ########################

### 3. Bubble chart and heatmap plots of GSEA genesets

# OPTIONAL: subset samples
GENE.plot.subset = subset(allsamples, allsamples$Drug %in% c("IHK44_1uM_2h",
                                                             "IHK44_1uM_6h",
                                                             "IHK44_100nM_2h",
                                                             "IHK44_100nM_6h")) 

# Bubble chart of ES values
ggplot(GENE.plot.subset, aes(x = Drug, y = NAME, size = abs(ES) - ES.min, colour = ES)) + 
  geom_point() +
  scale_colour_gradient2(low = "tomato",mid = "white", high = "lightskyblue") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Bubble chart of -log10(NOM.p.val) values
ggplot(GENE.plot.subset, aes(x = Drug, y = NAME, size = -log10.NOM.p.val, colour = ES)) + 
  geom_point() +
  scale_colour_gradient2(low = "blue",mid = "white", high = "orange") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
        
# Heatmap of enrichments
heatmap.data = GENE.plot.subset[,c("NAME","NES","Drug")]                                                        # select relevant columns
heatmap.data.wide = spread(heatmap.data, NAME, NES)                                                             # reshape df to wide format
heatmap.matrix = as.matrix(heatmap.data.wide[,-1])                                                              # convert df to matrix, exclude first column
rownames(heatmap.matrix) = heatmap.data.wide$Drug                                                               # set row names to sample name

# Create heatmap
pheatmap(heatmap.matrix,
         scale = "row", 
         cluster_rows = T, 
         cluster_cols = T,
         main = paste(GSEA.folder," GSEA heatmap", sep=""))


### 4. BG overlapping GSEA ranklist

## Define the geneset of interest
Geneset = "GRYDER_RH4_CRTFS_NOPANESSENTIALS"                                                                    #HALLMARK_MYC_TARGETS_V2
Geneset.allsamples = data.frame(RANK.IN.GENE.LIST = numeric(), 
                                RUNNING.ES = numeric(), 
                                SYMBOL = character(), 
                                Drug = character())                                                             # initialize empty df to store results


## Loop through each drug in sample.df    
lapply(sample.df$Drug, function(x) {
  ESplot.sample <- subset(sample.df, sample.df$Drug %in% x)                                                     # subset the data for the current drug
  ESplot.sample$geneset.path = paste(GSEA.folder, 
                                     "/", 
                                     paste(ESplot.sample$Drug, ESplot.sample$GSEA_Style, ESplot.sample$Digits, sep="."),
                                     "/", 
                                     Geneset, 
                                     ".tsv", 
                                     sep="")                                                                    # create file path
  Geneset.df = read.table(ESplot.sample$geneset.path, sep = "\t", header = T)                                   # read gene set data
  Geneset.df = Geneset.df[,c("RANK.IN.GENE.LIST", "RUNNING.ES", "SYMBOL")]                                      # select specific columns from the data
  Geneset.df$Drug = x                                                                                           # add the current drug name as a new column
  Geneset.allsamples <<- rbind(Geneset.allsamples, Geneset.df)                                                  # append the current gene set data to the main data frame
})

# OPTIONAL: save the 'allsamples' data frame to a .txt file
Geneset.allsamples = na.omit(Geneset.allsamples)                                                                # remove NA values
RH4_CRTFs_path = file.path(project.folder, paste0("RH4_CRTFs.GSEA_matrix.txt"))                                 # define file path
write.table(Geneset.allsamples, 
            file = RH4_CRTFs_path, 
            sep = "\t", 
            row.names = F, 
            col.names = T, 
            quote = F, 
            append = F)
        
# OPTIONAL: filter for specific drugs
Geneset.allsamples = subset(Geneset.allsamples, Geneset.allsamples$Drug %in% c("A485_1uM_6h", 
                                                                               "dCBP_1uM_6h", 
                                                                               "IHK44_1uM_6h"))                 # subselect drug/samples

# OPTIONAL: define color scheme                                                                                 # note: there should be n number of colors for n drugs # http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
Drugcolors = c("coral1", 
               "chartreuse4", 
               "orange")

# OPTIONAL: pick a gene of interest
Geneset.allsamples.pick = subset(Geneset.allsamples, Geneset.allsamples$SYMBOL %in% c("MYCN"))                  # subselect gene

# Plot 1: Enrichment Score plot
p1 = ggplot(Geneset.allsamples, aes(x = RANK.IN.GENE.LIST, y = RUNNING.ES, group = Drug, color = Drug)) + 
  geom_line(size = 1.2) + 
  geom_hline(yintercept = 0, lty = 'dashed') +
  scale_color_manual(values=Drugcolors) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = c(0.25, 0.25)) +
  ylab("Enrichment Score") + 
  xlab(paste("Genes Ranked by log2 fold change")) # +
  # ggtitle(paste("Geneset:", Geneset,sep=" ")) +                                                 # optional
  # geom_point(data = Geneset.allsamples.pick, color = "black") +                                 # optional
  # geom_text(data = Geneset.allsamples.pick, label = Geneset.allsamples.pick$SYMBOL)             # optional

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

# Plot 3: Combine previous two plots
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

#################### end of making GSEA plots ####################

####################################################
# Script: Correlation of Gene Expression and Methylation
#
# Description:
#      This script performs the correlation analysis between gene expression and methylation data. 
#      It processes an expression matrix and methylation matrix, performs correlation analysis 
#      on a gene-probe level, and generates correlation plots for selected genes.
#
# Requirements:
#      - Custom correlation functions from `main.functions.R` and `corr.met.exp.functions.R`.
#      - Configuration file (`config.tsv`) containing directory paths for working and data directories.
#      - Clinical data file (`clinical_data.txt`), expression data file (`Expression.NS.txt`), 
#        and methylation data file (`Methylation.NS.txt`).
#      - Methylation annotation file in CSV format.
####################################################
# Required environment: correlation.exp.met.yml.

# Load required libraries for data manipulation, statistical modeling, and plotting
library(tidyverse)
library(limma)
library(ComplexHeatmap)
library(circlize)
library(ggpubr)

# Load custom correlation functions for gene expression and methylation analysis
source("functions/main.functions.R")
source("functions/clustering/corr.met.exp.functions.R")


############### Configuration and Directory Setup ###############
# Read the configuration file (config.tsv) to retrieve the directory paths

config <- read.csv("config/config.tsv",
sep = "\t",
header = T)

# Assign configuration keys and values
keys <- config$key
values <- config$value
names(values) <- keys

# Define working directory and data directory based on configuration values
wkdir <- values["wkdir"]
datadir <- values["datadir"]


# Check if the output directory for correlation analysis exists, create if it doesn't
if (dir.exists(paste0(wkdir, "/clustering/corr_exp_met")) == F) {

    dir.create(paste0(wkdir, "/clustering/corr_exp_met"))

}


# Set output directory for correlation results
outdir <- paste0(wkdir, "/clustering/corr_exp_met")


############### Data Loading ###############
# Load clinical data from the specified file
base <- read.csv(paste0(datadir, "/clinical_data.txt"),
                 sep = "\t",
                 row.names = 1,
                 header = T)

# Load expression data from the specified file
exp <- read.table(paste0(wkdir,
"/clustering/deg/Expression.NS.txt"),
    sep = "\t",
    header = T,
    row.names = 1,
    check.names = F
)

# Load methylation data from the specified file
met <- read.table(paste0(wkdir,
"/clustering/dep/Methylation.NS.txt"),
    sep = "\t",
    header = T,
    row.names = 1,
    check.names = F
)

# Load methylation annotation file from the specified CSV file
epic <- read.csv(paste0(datadir,
"/methylome/anotation/infinium-methylationepic-v-1-0-b5-manifest-file.csv"),
                       sep = ",",
                       header = T,
                       skip = 7)

# Create methylation annotation using a custom function
met_anot <- create.met.anot(epic)

############### Correlation Analysis ###############
# Perform correlation analysis between gene expression and methylation
df.corr <- corr.met.exp.unique(exp, met, met_anot)

# Adjust p-values using the FDR method and merge with annotation data
df.corr$p.adj <- p.adjust(df.corr$p.val, method = "fdr")
df.corr <- merge(df.corr, met_anot, by.x = "CpG", by.y = "ID")
df.corr <- df.corr[order(df.corr$p.val),]

# Save the correlation results to a file
write.table(df.corr,
paste0(outdir,"/Correlation_exp_met.txt"),
sep = "\t",
row.names = T,
col.names = NA)




############### Correlation per Primary Tumor ###############
# Perform correlation analysis per primary tumor type

tissues <- unique(base$PRIMARY_TUMOR)

corr.tissue <- lapply(tissues, function(x) {
    
    # Filter samples by primary tumor type    
    df <- base[base$PRIMARY_TUMOR == x,]

    # Subset expression and methylation data by selected samples
    exp.df = exp[,rownames(df)]
    met.df = met[,rownames(df)]
 # Perform correlation analysis for the subset
    corr <- corr.met.exp.unique(exp.df, met.df, met_anot)
    corr$p.adj <- p.adjust(corr$p.val, method = "fdr")  # Adjust p-values
    corr$tissue <- rep(x, nrow(corr))  # Add tissue type information
    corr <- merge(corr, met_anot, by.x = "CpG", by.y = "ID")  # Merge with annotation
    corr <- corr[order(corr$p.val),]  # Sort by p-value
    return(corr)

})

# Combine correlation results for all tissues and save to a file
corr.tissue.df <- as.data.frame(do.call(rbind, corr.tissue))

write.table(corr.tissue.df,
paste0(outdir,"/Correlation_exp_met_primary.txt"),
sep = "\t", row.names = F, col.names = T)

############### Highest Correlated CpG Per Gene ###############
# Select the highest correlated CpG site per gene

df.corr <- df.corr[!duplicated(df.corr$gene),]


# Filter for significant correlations (adjusted p-value < 0.05)
df.corr <- df.corr[df.corr$p.adj < 0.05,]
df.corr.rep <- df.corr


# Update gene symbols in the expression matrix
symbol_genes <- alias2SymbolTable(rownames(exp))
rownames(exp)[is.na(symbol_genes) == F] <- symbol_genes[is.na(symbol_genes) == F]

# Subset the methylation matrix by the correlated CpG sites
met.corr <- met[df.corr$CpG,]
rownames(met.corr) <- df.corr$gene


############### Correlation Plot Generation ###############
# List of genes to generate correlation plots for

genes <- c("PTPRN", "PTPRN2", "CHGA", "PIGR", "CHGB")



# Generate and save correlation plots for the selected genes
correlation.plot(exp = exp,
met = met,
anot = df.corr.rep,
base = base,
genes = genes,
outdir = outdir
)


#
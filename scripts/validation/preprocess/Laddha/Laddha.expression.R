####################################################
# Script: RNA-seq Data Processing and Normalization for Laddha LUNET Dataset
#
# Description:
#     This script processes RNA-seq data from the Laddha LUNET dataset, normalizes 
#     it using CPM (Counts Per Million), annotates genes, and outputs a final 
#     expression matrix for downstream analysis.
#
# Requirements:
#     - `geo` environment for necessary packages.
#     - Configuration file (`config/config.tsv`) for directory paths.
#     - RNA-seq expression data file (`counts.tsv`).
#     - Custom normalization and annotation functions from `validation.functions.R`.
#
# Inputs:
#     - RNA-seq expression data file (`counts.tsv`).
#
# Outputs:
#     - Normalized and annotated expression matrix (`matrix_Laddha_expression.txt`).
#
# Required environment: `geo` (containing relevant packages for data processing).
####################################################

############### Load Required Libraries ###############
# Load libraries for data manipulation, normalization, and annotation

library(tidyverse)
library(edgeR)
library(limma)

# Load custom validation functions
source("functions/validation/validation.functions.R")

############### Configuration and Directory Setup ###############
# Load configuration file to get working and data directory paths
config <- read.csv("config/config.tsv",
sep = "\t",
header = T)
# Extract key-value pairs for configuration settings
keys <- config$key
values <- config$value
names(values) <- keys
# Define working directory and data directory from the configuration file
wkdir <- values["wkdir"]
datadir <- values["datadir"]

# Create output directories if they don't exist
if (dir.exists(paste0(wkdir, "/validation")) == F) {

    dir.create(paste0(wkdir, "/validation"))

}

if (dir.exists(paste0(wkdir, "/validation/Laddha")) == F) {

    dir.create(paste0(wkdir, "/validation/Laddha"))

}
# Define the output directory for results
outdir <- paste0(wkdir, "/validation/Laddha")

############### Data Loading ###############
# Load RNA-seq expression data from the counts file
expression <- read.table(paste0(datadir,"/external_data/Laddha_LUNET/counts.tsv"),
sep = "\t",
header = T,
row.names = 1)

############### Normalization ###############
# Normalize the RNA-seq data using CPM (Counts Per Million)
expression <- normalization(expression,"cpm")


############### Gene Annotation ###############
# Get the list of genes (Ensembl IDs) from the expression matrix
genes <- rownames(expression)

# Convert alias gene symbols to official gene symbols using limma
symbol_genes <- limma::alias2SymbolTable(genes)

# Create a new data frame with annotated gene symbols
df_expression <- data.frame("Expression" = genes ,
"Symbol" = symbol_genes)

# Merge the gene symbols back into the expression matrix
expression <- merge(df_expression,expression,by.x = "Expression",
by.y = "row.names", all.y = T)

# Replace missing symbols with the original expression values (if no official symbol exists)
expression$Symbol[is.na(expression$Symbol)] <- expression$Expression[is.na(expression$Symbol)] 


############### Variance Filtering and Data Cleanup ###############
# Calculate variance for each gene across all samples

variance <- apply(expression[,4:ncol(expression)],1,var)

# Order genes by variance (highest to lowest) and remove duplicates
expression <- expression[order(variance, decreasing = T),]
expression <- expression[!duplicated(expression$Symbol),]

# Remove rows with missing gene symbols
expression <- expression[complete.cases(expression$Symbol),]

# Set gene symbols as row names and retain only the expression data columns
rownames(expression) <- expression$Symbol
expression <- expression[,3:ncol(expression)]
#######################################################################

############### Save the Final Processed Expression Matrix ###############
# Write the processed and normalized expression matrix to a file
write.table(expression,
paste0(outdir,"/matrix_Laddha_expression.txt"),
sep = "\t",
row.names = T,
col.names = NA)

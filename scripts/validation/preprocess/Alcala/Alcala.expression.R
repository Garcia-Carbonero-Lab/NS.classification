####################################################
# Script: Normalization and Annotation of Alcala and Laddha Gene Expression Data
#
# Description:
#     This script processes gene expression data from the Alcala amd Laddha dataset. It performs 
#     normalization, gene annotation, and filtering to retain high-variance genes. 
#     The processed expression matrix is then saved for downstream analysis.
#
# Requirements:
#     - GEOquery: For handling GEO datasets (if applicable).
#     - edgeR: For normalization of gene expression data.
#     - limma: For gene symbol mapping and annotation.
#     - Custom functions from `validation.functions.R`.
#     - Configuration file (`config/config.tsv`) for specifying working and data directories.
#     - Gene expression data file (`read_counts_all.txt`) and metadata file (`Attributes.txt`).
####################################################
# Required environment: 'geo.yml' (containing prerequiste libraries).

# Load required libraries for data processing and gene expression analysis
library(tidyverse)
library(edgeR)
library(limma)

# Load the configuration file (config.tsv) to retrieve directory paths
config <- read.csv("config/config.tsv",
sep = "\t",
header = T)

# Define working directory and data directory from the configuration file
keys <- config$key
values <- config$value
names(values) <- keys
###############
wkdir <- values["wkdir"]
datadir <- values["datadir"]

# Load custom validation functions if needed
source("functions/validation/validation.functions.R")

############### Directory Setup ###############
# Create directories for storing validation results if they don't exist
if (dir.exists(paste0(wkdir, "/validation")) == F) {

    dir.create(paste0(wkdir, "/validation"))

}

if (dir.exists(paste0(wkdir, "/validation/Alcala")) == F) {

    dir.create(paste0(wkdir, "/validation/Alcala"))

}

# Set output directory for storing processed Alcala data
outdir <- paste0(wkdir, "/validation/Alcala")



############### Expression Data Loading ###############
# Load raw gene expression counts data (read counts) from the specified file
expression <- read.table(paste0(datadir,"/external_data/Alcala_LUNET/read_counts_all.txt"),
sep = " ",
header = T)

# Load metadata (Attributes) that describes the samples
base <- read.table(paste0(datadir,"/external_data/Alcala_LUNET/Attributes.txt"),
sep = "\t",
header = T)


############### Data Preprocessing ###############
# Set gene identifiers (gene_id) as row names in the expression matrix

rownames(expression) <- expression$gene_id

# Remove the first column (gene_id) and retain only the expression values
expression <- expression[,2:ncol(expression)]

# Filter samples for Neuroendocrine Tumors (NETs) based on the metadata
Nets <- base$Sample_ID[which(base$Histopathology %in% c("Atypical", "Typical", "Carcinoid"))]

# Subset the expression matrix to include only the selected NET samples
Nets <- Nets[Nets %in% colnames(expression)]

expression <- expression[,Nets]


############### Normalization ###############
# Normalize the expression data using Counts Per Million (CPM) normalization
expression <- normalization(expression,"cpm")


############### Gene Annotation ###############
# Remove version numbers from Ensembl gene IDs to facilitate annotation
rownames(expression) <- gsub("\\..*","",rownames(expression))

# Perform gene annotation using Ensembl gene IDs
genes <- rownames(expression)

genesID <-mygene::queryMany(genes, scopes="ensembl.gene", fields="symbol", species="human")
genesID <- as.data.frame(genesID[,c("query","symbol")])


# Merge the gene annotation (symbol) with the expression data
expression <- merge(genesID, expression,
by.x = "query", by.y = "row.names")

# Update gene symbols using the limma package's alias2SymbolTable function
genes <- expression$symbol

symbol_genes <- limma::alias2SymbolTable(genes)

# Create a dataframe mapping gene symbols to expression data
df_expression <- data.frame("Expression" = genes <- expression$symbol,
"Symbol" = symbol_genes)

# Merge the symbol gene data into the expression matrix
expression <- merge(df_expression,expression,by.x = "Expression",
by.y = "symbol", all.y = T)

# Fill missing symbols with their corresponding original expression identifiers
expression$Symbol[is.na(expression$Symbol)] <- expression$Expression[is.na(expression$Symbol)] 



############### Gene Filtering by Variance ###############
# Calculate the variance of each gene across samples and sort by variance

# Calculate the variance of each gene across samples and sort by variance
variance <- apply(expression[,4:ncol(expression)],1,var)
expression <- expression[order(variance, decreasing = T),]

# Remove duplicate gene symbols, keeping the one with the highest variance
expression <- expression[!duplicated(expression$Symbol),]

# Remove genes with missing or incomplete data
expression <- expression[complete.cases(expression$Symbol),]

# Set gene symbols as row names and drop unused columns
rownames(expression) <- expression$Symbol
expression <- expression[,4:ncol(expression)]

############### Save Processed Expression Data ###############
# Write the filtered and normalized expression matrix to a file
write.table(expression,
paste0(outdir, "/Alcala.expression.tsv"),
sep = "\t",
row.names = T,
col.names = NA)
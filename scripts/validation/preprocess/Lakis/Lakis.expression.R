####################################################
# Script: Processing and Normalizing Lakis PANET RNA-seq Data
#
# Description:
#     This script processes RNA-seq data from the Lakis PANET dataset. It reads expression data, 
#     normalizes it using CPM (Counts Per Million), annotates genes, 
#     filters patient data, and prepares the final expression matrix for downstream analysis.
#
# Requirements:
#     - Configuration file (`config/config.tsv`) containing paths to working and data directories.
#     - RNA-seq expression data file (`rnaseq_Lakis.tsv`).
#     - Patient metadata file (`base_Lakis.tsv`).
#     - Custom normalization and annotation functions from `validation.functions.R`.
####################################################

# Load required libraries for data manipulation and RNA-seq normalization
library(tidyverse)
library(edgeR)
library(limma)
# library(mygene)  # This package is commented out but required for gene annotation

# Load custom validation functions
source("functions/validation/validation.functions.R")

############### Configuration and Directory Setup ###############
# Read the configuration file (config.tsv) to retrieve the directory paths
config <- read.csv("config/config.tsv", sep = "\t", header = TRUE)

# Assign configuration keys and values
keys <- config$key
values <- config$value
names(values) <- keys

# Define working directory and data directory based on configuration values
wkdir <- values["wkdir"]
datadir <- values["datadir"]

# Create output directories if they don't exist
if (!dir.exists(paste0(wkdir, "/validation"))) {
    dir.create(paste0(wkdir, "/validation"))
}

if (!dir.exists(paste0(wkdir, "/validation/Lakis"))) {
    dir.create(paste0(wkdir, "/validation/Lakis"))
}

# Set output directory for results
outdir <- paste0(wkdir, "/validation/Lakis")

############### Data Loading and Reshaping ###############
# Load RNA-seq expression data from the specified file
expression <- read.table(paste0(datadir, "/external_data/Lakis_PANET/rnaseq_Lakis.tsv"),
                         sep = "\t", header = TRUE)

# Reshape the expression data: organize genes as rows and samples as columns
expression <- expression %>%
    select(icgc_sample_id, gene_id, raw_read_count) %>%  # Select relevant columns
    pivot_wider(names_from = icgc_sample_id, values_from = raw_read_count)  # Pivot to wide format

# Convert the expression matrix to a data frame and set rownames as gene IDs
expression <- as.data.frame(expression)
rownames(expression) <- expression$gene_id
expression <- expression[, -1]  # Remove gene_id column, retain only the expression data

############### Normalization ###############
# Normalize the RNA-seq data using CPM method (Counts Per Million)
expression <- normalization(expression, "cpm")

############### Gene Annotation ###############
# Get the list of genes (Ensembl IDs) from the expression matrix
genes <- rownames(expression)

# Annotate genes using the `mygene` package (annotating Ensembl gene IDs with symbols)
genesID <- mygene::queryMany(genes, scopes = "ensembl.gene", fields = "symbol", species = "human")
genesID <- as.data.frame(genesID[, c("query", "symbol")])

# Merge the gene annotations with the expression matrix
expression <- merge(genesID, expression, by.x = "query", by.y = "row.names")

# Use official gene symbols for downstream analysis
genes <- expression$symbol

# Convert alias gene symbols to official symbols using the limma package
symbol_genes <- limma::alias2SymbolTable(genes)

# Create a new data frame with annotated gene symbols
df_expression <- data.frame("Expression" = genes, "Symbol" = symbol_genes)

# Merge the gene symbols back into the expression matrix
expression <- merge(df_expression, expression, by.x = "Expression", by.y = "symbol", all.y = TRUE)

# Replace missing symbols with the original expression values (if no official symbol exists)
expression$Symbol[is.na(expression$Symbol)] <- expression$Expression[is.na(expression$Symbol)]

# Remove duplicate rows based on gene query IDs
expression <- expression[!duplicated(expression$query),]

############### Variance Filtering and Data Cleanup ###############
# Calculate variance for each gene across all samples
variance <- apply(expression[, 4:ncol(expression)], 1, var)

# Order genes by variance (highest to lowest) and remove duplicates
expression <- expression[order(variance, decreasing = TRUE),]
expression <- expression[!duplicated(expression$Symbol),]

# Remove rows with missing gene symbols
expression <- expression[complete.cases(expression$Symbol),]

# Set gene symbols as row names and retain only the expression data columns
rownames(expression) <- expression$Symbol
expression <- expression[, 4:ncol(expression)]

############### Patient Annotation and Data Subsetting ###############
# Load patient metadata file
base <- read.table(paste0(wkdir, "/validation/Lakis/base_Lakis_prev.tsv"), sep = "\t", header = TRUE, row.names = 1)

# Subset the base data to match the expression matrix columns (samples)
base <- base[colnames(expression),]

# Filter out specific differentiation types (if needed)
base <- base[!base$Differenciation %in% c("Poorly differentiated endocrine carcinoma", 
                                           "Mixed ductal endocrine carcinoma"),]

# Subset the expression matrix based on the filtered base data
expression <- expression[, rownames(base)]

# Rename expression matrix columns to reflect patient IDs from the metadata
colnames(expression) <- base$id_patient

############### Save the Final Processed Expression Matrix ###############
# Write the processed and normalized expression matrix to a file
write.table(expression,
            paste0(outdir, "/matrix_Lakis_expression.txt"),
            sep = "\t", row.names = TRUE, col.names = NA)

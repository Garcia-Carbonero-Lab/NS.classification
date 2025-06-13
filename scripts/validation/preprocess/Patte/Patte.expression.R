############################################################
# Script: Gene Expression Data Processing for Patte Dataset
#
# Description:
#     This script processes gene expression data for the Patte dataset. 
#     It normalizes the data. The final processed 
#     data is saved in tab-separated value (TSV) files for further analysis.
############################################################
# Required environment: 'geo.yml' (containing prerequiste libraries).

# Load required libraries for expression data processing and analysis
library(tidyverse)
library(edgeR)
library(limma)

source("functions/validation/validation.functions.R")

############### Load Configuration ###############
# Read configuration file to define paths
config <- read.csv("config/config.tsv",
sep = "\t",
header = T)

# Extract configuration keys and values
keys <- config$key
values <- config$value
names(values) <- keys

# Set working directory and data directory from config
wkdir <- values["wkdir"]
datadir <- values["datadir"]


############### Directory Setup ###############
# Create directories for validation results if they don't exist
if (dir.exists(paste0(wkdir, "/validation")) == F) {

    dir.create(paste0(wkdir, "/validation"))

}

if (dir.exists(paste0(wkdir, "/validation/Patte")) == F) {

    dir.create(paste0(wkdir, "/validation/Patte"))

}

# Set output directory for processed Karpathakis data
outdir <- paste0(wkdir, "/validation/Patte")

# Load expression data
expression <- readRDS(paste0(datadir,"/external_data/SiNETs_Patte/expression_data_genesymbol.rds"))
expression <- expression$counts
#Normalize data
expression <- normalization(expression,"cpm")

write.table(expression,
paste0(outdir, "/Patte.expression.tsv"),
sep = "\t",
row.names = T,
col.names = NA)
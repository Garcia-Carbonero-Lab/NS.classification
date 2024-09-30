# This script processes transcriptome and methylome data to filter features
# based on Mean Absolute Deviation (MAD) within specific tissue types.
# The script relies on preprocessed data and predefined functions. 
# The output is filtered expression and methylation matrices for downstream analysis.
#
# Requirements:
# - Preprocessed expression and methylation data with quality adjustments.
# - The script filters genes based on their MAD across multiple tissue types.
# - Configurations for directories are specified in "config/config.tsv".
# - Uses external filtering functions from "preprocess.R".
#
# Inputs:
# - expression.qual.txt: Quality-adjusted gene expression matrix.
# - methylation.qual.txt: Quality-adjusted methylation matrix (M-values).
# - clinical_data.txt: Clinical metadata containing sample information.
#
# Outputs:
# - expression_fs.txt: Filtered gene expression matrix (by MAD).
# - methylation_fs.txt: Filtered methylation matrix (by MAD).
#

# Required environment: 'preprocess.yml' (containing prerequiste libraries).


# Load required libraries
library(tidyverse)

# Load custom functions for data preprocessing
source("functions/preprocess/preprocess.R")

config <- read.csv("config/config.tsv",
    sep = "\t",
    header = T
)
# Load configuration file for setting working and data directories
keys <- config$key
values <- config$value
names(values) <- keys
# Set working and data directories from the config file
wkdir <- values["wkdir"]
datadir <- values["datadir"]


####################################################
# Load Quality-Adjusted Data
####################################################

# Load the quality-adjusted gene expression matrix

expression <- read.table(paste0(wkdir,
"/preprocess/transcriptome/expression.qual.txt"),
sep = '\t',
header = T,
row.names = 1,
check.names = F
)

# Load the quality-adjusted methylation matrix (M-values)
methylation <- read.table(paste0(wkdir,
"/preprocess/methylome/methylation.qual.txt"),
sep = '\t',
header = T,
row.names = 1,
check.names = F

)

# Load clinical data containing sample metadata
clinical.data <- read.table(paste0(datadir, "/clinical_data.txt"),
sep = '\t',
header = T,
row.names = 1
)

####################################################
# Create Directories for Filtered Data Output
####################################################

# Create a directory for storing the MAD-filtered data if it doesn't exist

if (dir.exists(paste0(wkdir, "/preprocess/mad")) == F) {

    dir.create(paste0(wkdir, "/preprocess/mad"))

}
# Set working directory to the newly created "mad" directory
setwd(paste0(wkdir, "/preprocess/mad"))



####################################################
# Filter Features by Mean Absolute Deviation (MAD)
####################################################

# Define a list of tissue types to consider for filtering
tissues <- c("COLORECTAL", "PANCREAS", "LUNG", "GASTRIC", "SMALL INTESTINE")


# Apply filtering on gene expression data, selecting the most variable genes (by MAD)
# across the defined tissue types
expression_fs <- filter_tissue(expression, clinical.data,
tissues = tissues, flag = "expression")

# Apply filtering on methylation data, selecting the most variable probes (by MAD)
# across the defined tissue types
methylation_fs <- filter_tissue(methylation, clinical.data,
tissues = tissues, flag = "methylation")



####################################################
# Save Filtered Data
####################################################

# Write the filtered gene expression matrix to a new fil
write.table(expression_fs,
paste0(wkdir,
"/preprocess/transcriptome/expression_fs.txt"),
sep = "\t",
row.names = T,
col.names = NA
)


# Write the filtered methylation matrix to a new file
write.table(methylation_fs,
paste0(wkdir,
"/preprocess/methylome/methylation_fs.txt"),
sep = "\t",
row.names = T,
col.names = NA
)

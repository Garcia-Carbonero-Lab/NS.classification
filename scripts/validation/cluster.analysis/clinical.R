
###############################################################################
# Script: Clinical Data Analysis for Subtype Classification
#
# Description:
#     This script analyzes clinical data by comparing feature distributions
#     across different subtypes. Two tables are created: one with NA values
#     and one without.
#
# Inputs:
#     - Validation dataset: Base_validation.subtype.tsv
#     - Configuration file: config.tsv
#
# Outputs:
#     - Clinical tables comparing feature distributions across subtypes:
#       - Clinical.txt (without NA values)
#       - Clinical.na.txt (with NA values)
#
# Required environment: `clinical.yml`
###############################################################################

# Load required libraries
library(ggplot2)
library(reshape2)
library(ggpubr)
library(rstatix)
library(tidyverse)

# Load custom clinical functions
source("functions/clustering/clinical.functions.R")

# Read configuration file to extract working and data directories
config <- read.csv("config/config.tsv", sep = "\t", header = TRUE)

# Extract key-value pairs for configuration settings
keys <- config$key
values <- config$value
names(values) <- keys

# Define working directory and data directory
wkdir <- values["wkdir"]
datadir <- values["datadir"]

# Create output directories if they do not exist
if (!dir.exists(paste0(wkdir, "/validation/clustering.analysis"))) {
    dir.create(paste0(wkdir, "/validation/clustering.analysis"))
}

if (!dir.exists(paste0(wkdir, "/validation/clustering.analysis/clinical"))) {
    dir.create(paste0(wkdir, "/validation/clustering.analysis/clinical"))
}

# Define the output directory for clinical analysis
outdir <- paste0(wkdir, "/validation/clustering.analysis/clinical")

# Load the validation dataset
base.validation <- read.table(paste0(wkdir, "/validation/base/Base_validation.subtype.tsv"),
                              sep = '\t', header = TRUE, row.names = 1)

# Modify the Stage column to create a binary Stage.IV column
base.validation$Stage.IV <- base.validation$Stage
base.validation$Stage.IV[base.validation$Stage %in% c("I","II","III")] <- "No"
base.validation$Stage.IV[base.validation$Stage == "IV"] <- "Yes"

# Define the clinical features to analyze
features <- c("primary.tumor", "Gender", "Grade", "Stage.IV", "CIMP.status", "Clusters", "tissue")

# Create a clinical table without NA values, grouped by subtype
table <- clinical.table(data = base.validation, features = features, group = "Subtype.validation")

# Create a clinical table with NA values, grouped by subtype
table_na <- clinical.table(data = base.validation, features = features, group = "Subtype.validation", na = TRUE)

# Write the clinical tables to output files
write.table(table, paste0(outdir, "/Clinical.txt"), sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
write.table(table_na, paste0(outdir, "/Clinical.na.txt"), sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

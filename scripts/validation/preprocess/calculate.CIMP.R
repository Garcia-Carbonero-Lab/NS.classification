###########################################################
# Script: CIMP Status Calculation and Integration into Clinical Data
#
# Description:
#     This script calculates the CIMP (CpG Island Methylator Phenotype) status 
#     for samples based on methylation data and integrates this information into 
#     the clinical dataset for further analysis. The result is a clinical table 
#     that includes the CIMP status for each sample.
#
# Requirements:
#     - Preprocessed methylation data and clinical data from multiple studies.
#     - A predefined function to calculate CIMP status (`addCIMPstatus()`).
#
# Inputs:
#     - Methylation matrix (`methylation_all.tsv`).
#     - Clinical dataset (`Base_validation.tsv`).
#
# Outputs:
#     - Updated clinical dataset with CIMP status (`Base_validation_cimp.tsv`).
#
# Required environment: `preprocess` (contains necessary libraries and functions for processing).
###########################################################



############### Load Configuration File ###############
# Load configuration file to retrieve working and data directory paths
config <- read.csv("config/config.tsv",
    sep = "\t",
    header = T
)

# Extract key-value pairs for configuration settings
keys <- config$key
values <- config$value
names(values) <- keys

# Define working directory and data directory from configuration file
wkdir <- values["wkdir"]
datadir <- values["datadir"]

############### Load Required Functions ###############
# Source the methylation processing functions required for CIMP calculation
source("functions/preprocess/methylation.R")



############### Load Clinical and Methylation Data ###############
# Load clinical data table (Base_validation.tsv)
clinical.data <- read.table(paste0(wkdir, "/validation/base/Base_validation.tsv"),
sep = '\t',
header = T,
row.names = 1
)

# Load the combined methylation matrix (methylation_all.tsv)
methylation <- read.table(paste0(wkdir, "/validation/whole.matrix/methylation_all.tsv"),
    sep = "\t",
    row.names = 1,
    header = T,
    check.names = F
)

############### Preprocessing of Methylation Data ###############
# Convert all methylation data values to numeric (ensure proper data type for calculation)
methylation <- apply(methylation, c(1,2), as.numeric)

# Remove rows with incomplete cases (i.e., missing data) from the methylation matrix
methylation <- methylation[complete.cases(methylation),]


############### CIMP Status Calculation ###############
# Calculate CIMP status for the samples using the predefined `addCIMPstatus` function
CIMP <- addCIMPstatus(methylation)



############### Integrate CIMP Status into Clinical Data ###############
# Merge the CIMP status with the clinical data
# The merge is done by matching row names (sample IDs) from both datasets

clinical.data <- merge(clinical.data, CIMP,
    by.x = "row.names",
    by.y = "row.names",
    all.x = T
)

############### Save Updated Clinical Data with CIMP Status ###############
# Write the updated clinical data with CIMP status to a new TSV file
write.table(clinical.data, paste0(wkdir,"/validation/base/Base_validation_cimp.tsv"),
    sep = "\t", row.names = F, col.names = T
)

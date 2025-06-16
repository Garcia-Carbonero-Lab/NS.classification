####################################################
# Script: Processing and Normalizing Lakis PANET Methylation Array Data
#
# Description:
#     This script processes methylation array data from the Lakis PANET dataset. It reads methylation 
#     values (beta values), reshapes the data, filters patient metadata, and calculates M-values for 
#     downstream analysis. The script outputs both the original beta matrix and the transformed M-value matrix.
#
# Requirements:
#     - GEOquery package for GEO dataset handling (if needed).
#     - Configuration file (`config/config.tsv`) for specifying working and data directories.
#     - Methylation data file (`meth_array.PAEN-AU.tsv`) and patient metadata file (`base_Lakis_prev.tsv`).
#     - Custom functions from `validation.functions.R` for any additional processing (if needed).
####################################################
# Required environment: 'geo.yml' (containing prerequiste libraries).

# Load required libraries for handling methylation data and data manipulation
library(GEOquery)
library(tidyverse)
# library(mygene)  # This package is commented out but might be needed for gene annotation

# Load custom validation functions (if required)
source("functions/validation/validation.functions.R")

############### Configuration and Directory Setup ###############
# Read the configuration file (config.tsv) to retrieve directory paths
config <- read.csv("config/config.tsv", sep = "\t", header = TRUE)

# Assign configuration keys and values
keys <- config$key
values <- config$value
names(values) <- keys

# Define working directory and data directory from the configuration file
wkdir <- values["wkdir"]
datadir <- values["datadir"]

# Check if output directories exist, create them if not
if (!dir.exists(paste0(wkdir, "/validation"))) {
    dir.create(paste0(wkdir, "/validation"))
}

if (!dir.exists(paste0(wkdir, "/validation/Lakis"))) {
    dir.create(paste0(wkdir, "/validation/Lakis"))
}

# Set output directory for validation results
outdir <- paste0(wkdir, "/validation/Lakis")

############### Methylation Data Loading and Reshaping ###############
# Load methylation data (beta values) from the specified file
methylation <- read.table(paste0(datadir, "/external_data/Lakis_PANET/meth_array.PAEN-AU.tsv"),
                          sep = "\t", header = TRUE)

# Reshape the methylation data: organize probes as rows and samples as columns
methylation <- methylation %>%
    select(icgc_sample_id, probe_id, methylation_value) %>%  # Select relevant columns
    pivot_wider(names_from = icgc_sample_id, values_from = methylation_value)  # Pivot to wide format

############### Patient Annotation and Data Filtering ###############
# Convert the reshaped methylation matrix into a data frame and set probe IDs as rownames
methylation <- as.data.frame(methylation)
rownames(methylation) <- methylation[, 1]  # Set the first column (probe_id) as rownames
methylation <- methylation[, -1]  # Remove the probe_id column

# Load patient metadata file
base <- read.table(paste0(wkdir, "/validation/Lakis/base_Lakis_prev.tsv"),
                   sep = "\t", header = TRUE, row.names = 1)

# Subset the metadata to match the samples in the methylation matrix
base <- base[colnames(methylation),]

# Filter out specific differentiation types and primary tumor types
base <- base[!base$Differenciation %in% c("Poorly differentiated endocrine carcinoma", 
                                          "Mixed ductal endocrine carcinoma"),]

base <- base[!base$primary.tumor %in% c("Normal - tissue adjacent to primary"),]

# Subset the methylation matrix based on the filtered metadata
methylation <- methylation[, rownames(base)]

# Rename the columns of the methylation matrix to reflect patient IDs from the metadata
colnames(methylation) <- base$id_patient

############### Saving the Beta Value Matrix ###############
# Write the processed beta value methylation matrix to a file
write.table(methylation,
            paste0(outdir, "/Methylation_Beta_Lakis.tsv"),
            sep = "\t", row.names = TRUE, col.names = NA)

############### M-value Calculation ###############
# M-value transformation: adjust extreme beta values and log-transform
methylation[methylation <= 0.001] <- 0.001  # Set a lower bound for beta values
methylation[methylation >= 0.999] <- 0.999  # Set an upper bound for beta values

# Calculate M-values using the formula M = log2(beta / (1 - beta))
methylation <- log((methylation / (1 - methylation)), 2)

############### Saving the M-value Matrix ###############
# Write the M-value transformed methylation matrix to a file
write.table(methylation,
            paste0(outdir, "/Methylation_M_Lakis.tsv"),
            sep = "\t", row.names = TRUE, col.names = NA)
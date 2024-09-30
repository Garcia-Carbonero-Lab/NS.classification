# Script: Methylation Data Preprocessing and Normalization
#
# Purpose:
#   This script processes raw methylation data using the ChAMP package and 
#   performs quality control, filtering, and normalization (BMIQ method). 
#   The resulting normalized beta and M-values are written to output files.
#
# Prerequisites:
#   - R environment with libraries: ChAMP, tidyverse.
#   - Configuration file: 'config/config.tsv' containing paths for directories.
#   - Methylation data (.IDAT files) located in the specified directory.
#   - Clinical data in 'clinical_data.txt'.
#   - Custom functions from 'methylation.process.R' for preprocessing.
#   - Requires 'methylation_m.yml' environment to be set up.
#
# Required environment: 'methylation_m.yml' (containing prerequiste libraries).


# Load required libraries
library(ChAMP)
library(tidyverse)

# Load custom preprocessing functions for methylation data
source("functions/preprocess/methylation.R")

# Load configuration file containing key-value pairs for directory paths
config <- read.csv("config/config.tsv",
    sep = "\t",
    header = T
)

# Extract key paths from the configuration file
keys <- config$key
values <- config$value
names(values) <- keys

#Set working directory and data directory paths
wkdir <- values["wkdir"]
datadir <- values["datadir"]


# Load clinical data from the specified file
clinical.data <- read.table(paste0(datadir, "/clinical_data.txt"),
    sep = "\t",
    header = T
)

# Import methylation data using ChAMP's champ.import function
# Methylation data should be in the specified 'methylome/raw' directory
myImport <- champ.import(paste0(datadir, "/methylome/raw"),
    arraytype = "EPIC"
)


# Calculate methylation quality for each sample using the custom 'qual.met' function
qual <- qual.met(myImport)

# Merge the quality metrics with the clinical data based on the 'CLINICAL_TRIAL_CODE'
clinical.data <- merge(clinical.data,
    qual,
    by.x = "CLINICAL_TRIAL_CODE",
    by.y = "row.names"
)

# Filter methylation data based on default quality control parameters

myLoad <- champ.filter(
    beta = myImport$beta, M = myImport$M,
    SampleCutoff = 0.15, pd = myImport$pd,
    arraytype = "EPIC",
    detP = myImport$detP
)

# Create the necessary directories for preprocessing if they don't already exist
if (dir.exists(paste0(wkdir, "/preprocess")) == F) {

    dir.create(paste0(wkdir, "/preprocess"))

}

if (dir.exists(paste0(wkdir, "/preprocess/methylome")) == F) {

    dir.create(paste0(wkdir, "/preprocess/methylome"))

}

# Set the working directory for preprocessing output
setwd(paste0(wkdir, "/preprocess/methylome"))

# Normalizaton don't support names with spaces


# Perform BMIQ normalization on the filtered beta values

myNormB <- champ.norm(
    beta = myLoad$beta, arraytype = "EPIC",
    cores = 10,
    method = "BMIQ"
)


# Remove rows with missing values from the normalized beta values
myNormB <- myNormB[complete.cases(myNormB), ]
myNormB <- apply(myNormB, c(1,2), as.numeric)


# Write the normalized beta values to a file
write.table(myNormB, "myNormB.txt", sep = "\t", row.names = T, col.names = NA)


# Update the clinical data file with new quality metrics and save it
write.table(clinical.data, paste0(datadir,"/clinical_data.txt"),
    sep = "\t", row.names = F, col.names = T
)


# Convert normalized beta values to M-values
# M-values are calculated as log2(beta / (1 - beta)) and bounded between 0.001 and 0.999
myNormM <- myNormB
myNormM[myNormM <= 0.001] <- 0.001
myNormM[myNormM >= 0.999] <- 0.999
myNormM <- log((myNormM / (1 - myNormM)), 2)


# Write the normalized M-values to a file
write.table(myNormM,
    "myNormM.txt",
    sep = "\t",
    col.names = T,
    row.names = T,
    quote = F
)

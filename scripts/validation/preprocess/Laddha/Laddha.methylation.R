####################################################
# Script: Processing Methylation Data for Laddha LUNET Dataset
#
# Description:
#     This script processes methylation data from the Laddha LUNET dataset. 
#     It generates the sample sheet from IDAT files, processes the methylation 
#     array data, and saves the results for downstream analysis.
#
# Requirements:
#     - `met.arrays` environment.
#     - Configuration file (`config/config.tsv`) for directory paths.
#     - Methylation data (IDAT files).
#     - Custom validation and data processing functions from `validation.functions.R`.
#     - IlluminaHumanMethylation450k annotation and manifest packages.
#
# Inputs:
#     - IDAT files from the Laddha LUNET dataset.
#
# Outputs:
#     - Processed methylation matrix in the output directory.
#
# Required environment: met.arrays.yml (containing relevant methylation packages).
####################################################


############### Load Required Libraries ###############
# Load libraries for handling methylation array data and data processing
library(minfi)
library(ENmix)
library(wateRmelon)
library(MASS)
library(broom)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)


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

# Load custom validation functions
source("functions/validation/validation.functions.R")

############### Configuration and Directory Setup ###############
# Create output directories if they don't exist

if (dir.exists(paste0(wkdir, "/validation")) == F) {

    dir.create(paste0(wkdir, "/validation"))

}

if (dir.exists(paste0(wkdir, "/validation/Laddha")) == F) {

    dir.create(paste0(wkdir, "/validation/Laddha"))

}

# Define the output directory for results
outdir <- paste0(wkdir, "/validation/Laddha")



############### Methylation Data Processing ###############
# Define path to the IDAT files
path <- paste0(datadir, "/external_data/Laddha_LUNET/IDAT")

# Generate a sample sheet for the methylation dataset based on IDAT files
generate.sample.sheet(path)


# Process the methylation data using the custom function from `validation.functions.R`
methylation_matrix(paste0(datadir, "/external_data/Laddha_LUNET"),
arraytype = "450K",
outdir)

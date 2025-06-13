####################################################
# Patte Methylation Data Processing Script
#
# Description:
#     This script processes Illumina EPIC methylation array data from the Alcala 
#     dataset. It generates a sample sheet, organizes IDAT files, performs normalization, 
#     and writes the final methylation matrix in M-value format for downstream analysis.
#
# Requirements:
#     - Libraries for handling methylation data: minfi, ENmix, wateRmelon, MASS, etc.
#     - Custom functions for methylation processing.
#     - Methylation data files in IDAT format.
#     - Configuration file (config.tsv) to define working and data directories.
# 
# Required environment: 'met.arrays.yml' (containing prerequisite libraries).
####################################################

#this iscript uses environment met.arrays.yml
library(minfi)
library(ENmix)
library(wateRmelon)
library(MASS)
library(broom)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)



####################################################
# Load configuration and set directories
####################################################
# Open the configuration file to set paths for working and data directories
config <- read.csv("config/config.tsv",
sep = "\t",
header = T)
# Extract key-value pairs for configuration settings
keys <- config$key
values <- config$value
names(values) <- keys

# Set working directory and data directory based on configuration values
wkdir <- values["wkdir"]
datadir <- values["datadir"]
source("functions/validation/validation.functions.R")

####################################################
# Create Directories for Results
####################################################
# Check if the validation directory exists, create it if it doesn't
if (dir.exists(paste0(wkdir, "/validation")) == F) {

    dir.create(paste0(wkdir, "/validation"))

}

# Check if the Patte validation directory exists, create it if it doesn't
if (dir.exists(paste0(wkdir, "/validation/Patte")) == F) {

    dir.create(paste0(wkdir, "/validation/Patte"))

}

# Define output directory for storing Patte methylation data
outdir <- paste0(wkdir, "/validation/Patte")

# Generate the sample sheet for Patte IDAT files
generate.sample.sheet(path)

# Load the generated sample sheet

targets <- read.metharray.sheet(path)


############### Methylation Matrix Generation ###############
# Function to generate methylation data matrix from IDAT files


methylation_matrix(paste0(datadir, "/external_data/SiNETs_Patte/IDAT"),
arraytype = "EPIC",
outdir)

# Load the M-value methylation matrix generated from IDAT files


M.value <- read.table(paste0(outdir,
"/M_matrix.csv"), sep=",", row.names = 1, header = T, check.names = F)  


############### Annotations and Final Output ###############
# Load sample annotations  to update column names


anotation <- read.csv(paste0(datadir, "/external_data/SiNETs_Patte/METADATA_Samples_TNEi_for_paper_parcial.csv"),
sep = "\t",
row.names = 1,
header = T,
check.names = FALSE
)

targets$met_id <- paste(targets$Slide,targets$Array, sep = "_")
rownames(targets) <- targets$met_id

anotation <- anotation[anotation$methyl.ID %in% targets$met_id,]

targets  <- targets[anotation$methyl.ID,]
anotation$Code_methylation <- targets$Sample_Name[match(anotation$methyl.ID, targets$met_id)]


anotation <- anotation[match(colnames(M.value), anotation$Code_methylation),]


colnames(M.value) <- rownames(anotation)

# Write the processed M-value methylation matrix to a TSV file
write.table(M.value,
paste0(outdir, "/Patte_M.tsv"),
sep = "\t",
row.names = T,
col.names = NA)
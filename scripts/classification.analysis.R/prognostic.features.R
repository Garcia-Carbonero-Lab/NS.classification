##################################################
# Script: Prognostic Feature Identification using Gene Expression and Methylation Data
#
# Description:
#      This script identifies prognostic features related to cancer survival 
#      using Cox proportional hazards models. It analyzes gene expression 
#      and methylation data to find features that are associated with survival 
#      and annotates methylation probes.
#
# Requirements:
#      - Clinical data file (`clinical_data.txt`) containing patient information.
#      - Gene expression data (`expression.qual.txt`).
#      - Methylation data (`Methylation.NS.txt`).
#      - Methylation annotation file (`infinium-methylationepic-v-1-0-b5-manifest-file.csv`).
###################################################
# Required environment: survival.yml

############### Library Setup ###############
# Load necessary libraries for survival analysis, plotting, and model fitting

library(survival)
library(survminer)
library(forestplot)
library(tidyverse)
library(aod)
library(forestmodel)
library(rms)


# Load custom functions for survival analysis and general functions
source("functions/clustering/survival.functions.R")
source("functions/main.functions.R")



############### Configuration File Loading ###############
# Read the configuration file (config.tsv) to obtain directory paths and other settings


config <- read.csv("config/config.tsv",
sep = "\t",
header = T)

# Extract key-value pairs from the config file for directory setup
keys <- config$key
values <- config$value
names(values) <- keys


# Set the working directory and data directory paths based on the config file
wkdir <- values["wkdir"]
datadir <- values["datadir"]


############### Directory Setup ###############
# Check if the directory for prognostic feature results exists, if not create it

if (dir.exists(paste0(wkdir, "/clustering/prognostic.features")) == F) {

    dir.create(paste0(wkdir, "/clustering/prognostic.features"))

}

# Set the output directory for prognostic feature results

outdir <- paste0(wkdir, "/clustering/prognostic.features")


############### Data Loading ###############
# Load clinical data from the specified file (clinical_data.txt)

base <- read.csv(paste0(datadir, "/clinical_data.txt"),
                 sep = "\t",
                 row.names = 1,
                 header = T)

# Load gene expression data from the specified file
expression <- read.table(paste0(wkdir,
"/preprocess/transcriptome/expression.qual.txt"),
    sep = "\t",
    header = T,
    row.names = 1,
    check.names = F)

# Load DEG filtered tables for different subtypes
NS1.filter <- read.table(
paste0(wkdir, "/clustering/deg/Table_NS1_filter.txt"),
sep = "\t",
row.names = 1,
header = T)

NS2.filter <- read.table(
paste0(wkdir, "/clustering/deg/Table_NS2_filter.txt"),
sep = "\t",
row.names = 1,
header = T)

NS3.filter <- read.table(
paste0(wkdir, "/clustering/deg/Table_NS3_filter.txt"),
sep = "\t",
row.names = 1,
header = T)

# Combine unique gene names from the three filtered DEG tables
genes <- unique(c(rownames(NS1.filter),rownames(NS2.filter),rownames(NS3.filter)))

# Filter the expression data to include only the selected genes
expression <- expression[genes,]

# Replace hyphens with underscores in gene row names
rownames(expression) <- gsub("-","_",row.names(expression))

############### Cox Proportional Hazards Model for Gene Expression ###############
# Perform Cox proportional hazards model with multiple covariates on gene expression data

multiple_cov <- cox.multiple(base = base,
                             data = expression,
                             OS = "OS.time",          # Overall survival time column
                             OS_event = "EXITUS",     # Event status column (0 = alive, 1 = dead)
                             covariates = c("STAGE_IV", "PRIMARY_TUMOR", "GENDER", "AGE"),  # Covariates
                             outdir = outdir,         # Output directory for results
                             flag = "Cox_Multiple_covariates",   # Flag to append to output file names
                             th = 0.05)                 # Threshold for variable selection

################################################################################



############### Methylation Data Analysis ###############
# Load methylation data from the specified file (Methylation.NS.txt)
#with differentially methylated probes (DMP)

methylation <- read.table(paste0(wkdir,
"/clustering/dmp/Methylation.NS.txt"),
sep = "\t",
row.names= 1,
header = T,
check.names = F)

# Load methylation annotation data (Infinium Methylation EPIC annotation file)

epic <- read.csv(paste0(datadir,
"/methylome/anotation/infinium-methylationepic-v-1-0-b5-manifest-file.csv"),
                       sep = ",",
                       header = T,
                       skip = 7)
    

# Transform methylation annotation matrix using custom function `create.met.anot`
met_anot <- create.met.anot(epic)

# Replace hyphens with underscores in methylation row names
rownames(methylation) <- gsub("-","_",row.names(methylation))

# Perform Cox proportional hazards model with multiple covariates on methylation data
multiple_cov <- cox.multiple(base = base,
                             data = methylation,
                             OS = "OS.time",          # Overall survival time column
                             OS_event = "EXITUS",     # Event status column (0 = alive, 1 = dead)
                             covariates = c("STAGE_IV", "PRIMARY_TUMOR", "GENDER", "AGE"),  # Covariates
                             outdir = outdir,         # Output directory for results
                             flag = "Cox_Multiple_covariates_methylation",  # Flag for methylation results
                             th = 0.05)                 # Threshold for variable selection


############### Prognostic Feature Extraction ###############
# Load the Cox results file for methylation features and filter for significant prognostic features (p < 0.05)
prognostic <- read.table(
paste0(wkdir, "/clustering/prognostic.features/Cox_Multiple_covariates_methylation.txt"),
sep = "\t",
row.names = 1,
header = T)

# Keep only rows where the p-value is less than 0.05
prognostic <- prognostic[prognostic$p < 0.05,]

# Merge the methylation annotation data with the prognostic features based on the probe ID

prognostic <- merge(met_anot, prognostic, by.x = "ID", by.y = "row.names")

# Save the annotated prognostic methylation features to a new file
write.table(prognostic,
paste0(outdir,"/Methylation.prognostic.anotated.txt"),
sep = "\t",
row.names= F,
col.names= T)

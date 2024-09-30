# Script: Transcriptome Preprocessing
# 
# Purpose:
#   This script processes transcriptome data from Affymetrix Clariom S Human arrays.
#   It reads the configuration file to set up directories and paths, loads required 
#   functions, and normalizes the gene expression data. The normalized expression matrix 
#   is then saved in a designated output directory. The environment 'arrays_m' should be 
#   preloaded before running this script.
# 
# Prerequisites:
#   - R environment with libraries: limma, makecdfenv, affy, vsn.
#   - Config file: 'config/config.tsv' containing key-value pairs for directories and settings.
#   - Affymetrix Clariom S Human array files in .CEL format and the CDF file.
#   - Custom functions from 'functions/preprocess/expression.R' to perform normalization.
# 
# Required environment: 'arrays_m' (containing library prerequisites).

# Load required libraries
library(limma)
library(makecdfenv)
library(affy)
library(vsn)



#Open config file
config <- read.csv("config/config.tsv",
sep = "\t",
header = T)
# Extract directory paths from the configuration file
keys <- config$key
values <- config$value
names(values) <- keys


# Set working directory and data directory paths
wkdir <- values["wkdir"]
datadir <- values["datadir"]


# Load preprocessing function for expression data
source("functions/preprocess/expression.R")

# Define target file and input directory paths
targetfile <- paste0(datadir, "/transcriptome/raw/targets.txt")
indir <- paste0(datadir, "/transcriptome/raw")


# Create necessary directories if they do not exist
if (dir.exists(paste0(wkdir, "/preprocess")) == F) {

    dir.create(paste0(wkdir, "/preprocess"))

}

if (dir.exists(paste0(wkdir, "/preprocess/transcriptome")) == F) {

    dir.create(paste0(wkdir, "/preprocess/transcriptome"))

}
# Define output directory for processed data
outdir <- paste0(wkdir, "/preprocess/transcriptome")

# Load the Clariom S Human CDF file for probe mapping
clariomshumancdf <- make.cdf.env("Clariom_S_Human.r1.Gene.CDF", indir)

# Call the normalization function to process and save the expression data
norm_exp(targetfile = targetfile # Path to the target file (.CEL file locations)
         , outputpath = outdir   # Output directory to save normalized data
         , flag = "RMA"          # Custom flag for the output files
         , inpath = indir)      # Path to input files and annotations
 
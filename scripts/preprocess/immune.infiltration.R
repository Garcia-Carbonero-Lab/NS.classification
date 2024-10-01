# This script processes gene expression data to compute immune cell infiltration scores
# using the MCPcounter method. The script reads a preprocessed expression matrix, 
# calculates MCPcounter scores, and saves the output in a designated directory.
#
# Requirements:
# - The script uses the MCPcounter package to estimate immune cell populations.
# - A quality-adjusted expression matrix must be available.
# - Configurations for directories are specified in "config/config.tsv".
#
# Inputs:
# - expression.qual.txt: Quality-adjusted gene expression matrix.
#
# Outputs:
# - MCPcounter.tsv: Table of estimated immune cell scores.
#
# Create directories if they don't exist and save results in appropriate locations.


# Required environment: 'preprocess.yml' (containing prerequiste libraries).


# Load required library for immune cell estimation

library(MCPcounter)

# Load configuration file to set working and data directories
config <- read.csv("config/config.tsv",
    sep = "\t",
    header = T
)
keys <- config$key
values <- config$value
names(values) <- keys


# Set working and data directories from the config file
wkdir <- values["wkdir"]
datadir <- values["datadir"]


####################################################
# Load Data
####################################################

# Load the quality-adjusted gene expression matrix
expression <- read.table(paste0(wkdir,
"/preprocess/transcriptome/expression.qual.txt"),
sep = '\t',
header = T,
row.names = 1,
check.names = F)


####################################################
# Create Directory for MCPcounter Output
####################################################

# Create a directory for storing the MCPcounter results if it doesn't exist
if (dir.exists(paste0(wkdir, "/immune")) == F) {

    dir.create(paste0(wkdir, "/immune"))

}


# Calculate immune cell scores using the MCPcounter algorithm
# Input expression data is assumed to use HUGO gene symbols

mcp <- MCPcounter.estimate(expression, featuresType = "HUGO_symbols")

# Adjust row names to replace spaces with dots (if any exist in the cell type labels)
rownames(mcp) <- gsub(" ",".",rownames(mcp))


####################################################
# Save MCPcounter Results
####################################################

# Write the MCPcounter output to a tab-separated file
write.table(mcp,
paste0(wkdir, "/immune/MCPcounter.tsv"),
sep = "\t",
row.names = T,
col.names = NA)


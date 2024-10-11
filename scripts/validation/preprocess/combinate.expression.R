###########################################################
# Script: Expression Data Combination for Multiple Studies
#
# Description:
#     This script prepares a combined expression dataset from multiple studies
#     (Alcala, Laddha, Lakis, Karpathakis) by loading, filtering, and merging
#     their expression matrices. The output is a unified expression matrix
#     containing samples common across the studies, with genes filtered for use 
#     in a MOFA (Multi-Omics Factor Analysis) model.
#
# Requirements:
#     - R environment with required packages: tidyverse.
#     - Expression data from multiple studies in tsv format.
#     - A list of genes to filter for MOFA model.
#
# Inputs:
#     - Expression data from the studies (Alcala, Lakis, Karpathakis).
#     - Base clinical data to filter samples.
#     - Pre-processed list of selected genes for MOFA filtering.
#
# Outputs:
#     - Combined expression matrix for all studies (`expression.combination.all.tsv`).
#     - Filtered expression matrix for MOFA model (`expresssion.fs.txt`).
#
# Required environment: preprocess
###########################################################


############### Load Required Libraries ###############
library(tidyverse)

############### Configuration and Directory Setup ###############
# Load configuration file to get working and data directory paths
config <- read.csv("config/config.tsv",
    sep = "\t",
    header = T
)

# Extract key-value pairs from the configuration file
keys <- config$key
values <- config$value
names(values) <- keys
# Define working directory and data directory from configuration file
wkdir <- values["wkdir"]
datadir <- values["datadir"]

# Create output directory for combined expression matrices if it doesn't exist
if (dir.exists(paste0(wkdir, "/validation/whole.matrix")) == F) {

    dir.create(paste0(wkdir, "/validation/whole.matrix"))

}

# Define output directory for writing results
outdir <- paste0(wkdir, "/validation/whole.matrix")


############### Load Base Clinical Metadata ###############
# Load the base clinical metadata, which contains information about the samples.

base <- read.table(paste0(wkdir, "/validation/base/Base_validation.tsv"),
sep = "\t",
row.names = 1,
header = T)

############### Load Expression Data from Multiple Studies ###############

# Load Alcala & Laddha expression data
exp_Alcala <- read.table(paste0(wkdir,"/validation/Alcala/Alcala.expression.tsv"),
sep = "\t",
row.names = 1,
header = T)

# Filter Alcala expression data to include only samples that exist in the base clinical metadata

exp_Alcala <- exp_Alcala[,colnames(exp_Alcala) %in% rownames(base)]

# Load Lakis expression data
exp_Lakis<- read.table(paste0(wkdir,"/validation/Lakis/matrix_Lakis_expression.txt"),
sep = "\t",
row.names = 1,
header = T)
# Filter Lakis expression data to include only samples present in the base clinical metadata
exp_Lakis <- exp_Lakis[,colnames(exp_Lakis) %in% rownames(base)]

# Load Karpathakis expression data
exp_Karpathakis<- read.table(paste0(wkdir,"/validation/Karpathakis/Karpathakis.expression.tsv"),
sep = "\t",
row.names = 1,
header = T,
check.names = F)
# Log-transform Karpathakis expression data to stabilize variance (log(x + 1))
exp_Karpathakis <- log(exp_Karpathakis +1)

# Filter Karpathakis expression data to include only samples present in the base clinical metadata
exp_Karpathakis <- exp_Karpathakis[,colnames(exp_Karpathakis) %in% rownames(base)]



############### Merge Expression Data from Multiple Studies ###############

# Merge Lakis and Alcala expression data by matching on gene symbols
combination <- merge(exp_Lakis, exp_Alcala,by.x = "row.names",by.y = "row.names")

# Set row names as the gene identifiers and remove redundant first column
rownames(combination) <- combination[,1]
combination <- combination[,-1]

# Merge the combined expression matrix (Lakis & Alcala) with Karpathakis data
combination <- merge(combination, exp_Karpathakis,by.x = "row.names",by.y = "row.names")

# Set row names again and remove the redundant first column
rownames(combination) <- combination[,1]
combination <- combination[,-1]

# Write the combined expression matrix for all studies to a file
write.table(combination, paste0(outdir,"/expression.combination.all.tsv"),
sep = "\t",
row.names = T,
col.names= NA)


############### Filter Genes Used in MOFA Model ###############

# Load the list of genes used in the MOFA model

exp <- read.table(paste0(wkdir,
"/preprocess/transcriptome/expression_fs.txt"),
    sep = "\t",
    header = T,
    row.names = 1,
    check.names = F
)

# Filter the combined expression matrix to include only the genes used in the MOFA model
combination <- combination[rownames(exp),]

# Write the filtered expression matrix for MOFA model to a file
write.table(combination,
paste0(outdir, "/expresssion.fs.txt"),
sep = "\t",
row.names = T,
col.names = NA)

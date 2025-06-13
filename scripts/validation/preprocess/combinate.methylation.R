###########################################################
# Script: Methylation Data Combination for Multiple Studies
#
# Description:
#     This script processes and combines methylation data from multiple studies
#     (Alcala, Laddha, Lakis, Patte) by loading, filtering, and merging
#     their methylation matrices. The output is a unified methylation matrix
#     containing samples common across the studies, with probes filtered for use
#     in a MOFA (Multi-Omics Factor Analysis) model.
#
# Requirements:
#     - R environment with required packages: tidyverse.
#     - Methylation data from multiple studies in tsv format.
#     - A list of methylation probes to filter for the MOFA model.
#
# Inputs:
#     - Methylation data from the studies (Alcala, Laddha, Lakis, Patte).
#     - Base clinical data to filter samples.
#     - Pre-processed list of selected probes for MOFA filtering.
#
# Outputs:
#     - Combined methylation matrix for all studies (`methylation_all.tsv`).
#     - Filtered methylation matrix for MOFA model (`methylation.fs.txt`).
#
# Required environment: preprocess
###########################################################

############### Load Required Libraries ###############
library(tidyverse)

##############################################33

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

# Create output directory for combined methylation matrices if it doesn't exist
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

############### Load Methylation Data from Multiple Studies ###############


# Load Alcala methylation data
met_Alcala <- read.table(paste0(wkdir,"/validation/Alcala/Alcala_M.tsv"),
sep = "\t",
row.names = 1,
header = T,
check.names = F)

# Filter Alcala methylation data to include only samples that exist in the base clinical metadata
colnames(met_Alcala) %in% rownames(base)
met_Alcala <- met_Alcala[,colnames(met_Alcala) %in% rownames(base)]

# Load Lakis methylation data
met_Lakis<- read.table(paste0(wkdir,"/validation/Lakis/Methylation_M_Lakis.tsv"),
sep = "\t",
row.names = 1,
header = T,
check.names = F)

# Filter Lakis methylation data to include only samples present in the base clinical metadata
met_Lakis <- met_Lakis[,colnames(met_Lakis) %in% rownames(base)]
# Load Laddha methylation data
met_Laddha<- read.table(paste0(wkdir,"/validation/Laddha/M_matrix_noinf.csv"),
sep = ",",
row.names = 1,
header = T)

# Match Laddha samples in the base clinical data based on 'geo' identifier
names <- rownames(base)[base$geo %in% colnames(met_Laddha)]
base_Laddha <- base[names,]
base_Laddha$names <- row.names(base_Laddha)
rownames(base_Laddha) <- base_Laddha$geo
base_Laddha <- base_Laddha[colnames(met_Laddha),]

# Update column names of Laddha methylation data with sample names from the base data
colnames(met_Laddha) <- base_Laddha$names

# Load Patte methylation data
met_Patte<- read.table(paste0(wkdir,"/validation/Patte/Patte_M.tsv"),
sep = "\t",
row.names = 1,
header = T)

# Filter Patte methylation data to include only samples present in the base clinical metadata
met_Patte <- met_Patte[,colnames(met_Patte) %in% rownames(base)]


############### Merge Methylation Data from Multiple Studies ###############

# Merge Alcala and Laddha methylation data by matching on CpG probe identifiers (row.names)
methylation <- merge(met_Alcala, met_Laddha,by.x = "row.names",by.y = "row.names")

# Set row names as CpG probe identifiers and remove redundant first column
rownames(methylation) <- methylation[,1]
methylation <- methylation[,-1]

# Merge the combined methylation matrix (Alcala & Laddha) with Lakis data
methylation <- merge(methylation, met_Lakis, by.x = "row.names",by.y = "row.names")

# Set row names again and remove the redundant first column
rownames(methylation) <- methylation[,1]
methylation <- methylation[,-1]

# Merge the combined matrix (Alcala, Laddha, Lakis) with Patte data
methylation <- merge(methylation, met_Patte, by.x = "row.names",by.y = "row.names")

# Set row names again and remove the redundant first column
rownames(methylation) <- methylation[,1]
methylation <- methylation[,-1]

# Write the combined methylation matrix for all studies to a file
write.table(methylation, paste0(outdir,"/methylation_all.tsv"),
sep = "\t",
row.names = T,
col.names= NA)

############### Filter Probes Used in MOFA Model ###############

# Load the list of probes used in the MOFA model
met <- read.table(paste0(wkdir,
"/preprocess/methylome/methylation_fs.txt"),
    sep = "\t",
    row.names = 1,
    header = T,
    check.names = F
)

# Filter the combined methylation matrix to include only the probes used in the MOFA model
methylation <- methylation[rownames(methylation) %in% rownames(met),]

# Write the filtered methylation matrix for MOFA model to a file
write.table(methylation,
paste0(outdir, "/methylation.fs.txt"),
sep = "\t",
row.names = T,
col.names = NA)

###############################################################################
# Script: Ensemble Clustering for MOFA Model Combination
#
# Description:
#     This script performs ensemble clustering on the MOFA combination model, 
#     using multiple clustering algorithms. It outputs the final subtype assignments 
#     for each sample after clustering, and updates the clinical metadata.
#
# Requirements:
#     - Pre-trained MOFA combination model.
#     - Required libraries: MOFA2, diceR, NMF, ggplot2, foreach, doParallel.
#     - Clustering algorithms: hc, diana, km, pam, gmm, block, cmeans, hdbscan.
#
# Inputs:
#     - MOFA model combination file.
#     - Clinical data with existing validation subtypes.
#     - Config file for directory paths.
#
# Outputs:
#     - Ensemble clustering results.
#     - Updated clinical metadata with the new subtype assignments.
#
# Required environment: `clustering.ensemble.yml`
###############################################################################

############### Load Required Libraries ###############
# Libraries for MOFA, clustering, and visualization

library(MOFA2)
library(diceR)
library(foreach)
library(doParallel)
library(NMF)
library(ggplot2)



############### Load Custom Functions ###############
# Source functions for classification and MOFA-related tasks
source("functions/clustering/classification.functions.R")
source("functions/MOFA/MOFA.functions.R")

############### Load Configuration File ###############
# Load configuration file to get working and data directory paths
config <- read.csv("config/config.tsv",
sep = "\t",
header = T)
# Extract key-value pairs for configuration settings
keys <- config$key
values <- config$value
names(values) <- keys

# Define working directory and data directory from configuration file
wkdir <- values["wkdir"]
datadir <- values["datadir"]

############### Create Output Directories ###############
# Create directories for storing clustering results if they don't already exist

if (dir.exists(paste0(wkdir, "/validation/clustering.combination")) == F) {

    dir.create(paste0(wkdir, "/validation/clustering.combination"))

}

if (dir.exists(paste0(wkdir, "/validation/clustering.combination/clustering_process")) == F) {

    dir.create(paste0(wkdir, "/validation/clustering.combination/clustering_process"))

}

############### Load Clinical Data ###############
# Load clinical metadata that includes existing validation subtypes

base.validation <- read.table(paste0(wkdir, "/validation/base/Base_validation.subtype.tsv"),
sep = '\t',
header = T,
row.names = 1
)

base.discovery <- read.csv(paste0(datadir, "/clinical_data.txt"), sep = "\t", row.names = 1, header = TRUE)

#Create base comparative

base.validation <- base.validation[,"Subtype.validation",drop = FALSE]
base.discovery <- base.discovery[,"Subtype",drop = FALSE]

base.comparation <- merge(base.validation,base.discovery,all = T,by= "row.names")
############### Define Factors to Remove ###############
# Specify factors to exclude from clustering
factor.rm <- c("Factor3", "Factor4","Factor5","Factor6","Factor7","Factor8")


############### Define Clustering Algorithms ###############
# Specify clustering algorithms for ensemble clustering

algorithms <- c('hc', 'diana', 'km', 'pam',
'gmm', 'block', 'cmeans', 'hdbscan')


outdir <- paste0(wkdir, "/validation/clustering.combination/clustering_process")

############### Load MOFA Combination Model ###############
# Load the pre-trained MOFA combination model
mod <- load_model(paste0(wkdir, "/validation/MOFA/model/Mofa.model.combination.hdf5"))

############### Perform Ensemble Clustering ###############
# Run ensemble clustering using the selected algorithms, excluding the specified factors
clustering.ensemble(mod,
factor.rm,
algorithms = algorithms,
outdir = outdir,
scale = F)

############### Load Clustering Results ###############
# Read the input data used for clustering and the combined similarity matrix
clust_input <- readRDS(paste0(outdir, "/clust_input.rds"))

############### Obtain Final Cluster Assignments ###############
# Obtain the final cluster assignments for each sample (k = 3 clusters)
class_scale <- clustering.ensemble.matrix(mod = mod,
factors.rm = factor.rm,
clust_input = clust_input,
k = 3,
outdir = outdir)

# Alternatively, read the pre-saved class matrix for 3 clusters
class_scale <- read.table(paste0(outdir, "/Class.matrix_3.txt"),
sep = "\t",
header = T,
row.names = 1
)

############### Update Clinical Metadata with New Subtypes ###############
# Ensure rownames in class_scale match with base data

rownames(base.comparation) <- base.comparation$Row.names

base.comparation <- base.comparation[row.names(class_scale),]
class_scale <- class_scale[rownames(base.comparation),,drop = F]

# Add the new subtype assignments from class_scale to the base clinical data
base.comparation$Subtype.combination <- class_scale$Subtype

# Write the updated clinical metadata with new subtypes to a file
write.table(base.comparation,
paste0(wkdir, "/validation/base/Base_comparation.subtype.tsv"),
sep = "\t",
row.names = T,
col.names = T)

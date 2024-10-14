###############################################################################
# Script: Ensemble Clustering for MOFA Validation Model
#
# Description:
#     This script runs ensemble clustering on the validation model of MOFA 
#     (Multi-Omics Factor Analysis). It updates the clinical metadata by adding 
#     the new subtypes discovered from the ensemble clustering.
#
# Requirements:
#     - Pre-trained MOFA validation model.
#     - Required libraries: MOFA2, diceR, foreach, doParallel, NMF, ggplot2.
#     - Clustering algorithms: hc, diana, km, pam, gmm, block, cmeans, hdbscan.
#
# Inputs:
#     - MOFA model validation file.
#     - Clinical data.
#     - Config file for directory paths.
#
# Outputs:
#     - Ensemble clustering results.
#     - Updated clinical metadata with the new validation subtypes.
#
# Required environment: `clustering.ensemble.yml`
###############################################################################

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
# Read configuration file to get working and data directory paths
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

############### Create Output Directories ###############
# Create directories for clustering results if they don't already exist

if (dir.exists(paste0(wkdir, "/validation/clustering.validation")) == F) {

    dir.create(paste0(wkdir, "/validation/clustering.validation"))

}

if (dir.exists(paste0(wkdir, "/validation/clustering.validation/clustering_process")) == F) {

    dir.create(paste0(wkdir, "/validation/clustering.validation/clustering_process"))

}

############### Load Clinical Data ###############
# Load clinical metadata that includes subtypes and CIMP status
base <- read.table(paste0(wkdir, "/validation/base/Base_validation_cimp.tsv"),
sep = '\t',
header = T,
row.names = 1
)


############### Define Factors to Exclude ###############
# Only use relevant factors (Factor1 and Factor3) from the discovery dataset for clustering
factor.rm <- c("Factor1", "Factor4","Factor5","Factor6","Factor7","Factor8")


############### Define Clustering Algorithms ###############
# Specify the clustering algorithms for ensemble clustering
algorithms <- c('hc', 'diana', 'km', 'pam',
'gmm', 'block', 'cmeans', 'hdbscan')


outdir <- paste0(wkdir, "/validation/clustering.validation/clustering_process")

############### Load MOFA Validation Model ###############
# Load the pre-trained MOFA validation model
mod <- load_model(paste0(wkdir, "/validation/MOFA/model/Mofa.model.validation.hdf5"))

############### Perform Ensemble Clustering ###############
# Perform ensemble clustering on the validation MOFA model
clustering.ensemble(mod,
factor.rm,
algorithms = algorithms,
outdir = outdir,
scale = F
)

############### Load Clustering Results ###############
# Read the input data used for clustering
clust_input <- readRDS(paste0(outdir, "/clust_input.rds"))

############### Obtain Final Cluster Assignments ###############
# Obtain the final cluster assignments for k=3 clusters
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
base <- base[row.names(class_scale),]
class_scale <- class_scale[rownames(base),,drop = F]

# Add the new subtype assignments from class_scale to the base clinical data
base$Subtype.validation <- class_scale$Subtype

# Write the updated clinical metadata with new validation subtypes to a file
write.table(base,
paste0(wkdir, "/validation/base/Base_validation.subtype.tsv"),
sep = "\t",
row.names = T,
col.names = T)

# This script performs ensemble clustering on a MOFA (Multi-Omics Factor Analysis) model 
# using various clustering algorithms. It integrates multiple methods to determine consensus 
# clustering results and assigns clinical subtypes based on the clustering outcomes.
#
# Requirements:
# - A pre-trained MOFA model that has been saved in HDF5 format.
# - Clinical metadata containing sample information.
# - Configuration for input and output directories specified in "config/config.tsv".
# - Utilizes external clustering functions from "ffunctions/clustering/classification.functions.R".
#
# Inputs:
# - clinical_data.txt: Clinical metadata with sample information.
# - Mofa.model.hdf5: MOFA model file containing integrated omics data.
#
# Outputs:
# - clinical_data.txt: Updated clinical metadata with assigned subtypes based on clustering.
# Required environment: 'clustering.ensemble.yml' (containing prerequiste libraries).

# Load required libraries

library(MOFA2)
library(ComplexHeatmap)
library(diceR)
library(circlize)
library(foreach)
library(doParallel)
library(NMF)
library(ggplot2)

# Load custom functions for clustering 
source("functions/clustering/classification.functions.R")
source("functions/MOFA/MOFA.functions.R")

# Read configuration file for setting up working and data directories
config <- read.csv("config/config.tsv",
sep = "\t",
header = T)

# Create named vectors for easy access to configuration values
keys <- config$key
values <- config$value
names(values) <- keys

# Set working and data directories based on the configuration
wkdir <- values["wkdir"]
datadir <- values["datadir"]

#We create output folder
if (dir.exists(paste0(wkdir, "/clustering")) == F) {

    dir.create(paste0(wkdir, "/clustering"))

}

if (dir.exists(paste0(wkdir, "/clustering/clustering_process")) == F) {

    dir.create(paste0(wkdir, "/clustering/clustering_process"))

}


outdir <- paste0(wkdir, "/clustering/clustering_process")



####################################################
# Load Clinical Data
####################################################

# Load clinical data containing sample information

base <- read.csv(paste0(datadir, "/clinical_data.txt"),
                 sep = "\t",
                 header = T)

####################################################
# Specify Factors and Algorithms for Clustering
####################################################

# List factors to be removed from the analysis
factor.rm <- c("Factor3", "Factor7")


# List of clustering algorithms to be used in the ensemble clustering
algorithms <- c('hc', 'diana', 'km', 'pam',
'gmm', 'block', 'cmeans', 'hdbscan')


####################################################
# Load the MOFA Model
####################################################

# Load the pre-trained MOFA model from the specified file
mod <- load_model(paste0(wkdir, "/MOFA/model/Mofa.model.hdf5"))




####################################################
# Create Similarity Matrix for Each Algorithm
####################################################

# Generate the similarity matrix for clustering based on the specified algorithms

clustering.ensemble(mod,
factor.rm,
algorithms = algorithms,
outdir = outdir,
scale = F)


####################################################
# Read Clustering Results
####################################################

# Read the clustering input data from the saved RDS file

clust_input <- readRDS(paste0(outdir, "/clust_input.rds"))


# Read the combined similarity matrix for each clustering algorithm
comb <- readRDS(paste0(outdir, "/comb.rds"))



####################################################
# Calculate Best Number of Clusters (K)
####################################################

# Calculate the consensus matrix for each clustering result

consensus <- lapply(comb, consensus_matrix)

# Compute cophenetic correlation coefficients for each consensus matrix
coph <- lapply(consensus,cophcor)
# Flatten the list of cophenetic correlations for plotting
coph <- unlist(coph)
# Plot cophenetic correlation values to determine the best K
coph.plot(coph,outdir)

# Determine the best K based on the maximum cophenetic correlation
k <- as.numeric(names(coph[coph == max(coph)]))


####################################################
# Obtain Class of Each Sample Based on Clustering
####################################################

# Generate the class matrix based on the selected best K
class_scale <- clustering.ensemble.matrix(mod = mod,
factors.rm = factor.rm,
clust_input = clust_input,
k = k,
outdir = outdir)

# Load the class scale matrix from the output file for further processing
class_scale <- read.table(paste0(outdir, "/Class.matrix_3.txt"),
sep = "\t",
header = T,
row.names = 1
)


# Adjust the class names to match the clinical data
class_scale <- class_scale[base$CLINICAL_TRIAL_CODE,,drop = F]

# Add the assigned subtypes to the clinical data
base$Subtype <- class_scale$Subtype

#Indicate Neuroendocrine Subtype names
base$Subtype[base$Subtype == 1] <- "NS1"
base$Subtype[base$Subtype == 2] <- "NS2"
base$Subtype[base$Subtype == 3] <- "NS3"

# Save the updated clinical data with the new subtypes back to the file
write.table(base,
paste0(datadir, "/clinical_data.txt"),
sep = "\t",
row.names = F,
col.names = T)

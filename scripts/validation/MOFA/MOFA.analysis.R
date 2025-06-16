####################################################
# Script to analyze MOFA model for multi-omics data
# 
# Requirements:
# - MOFA2 library.
# - Annotation files (for methylation data).
# - MOFA model 
# - Preprocessed clinical and expression/methylation data.
# - Configuration file (to define working/data directories).
# 
# Inputs:
# - MOFA model, preprocessed methylation, expression, and clinical data.
# 
# Outputs:
# - MOFA results, including variance explained, associations with clinical data, 
#   correlation plots, heatmaps, and survival analysis.
####################################################
# Required environment: 'MOFA.yml' (containing prerequiste libraries).

# Load required libraries for MOFA analysis and visualization
library(MOFA2)  # Multi-Omics Factor Analysis
library(ggplot2)  # Visualization
library(tidyverse)  # Data manipulation
library(pheatmap)  # Heatmap plotting
library(ggpubr)  # Publication-ready plots
library(ComplexHeatmap)  # Advanced heatmaps
library(circlize)  # Circular visualizations
library(survival)  # Survival analysis
library(survminer)  # Plotting survival curves
library(fgsea)  # Gene Set Enrichment Analysis
library(reshape2)  # Data reshaping
library(rstatix)  # Statistical tests

# Load custom functions for MOFA analysis
source("functions/MOFA/MOFA.functions.R")  
source("functions/main.functions.R")  

####################################################
# Read configuration file (working and data directories)
####################################################
# Configuration file defines paths for input/output data
config <- read.csv("config/config.tsv", sep = "\t", header = TRUE)

# Extract keys and values from the config file
keys <- config$key
values <- config$value
names(values) <- keys  # Set config keys as names of 'values' for easy access

# Define working and data directories based on config file
wkdir <- values["wkdir"]  # Working directory
datadir <- values["datadir"]  # Data directory

####################################################
# Create output folder for analysis results (if it doesn't exist)
####################################################
if (dir.exists(paste0(wkdir, "/MOFA/analysis_results")) == FALSE) {
    dir.create(paste0(wkdir, "/MOFA/analysis_results"))  # Create output directory for results
}

####################################################
# Prepare methylation annotation files
####################################################
# Load the EPIC hg38 manifest (annotation for methylation data)
# avaiable here https://emea.support.illumina.com/downloads/infinium-methylationepic-v1-0-product-files.html 
epic <- read.csv(paste0(datadir, "/methylome/anotation/infinium-methylationepic-v-1-0-b5-manifest-file.csv"),
                 sep = ",", header = TRUE, skip = 7)

# Prepare the methylation annotation using a custom function
met_anot <- create.met.anot(epic)

####################################################
# Load pre-trained MOFA model
####################################################
# Define output directory for the results
outdir <-paste0(wkdir, "/validation/MOFA")

# Load pre-trained MOFA model from HDF5 file
mod <- load_model(paste0(wkdir, "/validation/MOFA/model/Mofa.model.validation.hdf5"))

# Order clinical data to match the samples in the MOFA model
df <- base[mod@samples_metadata$sample, ]

# Extract factor values (Z) from the MOFA model and append them to the clinical data
Z <- obtain_Z(mod)
df <- cbind(df, Z)

# Merge updated clinical data with the MOFA model's metadata
mod@samples_metadata <- merge(mod@samples_metadata, df, by.x = 'sample', by.y = 'row.names')

# Annotate methylation data with additional information
mod <- anotate.metilation(mod, met_anot)

####################################################
# General Results: Variance Explained and Factor Importance
####################################################
# Create folder for general results (if not already created)
if (!dir.exists(paste0(outdir, "/general_results"))) {
    dir.create(paste0(outdir, "/general_results"))
}

# Write the clinical data with factors to a file for further analysis
write.table(df, paste0(outdir, "/general_results/Clinica_table_factors.txt"), sep = "\t", row.names = TRUE, col.names = NA)

# Plot the variance explained by each factor across different groups
pve <- plot_variance_explained(mod, x = 'group', y = 'factor')
pve <- theme_general(pve)  # Apply a general theme to the plot
ggsave(filename = paste0(outdir, "/general_results/Variance_groups.pdf"), device = "pdf", dpi = 300, width = 29, height = 15)

# Calculate the mean variance explained and save results
calculate.mean.variance(mod, paste0(outdir, "/general_results"))

# Calculate variance explained by each factor for each tissue type and save as table
ve <- calculate_variance_explained(mod)
ve <- Reduce(rbind, ve$r2_per_factor)  # Combine variance results into a single table
ve <- as.data.frame(ve)
ve$tissue <- c(rep("LUNG", 8), rep("PANCREAS", 8), rep("SMALL INTESTINE", 8))
write.table(ve, paste0(outdir, "/general_results/Variance_groups.txt"), sep = "\t", row.names = TRUE, col.names = NA)


####################################################
# Feature Importance Analysis
####################################################


if (dir.exists(paste0(outdir, "/features/")) == F) {

    dir.create(paste0(outdir, "/features/"))

}
#Create Gerneral_gensets.gmt using other gmt files

dir.gensets <- paste0(wkdir, "/genesets/")
gmts <- list.files(dir.gensets)
gmts <- paste0(dir.gensets, gmts)


system(paste0("cat ", paste0( gmts,collapse = " "), " > ", dir.gensets, "General_genesets.gmt"))


# Load pathway information for pathway analysis from a GMT file
pathways <- gmtPathways(paste0(wkdir, "/genesets/General_genesets.gmt"))


# Obtain important features for each omic layer from the MOFA model
# This will analyze which features are most relevant
# to the factors identified in the MOFA model

important.features(mod,
outdir = outdir,
gmt = pathways,
met_anot = met_anot)


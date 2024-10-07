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
# Load clinical data and MOFA model
####################################################
# Read clinical data (contains tumor type, patient characteristics, etc.)
base <- read.csv(paste0(datadir, "/clinical_data.txt"), sep = "\t", row.names = 1, header = TRUE)


# Standardize tumor location information
base$TUMOR_LOCATION <- base$BIOPSY_COMPLETE_LOCATION
base$TUMOR_LOCATION[grep("\\(", base$TUMOR_LOCATION, invert = TRUE)] <- "PRIMARY.TUMOR"
base$TUMOR_LOCATION <- gsub("\\(.*", "", base$TUMOR_LOCATION)

####################################################
# Load pre-trained MOFA model
####################################################
# Define output directory for the results
outdir <- paste0(wkdir, "/MOFA/analysis_results/")

# Load pre-trained MOFA model from HDF5 file
mod <- load_model(paste0(wkdir, "/MOFA/model/Mofa.model.hdf5"))

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
ve$tissue <- c(rep("COLORECTAL", 8), rep("GASTRIC", 8), rep("LUNG", 8), rep("PANCREAS", 8), rep("SMALL INTESTINE", 8))
write.table(ve, paste0(outdir, "/general_results/Variance_groups.txt"), sep = "\t", row.names = TRUE, col.names = NA)

####################################################
# Association with Discrete Variables (e.g., Clinical Data)
####################################################
# Create output directory for discrete variable associations (if not exists)

if (!dir.exists(paste0(outdir, "/discrete.variables"))) {
    dir.create(paste0(outdir, "/discrete.variables"))
}

# Define features to analyze and colors for each feature
features_cm <- c("PRIMARY_TUMOR", "GENDER", "GRADE", "STAGE_IV", "X5HIAA_HIGH", "CIMP.status", "TUMOR_LOCATION")
colors_fcm <- list(
    "PRIMARY_TUMOR" = c("COLORECTAL" = "orange", "LUNG" = "blue", "SMALL INTESTINE" = "yellow", "PANCREAS" = "green", "GASTRIC" = "pink"),
    "GENDER" = c("FEMALE" = "lightpink", "MALE" = "lightskyblue1"),
    "GRADE" = c("G1" = "green", "G2" = "yellow", "G3" = "red"),
    "STAGE_IV" = c("YES" = "red4", "NO" = "forestgreen"),
    "X5HIAA_HIGH" = c("YES" = "steelblue1", "NO" = "brown2"),
    "CIMP.status" = c("Low" = "blue", "Inter" = "green", "High" = "red"),
    "TUMOR_LOCATION" = c("PERITONEUM" = "lightsalmon", "PRIMARY.TUMOR" = "black", "LIVER" = "brown", "LYMPH NODE" = "lightskyblue1", "LUNG" = "blue", "OVARY" = "lightpink", "SKIN" = "gold", "PAROTID GLAND" = "magenta")
)

# Perform association analysis with discrete clinical features

discrete.features(mod = mod, p.th = 0.05, columns = features_cm, outdir = outdir, flag = "clinical", colors = colors_fcm, 
group = "PRIMARY_TUMOR", select = FALSE, levels = c("COLORECTAL", "GASTRIC", "LUNG", "PANCREAS", "SMALL INTESTINE"))




####################################################
# Immune Infiltration Analysis
####################################################

# Load MCP-counter immune infiltration data
mcp <- read.csv(paste0(wkdir, "/immune/MCPcounter.tsv"),
                 sep = "\t",
                 row.names = 1,
                 header = T,
                 check.names = F)


# Merge the immune infiltration data with the MOFA model's sample metadata
# This adds the immune population data to the model
mod@samples_metadata <- merge(mod@samples_metadata, t(mcp),
by.x = "sample", by.y = "row.names")



if (dir.exists(paste0(outdir, "/continous.variables/")) == F) {

    dir.create(paste0(outdir, "/continous.variables/"))

}

# Perform correlation analysis between MOFA factors and the immune population data
# Use `continous.features` to analyze associations between the MCP-counter data 
# and MOFA factors, creating plots for significant associations
continous.features(mod = mod,
p.th = 0.05,
columns = rownames(mcp),
outdir = outdir,
flag = "mcp",
group = "PRIMARY_TUMOR",
levels = c("COLORECTAL", "GASTRIC", "LUNG", "PANCREAS", "SMALL INTESTINE"))


####################################################
# Linear Modeling of Immune Data with MOFA Factors
####################################################

# Perform linear modeling between the MCP-counter immune data and MOFA factors
# This analysis identifies significant associations between immune cell populations
# and the latent factors from the MOFA model
immune_lm(as.data.frame(t(mcp)),
mod,
"mcp", outdir)

####################################################
# Survival Analysis
####################################################

# Create directory for storing survival analysis results if it doesn't exist
if (dir.exists(paste0(outdir, "/survival/")) == F) {

    dir.create(paste0(outdir, "/survival/"))

}



# Perform survival analysis using overall survival (OS) time and event status (EXITUS)
# The `surv.model` function relates MOFA factors to clinical outcomes (overall survival)
surv.model(mod,time = "OS.time", "EXITUS", flag = "overall_survival", outdir)








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


####################################################
# Gene Set Enrichment Analysis (GSEA)
####################################################

# Open and list GMT files for gene set enrichment analysis


# Load expression data for GSEA analysis
expression <- read.table(paste0(wkdir,
"/preprocess/transcriptome/expression.qual.txt"),
sep = '\t',
header = T,
row.names = 1,
check.names = F
)

# Load methylation data for GSEA analysis


methylation <- read.table(paste0(wkdir,
"/preprocess/methylome/methylation.qual.txt"),
sep = '\t',
header = T,
row.names = 1,
check.names = F)


# Perform GSEA on expression data using the loaded gene sets
gsea_factors(expression, mod, gmts = gmts, outdir,
type = "expression", met_anot = NULL, flag = "expression" )

# Perform GSEA on methylation data using the loaded gene sets
gsea_factors(methylation, mod, gmts = gmts, outdir,
type = "methylation", met_anot = met_anot, flag = "methylation")

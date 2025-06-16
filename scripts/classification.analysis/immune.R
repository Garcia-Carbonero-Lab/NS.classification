###################################################
# Script: Immune Cell Analysis using MCP-counter Data
#
# Description:
#      This script processes clinical and immune cell data to analyze the distribution 
#      of immune cell types across different subtypes. It uses the Kruskal-Wallis test 
#      to assess differences between subtypes and generates Wilcoxon plots 
#      to visualize the distribution of immune cells. The data is processed using the 
#      MCP-counter tool and results are saved in output directories.
#
# Requirements:
#      - Predefined environment with clustering functions from `common.functions.R`.
#      - Configuration file (`config.tsv`) with key-value pairs for directory paths.
#      - Clinical data file (`clinical_data.txt`) containing sample features.
#      - MCP-counter results (`MCPcounter.tsv`) containing immune cell proportions.
####################################################
# Required environment: deg.yml

############### Library Setup ###############
# Load necessary libraries for data analysis and visualization
library(ggplot2)
library(ggpubr)
library(rstatix)
library(tidyverse)
library(reshape2)

# Load custom functions for clustering analysis
source("functions/clustering/common.functions.R")


############### Configuration File Loading ###############
# Read the configuration file (config.tsv) to obtain directory paths and other settings
config <- read.csv("config/config.tsv",
sep = "\t",
header = T)

# Extract key-value pairs from the config file for directory setup
keys <- config$key
values <- config$value
names(values) <- keys
datadir <- values["datadir"]

############### Directory Setup ###############
# Set working directory and data directory paths based on the config filewkdir <- values["wkdir"]


# Check if the immune analysis output directory exists, if not create it
if (dir.exists(paste0(wkdir, "/clustering/immune")) == F) {

    dir.create(paste0(wkdir, "/clustering/immune"))

}
# Set the output directory for immune analysis results
outdir <- paste0(wkdir, "/clustering/immune")


############### Data Loading ###############
# Load clinical data from the specified file (clinical_data.txt)
base <- read.csv(paste0(datadir, "/clinical_data.txt"),
                 sep = "\t",
                 header = T,
                 row.names = 1)

# Load MCP-counter results (MCPcounter.tsv) which contain immune cell proportion data

mcp <- read.csv(paste0(wkdir, "/immune/MCPcounter.tsv"),
                 sep = "\t",
                 row.names = 1,
                 header = T,
                 check.names = F)

# Merge clinical data and MCP-counter immune cell data by the sample row names
df.mcp <- merge(base, t(mcp), by = "row.names")



############### Statistical Analysis ###############
# Study the differences in immune cell proportions across subtypes using the Kruskal-Wallis test
# The `kruskal` function evaluates immune cell differences between subtypes

kruskal.mcp <- kruskal(df.mcp,
                       group = "Subtype",         # Group by subtype
                       features = rownames(mcp),  # Immune cell features (MCP-counter results)
                       outdir = outdir,           # Output directory for Kruskal-Wallis results
                       flag = "mcp")              # Flag to append to output file names

############### Plot Generation: Wilcoxon Plots (Complete Data) ###############
# Check if the output directory for complete Wilcoxon plots exists, if not create it

if (dir.exists(paste0(wkdir, "/clustering/immune/complete.plots")) == F) {

    dir.create(paste0(wkdir, "/clustering/immune/complete.plots"))

}

rownames(base) <- base$Row.names



# Prepare Wilcoxon plots for immune cell distributions by subtype
# This uses the `wilcox_plots` function to generate box plots of immune cell data across subtypes

wilcox_plots(data = as.data.frame(t(mcp)),    # Immune cell data (transpose MCP data)
             base = base,                     # Clinical data
             group.x = "Subtype",             # Group by subtype (x-axis)
             group.fill = "Subtype",          # Fill color by subtype
             color.fill = c("gold", "cornflowerblue", "tomato1"),  # Colors for the subtypes
             outdir = paste0(outdir, "/complete.plots"),           # Output directory for plots
             flag = "MCP_Complete")                                # Flag to append to plot file names

############## Plot Generation: Wilcoxon Plots (Primary Tumor) ###############
# Check if the output directory for primary tumor Wilcoxon plots exists, if not create it

if (dir.exists(paste0(wkdir, "/clustering/immune/primary")) == F) {

    dir.create(paste0(wkdir, "/clustering/immune/primary"))

}

# Generate Wilcoxon plots for immune cell distributions by primary tumors
wilcox_plots(data = as.data.frame(t(mcp)),    # Immune cell data (transpose MCP data)
             base = base,                     # Clinical data
             group.x = "PRIMARY_TUMOR",       # Group by primary tumor (x-axis)
             group.fill = "Subtype",          # Fill color by subtype
             color.fill = c("gold", "cornflowerblue", "tomato1"),  # Colors for the subtypes
             outdir = paste0(outdir, "/primary"),                  # Output directory for plots
             flag = "MCP_primary")    
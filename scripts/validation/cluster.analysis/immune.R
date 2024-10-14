###############################################################################
# Script: Immune Cell Analysis using MCPcounter Data
#
# Description:
#     This script analyzes differences in immune cell populations across 
#     subtypes using MCPcounter data, performing statistical tests and generating 
#     visualizations.
#
# Inputs:
#     - Clinical data: Base_validation.subtype.tsv
#     - Immune data: MCPcounter.tsv
#     - Configuration file: config.tsv
#
# Outputs:
#     - Statistical results from Kruskal-Wallis test
#     - Wilcoxon plots for immune cell distributions
#
# Required environment: `deg`
#
###############################################################################

# Load necessary libraries for analysis and visualization
library(ggplot2)
library(ggpubr)
library(rstatix)
library(tidyverse)
library(reshape2)

# Load custom functions for analysis
source("functions/clustering/common.functions.R")

# Read configuration file for directories
config <- read.csv("config/config.tsv",
sep = "\t",
header = T)

keys <- config$key
values <- config$value
names(values) <- keys

# Define working and data directories
wkdir <- values["wkdir"]
datadir <- values["datadir"]

# Create output directory for immune analysis if it doesn't exist
if (dir.exists(paste0(wkdir, "/validation/clustering.analysis/immune")) == F) {

    dir.create(paste0(wkdir, "/validation/clustering.analysis/immune"))

}

outdir <-paste0(wkdir, "/validation/clustering.analysis/immune")

# Load clinical data
base <- read.table(paste0(wkdir, "/validation/base/Base_validation.subtype.tsv"),
                              sep = '\t', header = TRUE, row.names = 1)

# Load MCPcounter data
mcp <- read.csv(paste0(wkdir, "/validation//immune/MCPcounter.tsv"),
                 sep = "\t",
                 row.names = 1,
                 header = T,
                 check.names = F)

# Convert subtype validation column to factor
base$Subtype.validation <- as.factor(as.character(base$Subtype.validation))

# Merge MCPcounter data with clinical data
base <- merge(base, t(mcp), by = "row.names")


# Perform Kruskal-Wallis test to study differences of immune cells
kruskal <- kruskal(base,
"Subtype.validation",
rownames(mcp),
outdir = outdir,
flag = "mcp")

# Create a directory for complete plots if it doesn't exist
if (dir.exists(paste0(wkdir, "/validation/clustering.analysis/immune/complete.plots")) == F) {

    dir.create(paste0(wkdir, "/validation/clustering.analysis/immune/complete.plots"))

}
rownames(base) <- base$Row.names


# Generate Wilcoxon plots for differences of immune cells in each subtype
wilcox_plots(data = as.data.frame(t(mcp)),
base = base,
group.x = "Subtype.validation",
group.fill =  "Subtype.validation",
color.fill = c("gold", "cornflowerblue", "tomato1"),
outdir= paste0(outdir, "/complete.plots"),
flag = "MCP_Complete")

# Create a directory for primary tumor plots if it doesn't exist
if (dir.exists(paste0(wkdir, "/validation/clustering.analysis/immune/primary")) == F) {

    dir.create(paste0(wkdir, "/validation/clustering.analysis/immune/primary"))

}

# Generate Wilcoxon plots for differences of immune cells in each subtype by primary tumors
wilcox_plots(data = as.data.frame(t(mcp)),
base = base,
group.x = "primary.tumor",
group.fill =  "Subtype.validation",
color.fill = c("gold", "cornflowerblue","tomato1"),
outdir= paste0(outdir, "/primary"),
flag = "MCP_primary")

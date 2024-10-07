#######################################################
# Script: Differential Master Regulators Activation Analysis 
#
# Description:
#      This script performs a differential activation analysis
#      on master regulators for different subtypes of tumors. 
#      It compares specific subtypes (NS1, NS2, NS3) against 
#      all other subtypes, generating volcano plots to visualize 
#      the results and filtering significant results based on adjusted 
#      p-values. The results are saved to files for further analysis.
#
# Requirements:
#      - Clinical data file (`clinical_data.txt`)
#      - Master regulator results (`viper.results.txt`)
#      - Functions from the DEG environment and clustering
#######################################################
# Required environment: deg.yml.


############### Load Libraries ###############
# Load required libraries for  analysis, heatmap generation, and plotting
library(limma)
library(ggplot2)
library(rsample)
library(ggpubr)
library(rstatix)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(EnhancedVolcano)

############### Load Functions ###############
# Source custom functions required for analysis
source("functions/main.functions.R")
source("functions/clustering/common.functions.R")
source("functions/clustering/deg.functions.R")

############### Configuration File Loading ###############
# Read the configuration file (config.tsv) for directory paths and settings

config <- read.csv("config/config.tsv",
sep = "\t",
header = T)

# Extract key-value pairs from the config file
keys <- config$key
values <- config$value
names(values) <- keys

# Define working and data directories based on the config file
wkdir <- values["wkdir"]
datadir <- values["datadir"]

############### Directory Setup ###############
# Check if the clustering output directory exists, if not, create it
if (dir.exists(paste0(wkdir, "/clustering/master")) == F) {

    dir.create(paste0(wkdir, "/clustering/master"))

}

# Define the output directory
outdir <- paste0(wkdir, "/clustering/master")


############### Data Loading ###############
# Load clinical data (patient and tumor subtype info)
base <- read.csv(paste0(datadir, "/clinical_data.txt"),
                 sep = "\t",
                 header = T, row.names = 1)

# Load master regulator results (from VIPER analysis)
master <- read.table(paste0(wkdir,
"/master_regulators/viper/viper.results.txt"),
    sep = "\t",
    header = T,
    row.names = 1,
    check.names = F
)

# Clean and adjust the `PRIMARY_TUMOR` factor for downstream analyses
base$PRIMARY_TUMOR <- as.factor(gsub(" ", "_", base$PRIMARY_TUMOR))

############### Prepare Subtypes for  Comparison ###############
# Subtype NS1: All others compared to NS

base$NS1 <- as.character(base$Subtype)
base$NS1[base$NS1 != "NS1"] <- "Other"

# Subtype NS2: All others compared to NS2
base$NS2 <- as.character(base$Subtype)
base$NS2[base$NS2 != "NS2"] <- "Other"

# Subtype NS3: All others compared to NS3
base$NS3 <- as.character(base$Subtype)
base$NS3[base$NS3 != "NS3"] <- "Other"

###############  Analysis for Each Subtype ###############
# Run differential expression analysis for NS1 vs Others
deg.function(master,
base,
group = "NS1",
cov = "PRIMARY_TUMOR",
outdir = paste0(wkdir, "/clustering/master"),
flag = "master_NS1")

# Run differential expression analysis for NS2 vs Others
deg.function(master,
base,
group = "NS2",
cov = "PRIMARY_TUMOR",
outdir = paste0(wkdir, "/clustering/master"),
flag = "master_NS2")


# Run differential expression analysis for NS3 vs Others
deg.function(master,
base,
group = "NS3",
cov = "PRIMARY_TUMOR",
outdir = paste0(wkdir, "/clustering/master"),
flag = "master_NS3")




############### Load and Visualize  Results ###############
# Load the differential expression results for each subtype

NS1 <- read.table(paste0(wkdir,
"/clustering/master/TopTable_master_NS1_GroupNS1-GroupOther.txt"),
sep = "\t",
row.names = 1,
header = T
)

NS2 <- read.table(paste0(wkdir,
"/clustering/master/TopTable_master_NS2_GroupNS2-GroupOther.txt"),
sep = "\t",
row.names = 1,
header = T
)

NS3 <- read.table(paste0(wkdir,
"/clustering/master/TopTable_master_NS3_GroupNS3-GroupOther.txt"),
sep = "\t",
row.names = 1,
header = T
)

############### Volcano Plot Visualization for Each Subtype ###############
# NS1 vs Others: Generate a volcano plot

v.NS1 <- EnhancedVolcano(NS1,
lab = rownames(NS1),
x = "logFC",
y = "adj.P.Val",
title = "NS1 vs Others",
FCcutoff = 0,
pCutoff = 0.1,
subtitle = NULL,
boxedLabels = TRUE,
drawConnectors = TRUE,
labSize = 6.0,
xlim = c(min(NS1[["logFC"]], na.rm = TRUE) - 0.5, max(NS1[["logFC"]], na.rm = TRUE) +
0.5),
ylim = c(0, max(-log10(NS1[["adj.P.Val"]]), na.rm = TRUE) + 1))

ggsave(v.NS1,
filename = paste0(wkdir, "/clustering/master/Volcano_NS1.pdf"),
device = "pdf",
dpi = 500, width = 15,
height = 10
)

# NS2 vs Others: Generate a volcano plot
v.NS2 <- EnhancedVolcano(NS2,
lab = rownames(NS2),
x = "logFC",
y = "adj.P.Val",
title = "SubtypeNS2 vs Others",
FCcutoff = 0,
pCutoff = 0.05,
subtitle = NULL,
boxedLabels = TRUE,
drawConnectors = TRUE,
labSize = 6.0,
xlim = c(min(NS2[["logFC"]], na.rm = TRUE) - 0.5, max(NS2[["logFC"]], na.rm = TRUE) +
0.5),
ylim = c(0, max(-log10(NS2[["adj.P.Val"]]), na.rm = TRUE) + 1))

ggsave(v.NS2,
filename = paste0(wkdir, "/clustering/master/Volcano_NS2.pdf"),
device = "pdf",
dpi = 500, width = 15,
height = 10
)
# NS3 vs Others: Generate a volcano plot
v.NS3 <- EnhancedVolcano(NS3,
lab = rownames(NS3),
x = "logFC",
y = "adj.P.Val",
title = "SubtypeNS3 vs Subtype Others",
FCcutoff = 0,
pCutoff = 0.05,
subtitle = NULL,
boxedLabels = TRUE,
drawConnectors = TRUE,
labSize = 6.0,
xlim = c(min(NS3[["logFC"]], na.rm = TRUE) - 0.5, max(NS3[["logFC"]], na.rm = TRUE) +
0.5),
ylim = c(0, max(-log10(NS3[["adj.P.Val"]]), na.rm = TRUE) + 1))


ggsave(v.NS3,
filename = paste0(wkdir, "/clustering/master/Volcano_NS3.pdf"),
device = "pdf",
dpi = 500, width = 15,
height = 10
)

############### Filter Significant  ###############
# Extract genes with adjusted p-value < 0.05 from each analysis
NS1.genes <- rownames(NS1)[NS1$adj.P.Val < 0.05 ]

NS2.genes <- rownames(NS2)[NS2$adj.P.Val < 0.05 ]

NS3.genes <- rownames(NS3)[NS3$adj.P.Val < 0.05 ]

# Filter tables for significant genes
NS1_filter <- NS1[NS1$adj.P.Val < 0.05 ,]
NS2_filter <- NS2[NS2$adj.P.Val < 0.05 ,]
NS3_filter <- NS3[NS3$adj.P.Val < 0.05 ,]


############### Save Filtered DEG Results ###############
# Save the filtered DEG tables to files

write.table(NS1_filter,
paste0(outdir, "/Table_NS1_filter.txt"),
sep = "\t",
row.names = T,
col.names = NA)

write.table(NS2_filter,
paste0(outdir, "/Table_NS2_filter.txt"),
sep = "\t",
row.names = T,
col.names = NA)

write.table(NS3_filter,
paste0(outdir, "/Table_NS3_filter.txt"),
sep = "\t",
row.names = T,
col.names = NA)


############### Combine Significant Genes and Save ###############
# Combine significant genes across all subtypes
NS.genes <- c(NS1.genes,NS2.genes,NS3.genes)

# Filter the master table for the significant genes
master <- master[NS.genes,]

# Save the filtered master table with significant genes
write.table(master,
paste0(wkdir,
"/clustering/master/Master.NS.txt"),
sep = "\t",
row.names= T,
col.names= NA)
###############################################################################


###############################################################################
# Script: Differential Expression Analysis and Volcano Plot Generation
#
# Description:
#     This script performs differential expression analysis between subtypes 
#     and generates volcano plots to visualize significant DEGs.
#
# Inputs:
#     - Expression data: expression.combination.all.tsv
#     - Validation dataset: Base_validation.subtype.tsv
#     - Configuration file: config.tsv
#
# Outputs:
#     - DEG results for NS1, NS2, NS3
#     - Volcano plots (PDF) for each DEG comparison
#
# Required environment: `deg.yml`
#
###############################################################################

library(limma)
library(ggplot2)
library(rsample)
library(ggpubr)
library(rstatix)
library(tidyverse)
library(EnhancedVolcano)

# Load custom functions
source("functions/clustering/deg.functions.R")

# Read configuration file for working and data directories
config <- read.csv("config/config.tsv",
sep = "\t",
header = T)

keys <- config$key
values <- config$value
names(values) <- keys


wkdir <- values["wkdir"]
datadir <- values["datadir"]

# Create output directory for DEG analysis if it doesn't exist
if (dir.exists(paste0(wkdir, "/validation/clustering.analysis/deg")) == F) {

    dir.create(paste0(wkdir, "/validation/clustering.analysis/deg"))

}

outdir <- paste0(wkdir, "/validation/clustering.analysis/deg")

# Load validation data
base <- read.table(paste0(wkdir, "/validation/base/Base_validation.subtype.tsv"),
                              sep = '\t', header = TRUE, row.names = 1)
# Load gene expression data
exp <- read.table(paste0(wkdir, "/validation/whole.matrix/expression.combination.all.tsv"),
    sep = "\t",
    header = T,
    row.names = 1,
    check.names = F
)
#Define subtypes NS1, NS2, and NS3 for comparison against "Other"
base$NS1 <- base$Subtype.validation
base$NS1[base$NS1 !=  "NS1"] <- "Other"

base$NS2 <- base$Subtype.validation
base$NS2[base$NS2 !=  "NS2"] <- "Other"

base$NS3 <- base$Subtype.validation
base$NS3[base$NS3 !=  "NS3"] <- "Other"

# Perform DEG analysis for NS1, NS2, and NS3 using primary.tumor as covariate
deg.function(exp,
base,
group = "NS1",
cov = "primary.tumor",
outdir = outdir,
flag = "deg_NS1")

deg.function(exp,
base,
group = "NS2",
cov = "primary.tumor",
outdir = outdir,
flag = "deg_NS2")

deg.function(exp,
base,
group = "NS3",
cov = "primary.tumor",
outdir = outdir,
flag = "deg_NS3")

# Read DEG results
NS1 <- read.table(paste0(outdir,"/TopTable_deg_NS1_GroupNS1-GroupOther.txt"),
sep = "\t",
row.names = 1,
header = T
)

NS2 <- read.table(paste0(outdir,"/TopTable_deg_NS2_GroupNS2-GroupOther.txt"),
sep = "\t",
row.names = 1,
header = T
)

NS3 <- read.table(paste0(outdir, "/TopTable_deg_NS3_GroupNS3-GroupOther.txt"),
sep = "\t",
row.names = 1,
header = T
)

# Generate and save volcano plot for NS1 vs Others
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
filename = paste0(outdir, "/Volcano_NS1.pdf"),
device = "pdf",
dpi = 500, width = 15,
height = 10
)

# Generate and save volcano plot for NS2 vs Others
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
filename = paste0(outdir, "/Volcano_NS2.pdf"),
device = "pdf",
dpi = 500, width = 15,
height = 10
)
# Generate and save volcano plot for NS3 vs Others
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
filename = paste0(outdir, "/Volcano_NS3.pdf"),
device = "pdf",
dpi = 500, width = 15,
height = 10
)

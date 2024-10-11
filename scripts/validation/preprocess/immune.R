###########################################################
# Script: Immune Cell Infiltration Estimation using MCPcounter
#
# Description:
#     This script estimates the abundance of different immune and stromal cell populations 
#     using the MCPcounter algorithm. It processes the gene expression matrix from 
#     multiple studies, calculates immune cell infiltration scores, and outputs the results.
#
# Requirements:
#     - R environment with required packages: `MCPcounter`.
#     - Pre-processed expression data from multiple studies.
#
# Inputs:
#     - Combined expression matrix (`expression.combination.all.tsv`).
#
# Outputs:
#     - MCPcounter immune infiltration scores for all samples (`MCPcounter.tsv`).
#
# Required environment: `immune` (for processing immune cell infiltration with MCPcounter).
###########################################################

############### Load Required Libraries ###############
library(MCPcounter)


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




############### Load Pre-processed Expression Matrix ###############
# Load the combined expression matrix that includes data from multiple studies

expression <- read.table(paste0(wkdir, "/validation/whole.matrix/expression.combination.all.tsv"),
    sep = "\t",
    header = T,
    row.names = 1,
    check.names = F
)



############### Directory Creation for Immune Estimation Results ###############
# Create directory for storing immune infiltration results if it doesn't exist
if (dir.exists(paste0(wkdir, "/validation/immune")) == F) {

    dir.create(paste0(wkdir, "/validation/immune"))

}


############### Calculate MCPcounter Scores ###############
# Calculate MCPcounter immune and stromal cell population scores using the expression data
# Specify the feature type as "HUGO_symbols" (gene symbols) for the expression matrix

mcp <- MCPcounter.estimate(expression, featuresType = "HUGO_symbols")

# Ensure row names have no spaces (replace spaces with dots if necessary)
rownames(mcp) <- gsub(" ",".",rownames(mcp))


############### Save MCPcounter Results ###############
# Save the MCPcounter score matrix to the specified directory as a TSV file
write.table(mcp,
paste0(wkdir, "/validation/immune/MCPcounter.tsv"),
sep = "\t",
row.names = T,
col.names = NA)


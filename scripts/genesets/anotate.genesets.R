####################################################
# This script prepares gene set collections (such as Hallmarks, KEGG, Reactome, etc.) for further analysis.
# It reads gene set files in .gmt format and an expression matrix, ensures that the gene symbols in the gene sets
# match the symbols in the expression matrix, and writes the updated gene sets back into .gmt files.
#
# Requirements:
# - Gene sets in .gmt format (such as Hallmarks, KEGG, Reactome, etc.).
# - Expression matrix with gene symbols as row names.
# - Configuration file "config/config.tsv" to define working and data directories.
#
# Inputs:
# - Expression matrix in "preprocess/transcriptome/RMA_Expression_data_genes.tsv".
# - Gene set files in .gmt format (e.g., Hallmarks.gmt, KEGG.gmt).
#
# Outputs:
# - Updated .gmt files with gene symbols that match those in the expression matrix.
####################################################

# Required environment: 'geneset.yml' (containing prerequiste libraries).

# Load required libraries
library(limma)
library(GSEABase)

# Load configuration file to set working and data directories
config <- read.csv("config/config.tsv",
sep = "\t",
header = T)

# Extract keys and values from the configuration file
keys <- config$key
values <- config$value
names(values) <- keys

# Define working and data directories from the config file
wkdir <- values["wkdir"]
datadir <- values["datadir"]

# Source external functions for handling gene sets
source("functions/genesets/genesets.functions.R")


####################################################
# Create output directory if it doesn't exist
####################################################

if (dir.exists(paste0(wkdir, "/genesets")) == F) {

    dir.create(paste0(wkdir, "/genesets"))

}

####################################################
# Read expression matrix with gene symbols
####################################################

expression <- read.table(paste0(wkdir,
"/preprocess/transcriptome/RMA_Expression_data_genes.tsv"),
sep = '\t',
header = T,
row.names = 1,
check.names = F
)


####################################################
# Load gene sets from .gmt files
####################################################
hallmarks <- getGmt(paste0(wkdir, "/genesets/hallmarks.gmt"))
cells <- getGmt(paste0(wkdir, "/genesets/cells.gmt"))
gobp <- getGmt(paste0(wkdir, "/genesets/gobp.gmt"))
biocarta <- getGmt(paste0(wkdir, "/genesets/biocarta.gmt"))
kegg <- getGmt(paste0(wkdir, "/genesets/kegg.gmt"))
pid <- getGmt(paste0(wkdir, "/genesets/pid.gmt"))
reactome <- getGmt(paste0(wkdir, "/genesets/reactome.gmt"))



####################################################
# Annotate gene sets to match the gene symbols in the expression matrix
####################################################
# Match gene symbols in the gene sets to the symbols in the expression matrix

hallmarks <- anotate.gmt(hallmarks, expression)
cells <- anotate.gmt(cells, expression)
gobp <- anotate.gmt(gobp, expression)
biocarta <- anotate.gmt(biocarta, expression)
kegg <- anotate.gmt(kegg, expression)
pid <- anotate.gmt(pid, expression)
reactome <- anotate.gmt(reactome, expression)



####################################################
# Save the updated gene sets as .gmt files
####################################################
# Write the updated gene sets back to .gmt files for downstream analysis
write.sampleset(hallmarks,paste0(wkdir, "/genesets/hallmarks.gmt"))
write.sampleset(cells,paste0(wkdir, "/genesets/cells.gmt"))
write.sampleset(gobp,paste0(wkdir, "/genesets/gobp.gmt"))
write.sampleset(biocarta,paste0(wkdir, "/genesets/biocarta.gmt"))
write.sampleset(kegg,paste0(wkdir, "/genesets/kegg.gmt"))
write.sampleset(pid,paste0(wkdir, "/genesets/pid.gmt"))
write.sampleset(reactome,paste0(wkdir, "/genesets/reactome.gmt"))


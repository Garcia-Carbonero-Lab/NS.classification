# Preprocessing script for transcriptomic and methylomic data correction
#
# Requirements:
# - R environment with the necessary packages installed (umap, ggplot2, limma).
# - The script depends on external functions for preprocessing and PCA correlation.
# - A configuration file ("config.tsv") for setting directory paths and other parameters.
# - Normalized expression and methylation data files.
#
# The script performs PCA analysis and batch effect correction on the data before and 
# after quality correction. It generates quality control plots to visualize PCA results.

# Required environment: 'preprocess.yml' (containing prerequiste libraries).

# Load necessary libraries
library(umap)
library(ggplot2)
library(limma)

# Load custom functions from the preprocess and main functions scripts
source("functions/preprocess/preprocess.R")
source("functions/main.functions.R")


# Load the configuration file (contains paths and parameters)
config <- read.csv("config/config.tsv",
    sep = "\t",
    header = T
)

# Extract key-value pairs from the config file
keys <- config$key
values <- config$value
names(values) <- keys

# Define working and data directories from the config file
wkdir <- values["wkdir"]
datadir <- values["datadir"]

# Load normalized expression data (transcriptomic data)
expression <- read.table(paste0(wkdir,
"/preprocess/transcriptome/RMA_Expression_data_genes.tsv"),
sep = '\t',
header = T,
row.names = 1,
check.names = F
)

# Load normalized methylation data (M-values matrix)
methylation <- read.table(paste0(wkdir,
"/preprocess/methylome/myNormM.txt"),
sep = '\t',
header = T,
row.names = 1,
check.names = F

)

# Load clinical data (e.g., metadata, clinical information about the samples)
clinical.data <- read.table(paste0(datadir, "/clinical_data.txt"),
sep = '\t',
header = T,
row.names = 1
)

#######################################################
# Create necessary directories for quality plots if they do not exist

if (dir.exists(paste0(wkdir, "/preprocess")) == F) {

    dir.create(paste0(wkdir, "/preprocess"))

}

if (dir.exists(paste0(wkdir, "/preprocess/quality.plots")) == F) {

    dir.create(paste0(wkdir, "/preprocess/quality.plots"))
    dir.create(paste0(wkdir, "/preprocess/quality.plots/transcriptome"))
    dir.create(paste0(wkdir, "/preprocess/quality.plots/methylome"))

}



#######################################################
# PCA Analysis and Quality Correction - Before Correction
#######################################################

# Perform PCA on expression data and generate PCA correlation plots before quality correction
pca_exp <- pca_corr(data = expression,
base = clinical.data,
columns = "pos.vs.neg.auc",
outpath = paste0(wkdir, "/preprocess/quality.plots/transcriptome"),
flag = "before.expression"
)

# Perform PCA on methylation data and generate PCA correlation plots before quality correction
pca_met <- pca_corr(data = methylation,
base = clinical.data,
columns =  "qual.met",
outpath = paste0(wkdir, "/preprocess/quality.plots/methylome"),
flag = "before.methylation",
anotate1 = c(650, 450), # Annotation coordinates for the plot
anotate2 = c(650, 420)  # Annotation coordinates for the plot
)


#########################################################
# Batch Effect Correction - Eliminate the first principal component (PC1)
#########################################################

# Generate PCA for expression data and correct for batch effects using PC1
pcaexp <- generate_pca(expression)

pcaexp <- pcaexp[colnames(expression),,drop = F]

# Remove batch effect from expression data using PC1
expression_qual <- removeBatchEffect(x = expression,
covariates = pcaexp[, "PC1"])

# Save the batch-corrected expression data
write.table(expression_qual,
paste0(wkdir,
"/preprocess/transcriptome/expression.qual.txt"),
sep = "\t",
row.names = T,
col.names = NA
)


# Generate PCA for methylation data and correct for batch effects using PC1
pcamet <- generate_pca(methylation)

pcamet <- pcamet[colnames(methylation),,drop = F]

# Remove batch effect from methylation data using PC1
methylation_qual <- removeBatchEffect(x = methylation,
covariates = pcamet[, "PC1"])

# Save the batch-corrected methylation data
write.table(methylation_qual,
paste0(wkdir,
"/preprocess/methylome/methylation.qual.txt"),
sep = "\t",
row.names = T,
col.names = NA
)



#########################################################
# PCA Analysis and Quality Correction - After Correction
#########################################################

# Perform PCA on quality-corrected expression data and generate PCA correlation plots

pca_exp_qual <- pca_corr(data = expression_qual,
base = clinical.data,
columns = "pos.vs.neg.auc",
outpath = paste0(wkdir, "/preprocess/quality.plots/transcriptome"),
flag = "after.expression",
anotate1 = c(50, 90), # Annotation coordinates for the plot
anotate2 = c(50, 80)  # Annotation coordinates for the plot

)

# Perform PCA on quality-corrected methylation data and generate PCA correlation plots

pca_met_qual <- pca_corr(data = methylation_qual,
base = clinical.data,
columns = "qual.met",
outpath = paste0(wkdir, "/preprocess/quality.plots/methylome"),
flag = "after.methylation",
anotate1 = c(350, 500), # Annotation coordinates for the plot
anotate2 = c(350, 470)  # Annotation coordinates for the plot
)

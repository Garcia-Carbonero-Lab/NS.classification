####################################################
# This script sets up and runs a MOFA model using multi-omics data (e.g., expression and methylation).
#
# Requirements:
# - MOFA2 package.
# - Expression and methylation feature files.
# - Clinical data (grouping variable).
# - Configuration file (config/config.tsv) to define working and data directories.
#
# Inputs:
# - Preprocessed expression and methylation data.
# - Clinical data for grouping.
#
# Outputs:
# - MOFA model saved as an HDF5 file.
####################################################
# Required environment: 'MOFA.yml' (containing prerequiste libraries).

# Load MOFA2 library for multi-omics factor analysis
library(MOFA2)


####################################################
# Load configuration file and set working/data directories
####################################################

config <- read.csv("config/config.tsv",
    sep = "\t",
    header = T
)

keys <- config$key
values <- config$value
names(values) <- keys
# Define working and data directories from the config file
wkdir <- values["wkdir"]
datadir <- values["datadir"]


####################################################
# Create output directory for MOFA results if it doesn't exist
####################################################

if (dir.exists(paste0(wkdir, "/MOFA")) == F) {
    dir.create(paste0(wkdir, "/MOFA"))
}


####################################################
# Source MOFA-specific functions from external script (if applicable)
####################################################

source("functions/MOFA/MOFA.functions.R")


####################################################
# Load multi-omics data (Expression and Methylation features)
####################################################
# Load expression features from preprocessed file

exp <- read.table(paste0(wkdir,
"/preprocess/transcriptome/expression_fs.txt"),
    sep = "\t",
    header = T,
    row.names = 1,
    check.names = F
)

# Load methylation features from preprocessed file
met <- read.table(paste0(wkdir,
"/preprocess/methylome/methylation_fs.txt"),
    sep = "\t",
    row.names = 1,
    header = T,
    check.names = F
)


####################################################
# Load clinical data for grouping and ordering
####################################################
# Load clinical data, which contains sample information for grouping (e.g., PRIMARY_TUMOR)

base <- read.csv(paste0(datadir, "/clinical_data.txt"),
    sep = "\t",
    row.names = 1,
    header = T
)


####################################################
# Align features with clinical data (same sample order)
####################################################
# Ensure the expression and methylation data columns match the clinical data sample rows

exp <- exp[,row.names(base)]
met <- met[,row.names(base)]




####################################################
# Create list for MOFA input (multi-omics data)
####################################################

data.list <- list(
    "Expression" = as.matrix(exp),
    "Methylation" = as.matrix(met)
)


####################################################
# Create output directory for MOFA model if it doesn't exist
####################################################

if (dir.exists(paste0(wkdir, "/MOFA/model")) == F) {
    dir.create(paste0(wkdir, "/MOFA/model"))
}


####################################################
# Create and run MOFA model with the specified parameters
####################################################
# Run the MOFA model on the multi-omics data with 8 factors

Mofa.model(
    input = data.list, file = paste0(
        wkdir,
        "/MOFA/model/Mofa.model.hdf5"
    ),
    groups = base$PRIMARY_TUMOR,
    num_factors = 8,
    convergence_mode = "slow"
)

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

if (dir.exists(paste0(wkdir, "/validation/MOFA")) == F) {
    dir.create(paste0(wkdir, "/validation/MOFA"))
}

outdir <- paste0(wkdir, "/validation/MOFA")

####################################################
# Source MOFA-specific functions from external script (if applicable)
####################################################

source("functions/MOFA/MOFA.functions.R")


####################################################
# Load multi-omics data (Expression and Methylation features)
####################################################
# Load expression features from preprocessed file

exp.val <- read.table(paste0(wkdir, "/validation/whole.matrix/expresssion.fs.txt"),
    sep = "\t",
    header = T,
    row.names = 1,
    check.names = F
)

# Load methylation features from preprocessed file
met.val <- read.table(paste0(wkdir, "/validation/whole.matrix/methylation.fs.txt"),
    sep = "\t",
    row.names = 1,
    header = T,
    check.names = F
)


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


exp.comb <- merge(exp,exp.val,by = "row.names")
rownames(exp.comb) <- exp.comb[,1]
exp.comb <- exp.comb[,-1]

met.comb <- merge(met,met.val,by = "row.names")
rownames(met.comb) <- met.comb[,1]
met.comb <- met.comb[,-1]


####################################################
# Load clinical data for grouping and ordering
####################################################
# Load clinical data, which contains sample information

base.val <- read.table(paste0(wkdir, "/validation/base/Base_validation.tsv"),
sep = "\t",
row.names = 1,
header = T)

base <- read.csv(paste0(datadir, "/clinical_data.txt"),
    sep = "\t",
    row.names = 1,
    header = T
)

base.val <- base.val[,"primary.tumor",drop = FALSE]
base <- base[,"PRIMARY_TUMOR",drop = FALSE]

colnames(base) <- "primary.tumor"

base.val[,"primary.tumor"] <- as.factor(base.val[,"primary.tumor"])

levels(base.val[,"primary.tumor"]) <- paste0(levels(base.val[,"primary.tumor"]),".validation")

base.comb <- rbind(base,base.val)
####################################################
# Align features with clinical data (same sample order)
####################################################
# Ensure the expression and methylation data columns match the clinical data sample rows

exp.val <- exp.val[,row.names(base.val)]
met.val <- met.val[,row.names(base.val)]

exp.comb <- exp.comb[,row.names(base.comb)]
met.comb <- met.comb[,row.names(base.comb)]

####################################################
# Create list for MOFA input (multi-omics data)
####################################################

data.list.val <- list(
    "Expression" = as.matrix(exp.val),
    "Methylation" = as.matrix(met.val)
)

data.list.comb <- list(
    "Expression" = as.matrix(exp.comb),
    "Methylation" = as.matrix(met.comb)
)
####################################################
# Create output directory for MOFA model if it doesn't exist
####################################################

if (dir.exists(paste0(wkdir, "/validation/MOFA/model")) == F) {
    dir.create(paste0(wkdir, "/validation/MOFA/model"))
}


####################################################
# Create and run MOFA model with the specified parameters
####################################################
# Run the MOFA model on the multi-omics data with 8 factors

Mofa.model(
    input = data.list.val, file = paste0(
        wkdir,
        "/validation/MOFA/model/Mofa.model.validation.hdf5"
    ),
    groups = base.val$primary.tumor,
    num_factors = 8,
    convergence_mode = "slow"
)


Mofa.model(
    input = data.list.comb, file = paste0(
        wkdir,
        "/validation/MOFA/model/Mofa.model.combination.hdf5"
    ),
    groups = base.comb$primary.tumor,
    num_factors = 8,
    convergence_mode = "slow"
)

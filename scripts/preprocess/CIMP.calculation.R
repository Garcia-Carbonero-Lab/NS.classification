# Script for processing methylation data and calculating CIMP status
#
# Requirements:
# - R environment with required functions and data.
# - This script depends on a preprocessing environment and external functions
#   loaded from "methylation.process.R".
# - The configuration file ("config.tsv") is used to specify working directories.
# - Methylation data (quality-corrected) and clinical data files are required.
#
# The script calculates the CpG island methylator phenotype (CIMP) status for samples,
# merges it with the clinical data, and writes the updated clinical data back to disk.


# Load configuration file that contains paths and other parameters
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


# Load custom functions from the methylation processing file

source("functions/preprocess/methylation.R")


# Load clinical data (e.g., metadata, clinical information about the samples)
clinical.data <- read.table(paste0(datadir, "/clinical_data.txt"),
    sep = "\t",
    header = T
)

# Load quality-corrected methylation data (M-values matrix)
methylation <- read.table(paste0(wkdir,
"/preprocess/methylome/methylation.qual.txt"),
sep = '\t',
header = T,
row.names = 1,
check.names = F

)

# Convert methylation data to numeric format (if not already in numeric form)
methylation <- apply(methylation, c(1,2), as.numeric)


#######################################################
# Calculate CIMP (CpG Island Methylator Phenotype) Status
#######################################################

# Calculate CIMP status for each sample using the addCIMPstatus function
CIMP <- addCIMPstatus(methylation)



# Merge the calculated CIMP status with the clinical data
# Merge by "CLINICAL_TRIAL_CODE" in the clinical data and row names (sample IDs) in CIMP

clinical.data <- merge(clinical.data, CIMP,
    by.x = "CLINICAL_TRIAL_CODE",
    by.y = "row.names"

)

#######################################################
# Save the updated clinical data with the new CIMP status
#######################################################

# Write the updated clinical data to the original clinical data file

write.table(clinical.data, paste0(datadir,"/clinical_data.txt"),
    sep = "\t", row.names = F, col.names = T
)

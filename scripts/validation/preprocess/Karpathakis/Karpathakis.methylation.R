####################################################
# Karpathakis Methylation Data Processing Script
#
# This script processes methylation data from the Karpathakis SINET dataset.
# It performs the following tasks:
# 1. Loads configuration settings to define working and data directories.
# 2. Creates output directories if they do not already exist.
# 3. Reads the raw methylation data from an external file.
# 4. Transforms beta values into M values (logarithmic scale) for further analysis.
#
# Requirements:
# - GEOquery and limma libraries must be installed.
# - A configuration file "config/config.tsv" with paths for working and data directories.
# - Methylation data should be available in "matrix.txt" format.
#
# Inputs:
# - matrix.txt: Methylation data in beta values format from the Karpathakis SINET dataset.
#
# Outputs:
# - Methylation_beta_Karpathakis.tsv: Beta values output.
# - Methylation_M_Karpathakis.tsv: Log-transformed M values output.
#
# Required environment: 'geo.yml' (containing prerequisite libraries).
####################################################



####################################################
# Load configuration and set directories
####################################################
# Open the configuration file to set paths for working and data directories
config <- read.csv("config/config.tsv",
sep = "\t",
header = T)
# Extract key-value pairs for configuration settings
keys <- config$key
values <- config$value
names(values) <- keys

# Set working directory and data directory based on configuration values
wkdir <- values["wkdir"]
datadir <- values["datadir"]

####################################################
# Create Directories for Results
####################################################
# Check if the validation directory exists, create it if it doesn't
if (dir.exists(paste0(wkdir, "/validation")) == F) {

    dir.create(paste0(wkdir, "/validation"))

}

# Check if the Karpathakis validation directory exists, create it if it doesn't
if (dir.exists(paste0(wkdir, "/validation/Karpathakis")) == F) {

    dir.create(paste0(wkdir, "/validation/Karpathakis"))

}

# Define output directory for storing Karpathakis methylation data
outdir <- paste0(wkdir, "/validation/Karpathakis")


####################################################
# Read Methylation Data
####################################################
# Read the methylation data file from the specified directory
# The data is in tab-separated format and we skip the first 90 lines to exclude metadata
methylation <- read.csv(paste0(datadir, "/external_data/Karpathakis_SINET/matrix.txt"),
sep = "\t",
row.names = 1,
header = T,
skip = 90)

# Write the beta values to an output file
write.table(methylation,
paste0(outdir, "/Methylation_beta_Karpathakis.tsv"),
sep = "\t",
row.names = T,
col.names = NA)


####################################################
# Transform Beta Values to M Values
####################################################
# Ensure beta values fall within a valid range to prevent infinite results during log transformation

methylation[methylation <= 0.001] <- 0.001
methylation[methylation >= 0.999] <- 0.999


# Convert the data to numeric format and apply the log transformation
# The formula M = log2(beta / (1 - beta)) is used to convert beta values into M values

methylation <- apply(methylation,c(1,2),as.numeric)
methylation <- log((methylation / (1 - methylation)), 2)

# Write the transformed M values to an output file
write.table(methylation,
paste0(outdir, "/Methylation_M_Karpathakis.tsv"),
sep = "\t",
row.names = T,
col.names = NA)


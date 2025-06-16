####################################################
# Script: Clinical Data Processing
#
# Description:
#      This script processes clinical data to evaluate the distribution of clinical features across
#      different subtypes. It generates tables with feature distributions and saves
#      the results in output files. The script leverages the `clinical.table` function to perform Fisher
#      tests and calculate percentages for each feature across the specified groups.
#
# Requirements:
#      - Predefined environment with clinical functions from `clinical.functions.R`.
#      - Configuration file (`config.tsv`) with key-value pairs for directory paths.
#      - Clinical data file (`clinical_data.txt`) containing sample features.
####################################################
# Required environment: any with R base installed.

# Load required clinical functions
source("functions/clustering/clinical.functions.R")

# Read configuration file (config.tsv) to obtain necessary directory paths
config <- read.csv("config/config.tsv",
                   sep = "\t",
                   header = TRUE)

# Assign configuration keys and values
keys <- config$key
values <- config$value
names(values) <- keys

############### Directory Setup ###############
# Working directory and data directory paths from config
wkdir <- values["wkdir"]
datadir <- values["datadir"]

# Check if clinical output directory exists, if not create it
if (!dir.exists(paste0(wkdir, "/clustering/clinical"))) {
    dir.create(paste0(wkdir, "/clustering/clinical"))
}

# Set output directory
outdir <- paste0(wkdir, "/clustering/clinical")

############### Data Loading ###############
# Read the clinical data from the specified file
base <- read.csv(paste0(datadir, "/clinical_data.txt"),
                 sep = "\t",
                 row.names = 1,
                 header = TRUE)

# Define the clinical features to analyze
features <- c("AGE.median",
              "BIOPSY_SITES", 
              "BIOPSY_METASTASIS_LOCATION",
              "GENDER", 
              "PRIMARY_TUMOR",
              "GRADE",
              "BIOPSY_LOCATION", 
              "X5HIAA_HIGH",
              "CIMP.status", 
              "STAGE_IV")

############### Clinical Feature Analysis ###############
# Create a table showing the distribution of clinical features across subtypes
table <- clinical.table(data = base,
                        features = features,
                        group = "Subtype")

# Create a similar table but include NA values in the analysis
table_na <- clinical.table(data = base,
                           features = features,
                           group = "Subtype",
                           na = TRUE)

############### Output Results ###############
# Save the generated tables into output files

# Save table showing clinical features distribution across subtypes
write.table(table,
            paste0(outdir, "/Clinical.txt"),
            sep = "\t",
            row.names = TRUE,
            col.names = NA,
            quote = FALSE)

# Save table with NA values, showing clinical features distribution across subtypes
write.table(table_na,
            paste0(outdir, "/Clinical_subtypes_na.txt"),
            sep = "\t",
            row.names = TRUE,
            col.names = NA,
            quote = FALSE)

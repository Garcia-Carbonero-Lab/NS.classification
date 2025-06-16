####################################################
# Lakis Clinical Data Processing Script
#
# This script processes clinical data from the Lakis PANET dataset.
# It merges donor, specimen, and sample clinical information, extracts
# relevant columns, and reformats the data for downstream analysis.
#
# Requirements:
# - A configuration file "config/config.tsv" with paths for working and data directories.
# - Clinical data files from the Lakis PANET dataset should be available.
#
# Inputs:
# - clinical_data_donor.tsv: Clinical data for donors.
# - specimen.PAEN-AU.tsv: Specimen data.
# - clinical.data.tsv: Clinical sample data.
#
# Outputs:
# - base_Lakis_prev.tsv: A merged and reformatted clinical dataset ready for further analysis.
#
# Required environment: any with R installed
####################################################

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

# Check if the Lakis validation directory exists, create it if it doesn't
if (dir.exists(paste0(wkdir, "/validation/Lakis")) == F) {

    dir.create(paste0(wkdir, "/validation/Lakis"))

}

# Define output directory for storing Lakis clinical data
outdir <- paste0(wkdir, "/validation/Lakis")


####################################################
# Read Clinical Data from Lakis PANET Dataset
####################################################
# Load clinical data for donors, specimens, and samples from external files

base_icgc_donor <- read.table(paste0(datadir, "/external_data/Lakis_PANET/clinical_data_donor.tsv"),
sep = "\t",
header= T)


base_icgc_specimen <- read.table(paste0(datadir, "/external_data/Lakis_PANET/specimen.PAEN-AU.tsv"),
sep = "\t",
header= T)

base_icgc_sample<- read.table(paste0(datadir, "/external_data/Lakis_PANET/clinical.data.tsv"),
sep = "\t",
header= T)

####################################################
# Merge Clinical Data
####################################################
# Merge the donor and specimen data on the 'icgc_donor_id' field
base <- merge(base_icgc_donor,base_icgc_specimen,by.x = "icgc_donor_id", by.y = "icgc_donor_id")

# Merge the combined data with the sample data on the 'icgc_specimen_id' field
base <- merge(base,base_icgc_sample,by.x = "icgc_specimen_id", by.y = "icgc_specimen_id")



####################################################
# Select and Rename Relevant Columns
####################################################
# Select relevant columns from the merged dataset for downstream analysis
columns <- c("icgc_sample_id","icgc_donor_id.x", "donor_sex", 
"donor_vital_status", "donor_age_at_diagnosis", "donor_tumour_stage_at_diagnosis",
"donor_survival_time", "percentage_cellularity.x", "disease_status_last_followup",
"specimen_type", "tumour_histological_type", "tumour_grade", "tumour_stage")

base <- base[,columns]

# Rename the selected columns for consistency and clarity
colnames(base) <- c("id_sample", "id_patient", "Gender", "Status", "Age", "Stage"
, "Overall.survival", "Celularity", "dissease", "primary.tumor",
"Differenciation","Grade", "Stage")


####################################################
# Write Processed Data to Output File
####################################################
# Write the cleaned and reformatted clinical data to an output file for further analysis
write.table(base,
paste0(outdir, "/base_Lakis_prev.tsv"),
sep = "\t",
row.names = F,
col.names = T)


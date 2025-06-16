############################################################
# Script: Methylation Data Processing for Alcala Dataset
#
# Description:
#     This script processes Illumina EPIC methylation array data from the Alcala 
#     dataset. It generates a sample sheet, organizes IDAT files, performs normalization, 
#     and writes the final methylation matrix in M-value format for downstream analysis.
#
# Requirements:
#     - Libraries for handling methylation data: minfi, ENmix, wateRmelon, MASS, etc.
#     - Custom functions for methylation processing.
#     - Methylation data files in IDAT format.
#     - Configuration file (config.tsv) to define working and data directories.
############################################################

# Required environment: 'met.arrays.yml' (containing prerequiste libraries).

#this iscript uses environment met.arrays.yml
library(minfi)
library(ENmix)
library(wateRmelon)
library(MASS)
library(broom)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)


# Load custom functions, if applicable
source("functions/validation/validation.functions.R")

############### Load Configuration ###############
# Read the configuration file to define paths
config <- read.csv("config/config.tsv",
sep = "\t",
header = T)

# Extract configuration keys and values
keys <- config$key
values <- config$value
names(values) <- keys

# Define working and data directories from configuration file
wkdir <- values["wkdir"]
datadir <- values["datadir"]


############### Directory Setup ###############
# Create directories for validation results if they don't already exist
if (dir.exists(paste0(wkdir, "/validation")) == F) {

    dir.create(paste0(wkdir, "/validation"))

}

if (dir.exists(paste0(wkdir, "/validation/Alcala")) == F) {

    dir.create(paste0(wkdir, "/validation/Alcala"))

}

# Set output directory for processed Alcala data
outdir <- paste0(wkdir, "/validation/Alcala")



############### Define IDAT File Path ###############
# Specify the path to the IDAT files (Illumina array raw data)
path <- paste0(datadir, "/external_data/Alcala_LUNET/IDAT")

############### Generate Sample Sheet ###############
# Function to generate a sample sheet from IDAT file names
generate.sample.sheet.alcala <- function(path) {

# List all IDAT files in the specified directory
files <- list.files(path)


# Filter to keep only files with ".idat" extension    
files <- files[grep("\\.idat",files)]


# Split filenames to extract relevant parts (e.g., Sentrix ID)
flist <- strsplit(files, "_")
names <- gsub(".idat","",files)
names <- gsub("_Grn","",names)
names <- gsub("_Red","",names)

# Extract Sentrix IDs and position
s_entrix_id <- unlist(lapply(flist,function(x){x[[1]]}))
s_entrix <- unlist(lapply(flist,function(x){x[[2]]}))
s_entrix_id <- unlist(s_entrix_id)

# Organize files by Sentrix ID, moving them into respective folders
for(i in unique(s_entrix_id)) {

#create s_entrix_folders
 if (dir.exists(paste0(path,"/",i)) == F) {

    dir.create(paste0(path,"/",i))

}

# Move files to corresponding folders and rename them

    files.i <- files[grep(i,files)]
    lapply(files.i, function(x){
        system(paste0("mv ", path,"/",x, " ",
        path,"/",i,"/",x),
        )
  
    })
}


#create sample_sheet

sample_sheet <- data.frame("Sample_Name" = names,
"Sample_Well" = rep("",length(names)),
"Sample_Plate" = rep("",length(names)),
"Sample_Group" = rep("",length(names)),
"Pool_ID" = rep("",length(names)),
"Sentrix_ID"   = s_entrix_id,
"Sentrix_Position" = s_entrix
)

# Remove duplicate entries
sample_sheet <- sample_sheet[!duplicated(sample_sheet),]


# Write the sample sheet to a CSV file if it doesn't already exist
if(file.exists(paste0(path,"/Sample_sheet.csv")) == F){

write.table(sample_sheet,
paste0(path,"/Sample_sheet.csv"),
sep = ",",
row.names = F,
col.names = T
)
}
}

# Generate the sample sheet for Alcala IDAT files
generate.sample.sheet.alcala(path)

############### Fix Color Problem with Sample ###############
# Load the generated sample sheet


ss <- read.table(paste0(datadir, "/external_data/Alcala_LUNET/IDAT/Sample_sheet.csv"),
sep = ",",
header = T)


# Remove problematic sample (only one color channel is available for it)
ss <- ss[grep("201414140136_R01C01",ss$Sample_Name, invert = T),]


# Write the corrected sample sheet back to the file
write.table(ss, paste0(datadir, "/external_data/Alcala_LUNET/IDAT/Sample_sheet.csv"),
 sep=",", row.names = F, col.names=T)



############### Methylation Matrix Generation ###############
# Function to generate methylation data matrix from IDAT files

methylation_matrix(paste0(datadir, "/external_data/Alcala_LUNET/IDAT"),
arraytype = "EPIC",
outdir)


# Load the M-value methylation matrix generated from IDAT files
M.value <- read.table(paste0(outdir,
"/M_matrix.csv"), sep=",", row.names = 1, header = T, check.names = F)  


############### Annotations and Final Output ###############
# Load sample annotations (e.g., Ms_IDs) to update column names
anotation <- read.csv(paste0(datadir, "/external_data/Alcala_LUNET/methylation_targets_MS_IDs.csv"),
sep = ",",
row.names = 1,
header = T,
check.names = FALSE
)

# Match annotation to the M-value matrix by column names
anotation <- anotation[colnames(M.value),]

# Update the M-value matrix with the appropriate sample IDs from annotation
colnames(M.value) <- anotation$Ms_IDs

# Write the processed M-value methylation matrix to a TSV file
write.table(M.value,
paste0(outdir, "/Alcala_M.tsv"),
sep = "\t",
row.names = T,
col.names = NA)

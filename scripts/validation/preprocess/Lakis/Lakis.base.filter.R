####################################################
# Lakis Data Integration Script
#
# This script integrates clinical, expression, methylation, and mutation data 
# from the Lakis PANET dataset. It selects relevant patient data, processes 
# gene mutation information, and merges it into a final dataset.
#
# Requirements:
# - Tidyverse package for data manipulation.
# - mygene package for gene ID conversion.
# - Expression and methylation matrices.
# - Clinical data, somatic mutation data from Lakis PANET.
#
# Inputs:
# - base_Lakis_prev.tsv: Clinical dataset.
# - matrix_Lakis_expression.txt: Expression matrix for Lakis dataset.
# - Methylation_M_Lakis.tsv: Methylation matrix (M-values).
# - simple_somatic_mutation.open.PAEN-AU.tsv: Somatic mutation data.
#
# Outputs:
# - base_Lakis.tsv: Final integrated dataset combining clinical, expression, 
#   methylation, and mutation information.
#
# Required environment: 'geo.yml' (containing prerequisite libraries).
####################################################


library(tidyverse)


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

outdir <- paste0(wkdir, "/validation/Lakis")



####################################################
# Load Clinical Data
####################################################
# Load the clinical base data, which was previously processed and saved as 'base_Lakis_prev.tsv'
base <- read.table(paste0(outdir, "/base_Lakis_prev.tsv"),
sep = "\t",
header = T)


####################################################
# Load Expression and Methylation Data
####################################################
# Load expression matrix for Lakis dataset
expression <- read.table(paste0(outdir, "/matrix_Lakis_expression.txt"),
sep = "\t",
row.names = 1,
header = T)

# Load methylation matrix for Lakis dataset (M-values)
methylation <- read.table(paste0(outdir, "/Methylation_M_Lakis.tsv"),
sep = "\t",
row.names = 1,
header = T)



####################################################
# Filter Base Data for Patients Present in Expression/Methylation Datasets
####################################################
# Identify unique patient IDs from expression and methylation datasets
patients <- unique(c(colnames(expression), colnames(methylation)))

# Select clinical data for patients that are present in either the expression or methylation datasets
base <- base[rownames(base)[base$id_patient %in% patients],]

# Order the clinical data by patient ID for consistent processing
base <- base[order(base$id_patient),]

# Remove entries for normal tissues or samples not derived from primary tumors
base <- base[!base$primary.tumor %in% c(
    "Normal - tissue adjacent to primary", "Normal - blood derived",
    "Normal - solid tissue"

),]


# Sort the clinical data by 'Celularity' (percentage of tumor cellularity)
base <- base[order(base$Celularity),]

# Remove duplicate entries, keeping only one record per patient
base <- base[!duplicated(base$id_patient),]
############################################################


####################################################
# Process Somatic Mutation Data
####################################################
# Load somatic mutation data for the Lakis dataset

mutations <- read.table(paste0(datadir, "/external_data/Lakis_PANET/simple_somatic_mutation.open.PAEN-AU.tsv"),
sep = "\t",
header = T)
# Remove duplicate mutation records
mutations <- mutations[!duplicated(mutations$icgc_mutation_id),]

# Extract unique gene names from mutation data
genes <- unique((mutations$gene_affected))

# Convert gene Ensembl IDs to gene symbols using the mygene package
genesID <-mygene::queryMany(genes, scopes="ensembl.gene", fields="symbol", species="human")
genesID <- as.data.frame(genesID[,c("query","symbol")])

# Merge mutation data with gene symbols
mutations <- merge(mutations,genesID, by.x = "gene_affected", by.y = "query", all.x = T)


####################################################
# Filter for Mutations in Genes of Interest
####################################################
# Define a list of genes of interest for analysis
interest <- c("MEN1", "ATRX", "DAXX", "PIK3CA", "PTEN", "TSC2", "VHL", "TP53", "RB1",
"CDKN1B", "APC", "KRAS", "BRAF", "PSIP1", "ARID1A", "SMARCA1", " SMARCA2", "SMARCA4",
"ARID2", "SMARCC2", "SMARCB1", "BCL11A", "SETD1B", "SETDB1", "NSD1", "KDM4A", "PHF8",
"JMJD1", "CBX6", "EZH1", "YY1", "BRWD3", "HDAC4", "STAG1", "NIPBL", "DICER", "ERCC6L",
"KMT2A", "KMT2C", "KMT2D", "ARID1A", "ARID1B", "ARID2", "NOTCH1", "NOTCH2", "NOTCH3",
"NOTCH4", "NF1", "INSR", "TSC2", "RICTOR", "PTPRN", "BEX1", "NOVA1")

# Filter mutation data for the genes of interest
mutations <- mutations[mutations$symbol %in% interest,]

# Add a 'mutated' column to indicate mutation status
mutations$mutated <- rep("mutated",nrow(mutations))

# Retain only relevant columns: sample ID, gene symbol, and mutation status
mutations <- mutations[,c("icgc_sample_id", "symbol", "mutated")]

# Remove duplicate mutation records
mutations <- mutations[!duplicated(mutations),]


####################################################
# Transform Mutation Data to Wide Format
####################################################
# Pivot mutation data to wide format, with gene symbols as column headers
mutations_wider <- mutations %>% pivot_wider(

    names_from = icgc_sample_id,
    values_from = mutated
)

# Transpose the matrix so that samples are rows and genes are columns
mutations_wider <- t(mutations_wider)

# Set column names and remove the first row (gene symbols)
colnames(mutations_wider) <- mutations_wider[1,]
mutations_wider <- mutations_wider[-1,]

# Replace missing mutation data with "WT" (wild type)
mutations_wider[is.na(mutations_wider)] <- "WT"


####################################################
# Merge Mutation Data with Clinical Base Data
####################################################
# Merge the mutation data with the clinical base data by sample ID
base <- merge(base,mutations_wider, by.x = "id_sample", 
by.y = "row.names", all.x = TRUE)

base <- base[,-6]


####################################################
# Write Final Integrated Dataset to Output File
####################################################
# Write the integrated clinical, expression, methylation, and mutation data to a TSV file
write.table(base,
paste0(outdir, "/base_Lakis.tsv"),
sep = "\t",
row.names = F,
col.names = T)


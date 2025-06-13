###########################################################
# Script: Validation Cohort Preparation for Multiple Studies
#
# Description:
#     This script prepares a combined clinical dataset from multiple studies
#     (Lakis, Alcala, Laddha, Patte) by loading, cleaning, and standardizing
#     their clinical metadata. The output is a unified table of samples with
#     associated attributes, ready for further validation analysis.
#
# Requirements:
#     - GEO environment with required packages.
#     - Clinical data for multiple studies in tsv files.
#
# Inputs:
#     - Clinical metadata for each study.
#     - Pre-processed expression and methylation data for filtering.
#
# Outputs:
#     - Unified clinical dataset in `Base_validation.tsv`.
#
# Required environment: `geo` (contains necessary libraries for processing).
###########################################################

############### Load Required Libraries ###############
library(tidyverse)
library(GEOquery)

############### Configuration and Directory Setup ###############
# Load configuration file to get working and data directory paths
config <- read.csv("config/config.tsv",
sep = "\t",
header = T)

# Extract key-value pairs for configuration settings
keys <- config$key
values <- config$value
names(values) <- keys

# Define working directory and data directory from configuration file
wkdir <- values["wkdir"]
datadir <- values["datadir"]

# Create necessary directories if they don't exist
if (dir.exists(paste0(wkdir, "/validation")) == F) {

    dir.create(paste0(wkdir, "/validation"))

}

if (dir.exists(paste0(wkdir, "/validation/base")) == F) {

    dir.create(paste0(wkdir, "/validation/base"))

}

# Define output directory for results
outdir <- paste0(wkdir, "/validation/base")


##################################################################
# Lakis Study Data Processing

# Load Lakis clinical base file

base_Lakis <- read.table(paste0(wkdir, "/validation/Lakis/base_Lakis.tsv"),
sep = "\t", row.names = 1, header = T)

# Reformat rownames and drop first column
rownames(base_Lakis) <- base_Lakis[,1]
base_Lakis <- base_Lakis[,-1]

# Convert relevant columns to factors and recode levels
base_Lakis$Gender <- as.factor(base_Lakis$Gender)
levels(base_Lakis$Gender) <- c("Female", "Male")

base_Lakis$tissue <- as.factor(base_Lakis$primary.tumor)
levels(base_Lakis$tissue) <- c("Primary.tumor", "Primary.tumor")


base_Lakis$primary.tumor <- as.factor(base_Lakis$primary.tumor)
levels(base_Lakis$primary.tumor) <- c("pancreas", "pancreas")

base_Lakis$Differenciation <- as.factor(base_Lakis$Differenciation)
levels(base_Lakis$Differenciation) <- c("NET", "NET", "NET","NET")

base_Lakis$Grade <- as.factor(base_Lakis$Grade)
levels(base_Lakis$Grade) <- c(NA, "G1", "G2")

base_Lakis$Stage <- as.factor(base_Lakis$Stage)
levels(base_Lakis$Stage) <- c(NA, "I", "I", "I","II","II","III","III",
"III","II","II","II","II", "IV","III","III","III")

# Adjust survival status and survival time to months
base_Lakis$Overall.survival <- base_Lakis$Overall.survival / 30.41

base_Lakis$Status <- ifelse(base_Lakis$Status == "alive",0,1)

# Label study source
base_Lakis$paper <- rep("Lakis_PANET", nrow(base_Lakis))

# Filter samples based on common availability in expression and methylation data
expression.Lakis <- read.table(paste0(wkdir, "/validation/Lakis/matrix_Lakis_expression.txt"),
sep = "\t",row.names = 1, header = T,nrows = 3)

methylation.Lakis <- read.table(paste0(wkdir, "/validation/Lakis/Methylation_M_Lakis.tsv"),
sep = "\t",
row.names = 1,
header = T,
nrows = 3
)

patients <- intersect(colnames(expression.Lakis),colnames(methylation.Lakis))

# Initialize base dataset
base_Lakis <- base_Lakis[rownames(base_Lakis) %in% patients,]

base <- base_Lakis

#############################################################################
# Alcala Study Data Processing

# Load Alcala clinical base file

base_Alcala <- read.table(paste0(datadir, "/external_data/Alcala_LUNET/Attributes.txt"),
sep = "\t", row.names = 1, header = T)

# Select specific columns and rename them
columns <- c("Age", "Sex", "Smoking_status", "Survival_months", "Survival_censor",
"Histopathology", "Stage_UICC", "cluster_LNET")

base_Alcala <- base_Alcala[,columns]

colnames(base_Alcala) <- c("Age", "Gender", "Smoking",
"Overall.survival","Status", "Differenciation", "Stage","Clusters")

# Convert relevant columns to factors and recode levels
base_Alcala$Gender <- as.factor(base_Alcala$Gender)
levels(base_Alcala$Gender) <- c("Female", "Male")

base_Alcala$Grade <- as.factor(base_Alcala$Differenciation)
levels(base_Alcala$Grade) <- c("G2",NA,"G3","G3", "G3", "G3","G3","G1")

base_Alcala$Differenciation <- as.factor(base_Alcala$Differenciation)
levels(base_Alcala$Differenciation) <- c("NET.atypical", "NET", "Carcinoma","Carcinoma",
"Carcinoma", "Carcinoma","Carcinoma","NET.typical")
# Remove carcinoma samples and set additional variables
base_Alcala <- base_Alcala[base_Alcala$Differenciation != "Carcinoma",]

base_Alcala$primary.tumor <- rep("lung", nrow(base_Alcala))

base_Alcala$tissue <- rep("Primary.tumor", nrow(base_Alcala))

base_Alcala$paper <- rep("Alcala_LUNET", nrow(base_Alcala))
# Recode survival status
base_Alcala$Status <- as.factor(base_Alcala$Status)
levels(base_Alcala$Status) <- c(0,1)

base_Alcala$Stage <- as.factor(base_Alcala$Stage)
levels(base_Alcala$Stage) <- c("I","I","I","I","I","I","II","II","III","III","IV","IV")

base_Alcala$Clusters <- as.factor(base_Alcala$Clusters)
levels(base_Alcala$Clusters) <- c("Carcinoid-A1", "Carcinoid-A2","Carcinoid-B")
# Filter samples based on availability in expression and methylation data
exp.Alcala <- read.table(paste0(wkdir, "/validation/Alcala/Alcala.expression.tsv"),
sep = "\t",row.names = 1, header = T,nrows = 3)

met.Alcala <- read.table(paste0(wkdir, "/validation/Alcala/Alcala_M.tsv"),
sep = "\t",row.names = 1, header = T,nrows = 3)

patients <- intersect(colnames(exp.Alcala),colnames(met.Alcala))

base_Alcala <- base_Alcala[rownames(base_Alcala) %in% patients,]

# Harmonize column structure with base dataset
columns_prev <- colnames(base)[!colnames(base) %in% colnames(base_Alcala)]

for (column in columns_prev){

base_Alcala[,column] <- rep(NA, nrow(base_Alcala))


}

columns_next <- colnames(base_Alcala)[!colnames(base_Alcala) %in% colnames(base)]

for (column in columns_next){

base[,column] <- rep(NA, nrow(base))


}

base_Alcala <- base_Alcala[,colnames(base)]




# Append Alcala study data to base
base <- rbind(base,base_Alcala)
############################################################################################
############################################################################################
# Laddha Study Data Preparation (Commented Out)
# Annotations to load clinical data for Laddha (commented due to download issues)
# Clinical data is already prepared in the data folder.

# GSE118132 <- getGEO("GSE118132", GSEMatrix =TRUE, getGPL=FALSE)
#Timeout was reached: [] Operation timed out after 120000 milliseconds
# and several times are needed to download phenodata. Thus clinical data
# are already prepare in data folder  

#GSE118131 <- getGEO("GSE118131", GSEMatrix =TRUE)
#GSE118131 <- GSE118131[[1]]
#base_Laddha <- pData(phenoData(GSE118131))

#srr <- c("SRR7646229" ,"SRR7646230" ,"SRR7646231" ,"SRR7646232"
# ,"SRR7646233" ,"SRR7646234" ,"SRR7646235" ,"SRR7646236" ,"SRR7646237"
#  ,"SRR7646238" ,"SRR7646239" ,"SRR7646240" ,"SRR7646241" ,"SRR7646242" 
#  ,"SRR7646243" ,"SRR7646244" ,"SRR7646245" ,"SRR7646246" ,"SRR7646247" 
#  ,"SRR7646248" ,"SRR7646249" ,"SRR7646250" ,"SRR7646251" ,"SRR7646252" ,"SRR7646253"
# ,"SRR7646254" ,"SRR7646255" ,"SRR7646256" ,"SRR7646257" ,"SRR7646258")

#base_Laddha$srr <- srr 

#columns <- c("title","srr","age (years):ch1", "cell type:ch1", "diagnosis:ch1"
#, "gender:ch1", "lcs subtype:ch1", "phenotype:ch1", "tumor location:ch1")


#base_Laddha <- base_Laddha[,columns]

#GSE118132 <- getGEO("GSE118132", GSEMatrix =TRUE, getGPL=FALSE)
#GSE118132 <- GSE118132[[1]]
#base_Laddha2 <- pData(phenoData(GSE118132))
#base_Laddha2 <- base_Laddha2[,c("title","geo_accession")]

#base_Laddha <- merge(base_Laddha, base_Laddha2,by = "title",all = F)
#base_Laddha <- base_Laddha[,-1]
#rownames(base_Laddha) <- base_Laddha$srr
#base$geo <- rep(NA,nrow(base))
#base[base_Laddha$srr,"geo"] <- base_Laddha$geo_accession
#base_Laddha <- base_Laddha[,-1]

#colnames(base_Laddha) <- c("Age", "tissue", "Differenciation", "Gender", "Clusters", "primary.tumor", "tumor.location","geo") 

#base_Laddha$tissue <- rep("Primary.tumor",nrow(base_Laddha))
#base_Laddha$Differenciation <- as.factor(base_Laddha$Differenciation)
#levels(base_Laddha$Differenciation) <- c("NET.atypical", "NET.typical")

#base_Laddha$Gender <- as.factor(base_Laddha$Gender)
#levels(base_Laddha$Gender) <- c("Female", "Male")


#base_Laddha$Clusters <- as.factor(base_Laddha$Clusters)
#levels(base_Laddha$Clusters) <- c("LC1", "LC2", "LC3")

#base_Laddha$primary.tumor <- rep("lung", nrow(base_Laddha))
#base_Laddha$paper <- rep("Laddha_LUNET",nrow(base_Laddha))

#write.table(base_Laddha, "data/external_data/Laddha_LUNET/base_Laddha.tsv",
#sep = "\t", row.names = T, col.names = NA)

# Load pre-prepared clinical data for Laddha
base_Laddha <- read.table("data/external_data/Laddha_LUNET/base_Laddha.tsv",
sep = "\t",
row.names = 1,
header = T)

# Harmonize column structure with base dataset
columns_prev <- colnames(base)[!colnames(base) %in% colnames(base_Laddha)]

for (column in columns_prev){

base_Laddha[,column] <- rep(NA, nrow(base_Laddha))


}

columns_next <- colnames(base_Laddha)[!colnames(base_Laddha) %in% colnames(base)]

for (column in columns_next){


base[,column] <- rep(NA, nrow(base))


}

base_Laddha <- base_Laddha[,colnames(base)]
# Append Laddha study data to base
base <- rbind(base,base_Laddha)

rownames(base)[base$paper == "Laddha_LUNET"] <- rownames(base_Laddha)

######################################################################################
# Patte Study Data Processing

# Load Patte sample-level clinical data
base_Patte <- read.table(paste0(datadir,"/external_data/SiNETs_Patte/METADATA_Samples_TNEi_for_paper_parcial.csv"),
sep = "\t",
row.names = 1,
header = T,
check.names = FALSE
)

# Filter and clean based on tissue and duplication
base_Patte <- base_Patte[base_Patte$tissue %in% c("T","MN","Liver","LN"),]
base_Patte$tissue <- factor(base_Patte$tissue, levels = c("T","MN","Liver","LN"))
base_Patte <- base_Patte[order(base_Patte$tissue),]

base_Patte <- base_Patte[!duplicated(base_Patte$id.patient),]

# Select and rename key variables
base_Patte <- base_Patte[,c( "id.patient", "tissue", "grade","Ki67","pT", "chr18del", "chr4gain", "chr5gain", "chr7gain", "chr10gain", "chr14gain", "chr20gain")]
colnames(base_Patte) <- c("id_patient", "tissue", "Grade", "Ki67", "Stage", "chr18del", "chr4gain", "chr5gain", "chr7gain", "chr10gain", "chr14gain", "chr20gain")
base_Patte$Sample <- rownames(base_Patte)

# Load and merge patient-level metadata
base_patient <- read.csv(paste0(datadir, "/external_data/SiNETs_Patte/Metadata.csv"),
sep = "\t",
row.names = 1,
header = T,
check.names = FALSE,
skip = 1

)

colnames(base_patient)

base_patient <- base_patient[,c("diagnostic.age", "primary.site", "OS_status", "time.surv")]

colnames(base_patient) <- c("Age", "primary.tumor", "Status", "Overall.survival")

# Merge sample- and patient-level metadata
base_Patte <- merge(base_Patte, base_patient, by.x = "id_patient", by.y = "row.names", all = T)

rownames(base_Patte) <- base_Patte$Sample

# Format variables
base_Patte$Overall.survival <- base_Patte$Overall.survival / 30.41
base_Patte$primary.tumor <- "small.intestine"


levels(base_Patte$tissue) <- c("Primary.tumor", "Liver", "Lymph.node", "Mesenteric")

base_Patte$Stage <- as.factor(base_Patte$Stage)
levels(base_Patte$Stage) <- c("I", "II", "III", "IV")

base_Patte$Status <- as.factor(base_Patte$Status)
base_Patte$paper <- rep("Patte_SINET", nrow(base_Patte))



# Filter samples based on availability in expression and methylation data
expression.Patte <- read.table(paste0(wkdir, "/validation/Patte/Patte.expression.tsv"),
sep = "\t",row.names = 1, header = T,nrows = 3)

methylation.Patte <- read.table(paste0(wkdir, "/validation/Patte/Patte_M.tsv"),
sep = "\t",row.names = 1, header = T,nrows = 3)

patients <- intersect(colnames(expression.Patte),colnames(methylation.Patte))

patients <- patients[patients %in% rownames(base_Patte)]

base_Patte <- base_Patte[patients,]

# Harmonize column structure with base dataset

columns_prev <- colnames(base)[!colnames(base) %in% colnames(base_Gibert)]

for (column in columns_prev){

base_Gibert[,column] <- rep(NA, nrow(base_Gibert))


}

columns_next <- colnames(base_Gibert)[!colnames(base_Gibert) %in% colnames(base)]

for (column in columns_next){

base[,column] <- rep(NA, nrow(base))


}


base_Patte <- base_Patte[,colnames(base)]

# Append Patte study data to base
base <- rbind(base,base_Patte)

# Adjust CNV columns if applicable
base$chr14gain[base$cnv == "gain14"] <- "YES"
base$chr18del[base$cnv == "loss18"] <- "YES"



##############################################################################################
# Output final combined dataset
write.table(base, paste0(outdir, "/Base_validation.tsv"),
sep = "\t",
row.names = T,
col.names = NA)

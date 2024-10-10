############################################################
# Script: Gene Expression Data Processing for Karpathakis Dataset
#
# Description:
#     This script processes gene expression data for the Karpathakis dataset. 
#     It normalizes the data, selects genes with high variance, adjusts sample names, 
#     and filters out irrelevant or duplicate patient samples. The final processed 
#     data is saved in tab-separated value (TSV) files for further analysis.
############################################################
# Required environment: 'geo.yml' (containing prerequiste libraries).

# Load required libraries for expression data processing and analysis
library(edgeR)
library(limma)


############### Load Configuration ###############
# Read configuration file to define paths
config <- read.csv("config/config.tsv",
sep = "\t",
header = T)

# Extract configuration keys and values
keys <- config$key
values <- config$value
names(values) <- keys

# Set working directory and data directory from config
wkdir <- values["wkdir"]
datadir <- values["datadir"]


############### Directory Setup ###############
# Create directories for validation results if they don't exist
if (dir.exists(paste0(wkdir, "/validation")) == F) {

    dir.create(paste0(wkdir, "/validation"))

}

if (dir.exists(paste0(wkdir, "/validation/Karpathakis")) == F) {

    dir.create(paste0(wkdir, "/validation/Karpathakis"))

}

# Set output directory for processed Karpathakis data
outdir <- paste0(wkdir, "/validation/Karpathakis")

############### Load Data ###############
# Load baseline patient information

base <- read.table(paste0(datadir,"/external_data/Karpathakis_SINET/Base_Karpathakis.tsv"),
sep = "\t",
row.names = 1,
header = T)

# Load gene expression data
base_exp <- read.table(paste0(datadir,"//external_data/Karpathakis_SINET/data_expression.tsv"),
sep = "\t",
row.names = 1,
header = T,
skip =1)

# Load additional expression table with p-values
exp <- read.csv(paste0(datadir,"/external_data/Karpathakis_SINET/AKGINETDASL_allTvsallN_expressiontable_normback_geneandpval_030414.txt"),
sep = "\t",
header = T,
check.names=F)


############### Preprocessing ###############
# Extract p-values and clean up column names
p.val <- exp[,grep("Detection Pval",colnames(exp))]
colnames(p.val) <- gsub("\\.Detection Pval","",colnames(p.val))


# Compute the median p-value for each sample

m.p.value <- apply(p.val,2,median)

# Extract signal intensities and clean up column names
exp <- exp[,c(1,grep("Signal",colnames(exp)))]
colnames(exp) <- gsub("\\.AVG_Signal","",colnames(exp))

# Select genes based on variance and remove duplicates
genes <- exp$SYMBOL

symbol_genes <- limma::alias2SymbolTable(genes)

df_expression <- data.frame("Expression" = genes <- exp$SYMBOL,
"Symbol" = symbol_genes)

# Merge gene symbols into the expression matrix
exp <- merge(df_expression,exp,by.x = "Expression",
by.y = "SYMBOL", all.y = T)

exp$Symbol[is.na(exp$Symbol)] <- exp$Expression[is.na(exp$Symbol)] 

# Select the most variable genes
var <- apply(exp[3:ncol(exp)],1,var)
exp <- exp[order(var,decreasing = T),]
exp <- exp[!duplicated(exp$Symbol),]
rownames(exp) <- exp$Symbol
exp <- exp[,3:ncol(exp)] # Remove non-expression columns


############### Adjust Sample IDs ###############
# Clean up patient IDs in the base dataset

base$id_patient <- gsub("Genomic DNA from","",base$id_patient)
base$id_patient <- gsub(" ","",base$id_patient)

# Clean up sample names in the expression data
colnames(exp) <- gsub("_A","A",colnames(exp))
colnames(exp) <- gsub("_B","B",colnames(exp))
colnames(exp) <- gsub("_06","06",colnames(exp))
colnames(exp) <- gsub("_H","H",colnames(exp))
colnames(exp) <- gsub("_1","1",colnames(exp))
colnames(exp) <- gsub("_IIA","IIA",colnames(exp))
colnames(exp) <- gsub("_VIC","VIC",colnames(exp))

# Additional specific patient ID corrections
colnames(exp) <- gsub("033349C_T1","033349CT1_T",colnames(exp))
colnames(exp) <- gsub("033349C_T2","033349CT2",colnames(exp))
colnames(exp) <- gsub("B052175A2_T","B052175A23081_T",colnames(exp))
colnames(exp) <- gsub("PH1015539A03_T","PH1015539A03691_T",colnames(exp))
colnames(exp) <- gsub("9572A6_T","9572A6308_T",colnames(exp))
colnames(exp) <- gsub("10D22894A1_T","10D22894A169_T",colnames(exp))
colnames(exp) <- gsub("1166106G_T","1166106G69_T",colnames(exp))
colnames(exp) <- gsub("1081741A_T","1081741A308_T",colnames(exp))
colnames(exp) <- gsub("165092_T","0165092_T",colnames(exp))
colnames(exp) <- gsub("08C106462_T","08C10646269_T",colnames(exp))
colnames(exp) <- gsub("81755_T","8175834_T",colnames(exp))
colnames(exp) <- gsub("56235_T","5623834_T",colnames(exp))
colnames(exp) <- gsub("1186316_T","1186316308A_T",colnames(exp))
colnames(exp) <- gsub("PH1015539A03_N","PH1015539A03692_N",colnames(exp))
colnames(exp) <- gsub("066042A2_T1","066042A2_T",colnames(exp))


############### Update Metadata ###############
# Align patient metadata with expression data
base_exp <- base_exp[colnames(exp),]
base_exp$p.val <- m.p.value
base_exp$ID <- rownames(base_exp)

# More specific corrections for column names
colnames(exp) <- gsub("1166106G69_T","116106G69_T",colnames(exp))
colnames(exp) <- gsub("68279A3_T","068279A3_T",colnames(exp))
rownames(base_exp) <- colnames(exp)

# Match and filter patients
base <- base[base$id_patient %in% colnames(exp),]
exp <- exp[,colnames(exp) %in% base$id_patient]
base_exp <- base_exp[rownames(base_exp) %in% base$id_patient,]

# Adjust row names in the base data
base$geo.names <- rownames(base)
rownames(base) <- base$id_patient
base <- base[rownames(base_exp),]

colnames(exp) <- base$geo.names
base$p.val <- base_exp$p.val

# Filter out irrelevant or duplicate patient samples

rep_pat <- c("0414821A4_N","0612608A8_N","066042A2_T","066042A6_T","1015089A3_N",
"12A5867A1_N","9913069VIC_N","H112177A2_N","PH1015539A03692_N")

base <- base[! base$id_patient %in% rep_pat,]

# Reorder expression data based on patient IDs
exp <- exp[,base$geo.names]



############### Save Processed Data ###############
# Write the processed expression data to a TSV file
write.table(exp,
paste0(outdir,"/Karpathakis.expression.tsv"),
sep ="\t",
row.names = T,
col.names = NA)

# Write the updated patient metadata to a TSV file
write.table(base,
paste0(outdir,"/base_Karpathakis.expression.tsv"),
sep = "\t",
row.names = T,
col.names = NA)

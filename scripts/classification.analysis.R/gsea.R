####################################################
# Script: GSEA Analysis for NS Subtypes
#
# Description:
#      This script performs Gene Set Enrichment Analysis (GSEA) for three different subtypes (NS1, NS2, NS3)
#      using transcriptome expression data. The script generates dot plots for each subtype comparison
#      against other samples and saves the analysis results in designated output files.
#
# Requirements:
#      - Predefined environment with GSEA and related functions from `main.functions.R`, `gsea.functions.R`, and `common.functions.R`.
#      - Configuration file (`config.tsv`) containing directory paths for working and data directories.
#      - Clinical data file (`clinical_data.txt`) and gene expression data file (`expression.qual.txt`).
#      - Gene sets in GMT format for enrichment analysis.
####################################################
# Required environment: gsea.yml.

# Load required libraries
library(fgsea)
library(GSVA)
library(ggplot2)
library(GSEABase)
library(ComplexHeatmap)
library(circlize)

# Load custom functions for GSEA and common clustering operations
source("functions/main.functions.R")
source("functions/clustering/gsea.functions.R")
source("functions/clustering/common.functions.R")


############### Configuration and Directory Setup ###############
# Read the configuration file (config.tsv) to retrieve the directory paths

config <- read.csv("config/config.tsv",
sep = "\t",
header = T)
# Assign configuration keys and values
keys <- config$key
values <- config$value
names(values) <- keys

# Define working directory and data directory based on configuration values

wkdir <- values["wkdir"]
datadir <- values["datadir"]

# Check if the output directory for GSEA analysis exists, create if it doesn't
if (dir.exists(paste0(wkdir, "/clustering/gsea.NS")) == F) {

    dir.create(paste0(wkdir, "/clustering/gsea.NS"))

}

# Set output directory for GSEA results
outdir <- paste0(wkdir, "/clustering/gsea.NS")



############### Data Loading ###############
# Load clinical data from the specified file
base <- read.csv(paste0(datadir, "/clinical_data.txt"),
                 sep = "\t",
                 header = T,
                 row.names =1)

# Load transcriptome expression data from the specified file
exp <- read.table(paste0(wkdir,
"/preprocess/transcriptome/expression.qual.txt"),
    sep = "\t",
    header = T,
    row.names = 1,
    check.names = F
)



############### Differential Expression Data Loading ###############
# Load differential expression data for NS1 vs Other group comparison

NS1 <- read.table(paste0(wkdir,
"/clustering/deg/TopTable_deg_NS1_GroupNS1-GroupOther.txt"),
sep = "\t",
row.names = 1,
header = T
)
# Create ranked list based on t-statistics for NS1 vs Others
rnk1 <- NS1$t 
names(rnk1) <- rownames(NS1)
rnk1 <- rnk1[order(rnk1, decreasing = T)]

# Repeat the same for NS2 vs Others
NS2 <- read.table(paste0(wkdir,
"/clustering/deg/TopTable_deg_NS2_GroupNS2-GroupOther.txt"),
sep = "\t",
row.names = 1,
header = T
)


rnk2 <- NS2$t 
names(rnk2) <- rownames(NS2)
rnk2 <- rnk2[order(rnk2, decreasing = T)]


# Repeat the same for NS3 vs Others
NS3 <- read.table(paste0(wkdir,
"/clustering/deg/TopTable_deg_NS3_GroupNS3-GroupOther.txt"),
sep = "\t",
row.names = 1,
header = T
)


rnk3 <- NS3$t 
names(rnk3) <- rownames(NS3)
rnk3 <- rnk3[order(rnk3, decreasing = T)]


############### GSEA Gene Set Files Loading ###############
# Load GMT files containing gene sets for GSEA analysis
gmts <- list.files(paste0(wkdir,"/genesets"))
gmts <- gmts[grep(".gmt", gmts)]
gmts <- paste0(wkdir, "/genesets/", gmts)



############### GSEA Analysis for Each Subtype ###############
# Loop through each GMT file and perform GSEA for NS1, NS2, and NS3 subtypes
for (gmt in gmts) {

    # Extract the gene set name from the GMT file path
        name <- gsub(".*/|.gmt", "", gmt)
    # Perform GSEA for NS1 vs Others
        gsea1 <- gsea(rnk = rnk1,
                       gmt = gmt
                       , title = paste0(outdir, "/", name,
                    "_dotplot_NS1vsOthers.pdf"))

         gsea1$group1 <- "NS1"
         gsea1$group2 <- "Other"

    # Save GSEA results for NS1 vs Others   
        write.table(gsea1
                    , paste0(outdir, "/", name,
                    "_fgsea_G1vsG2.txt")
                    ,sep = '\t'
                    ,row.names = F
                    ,col.names = T)


    # Perform GSEA for NS2 vs Others
                gsea2 <- gsea(rnk = rnk2,
                       gmt = gmt
                       , title = paste0(outdir, "/", name,
                    "_dotplot_NS2vsOthers.pdf"))

         gsea2$group1 <- "NS2"
         gsea2$group2 <- "Other"

    # Save GSEA results for NS2 vs Others
        write.table(gsea2
                    , paste0(outdir, "/", name,
                    "_fgsea_NS2vsOthers.txt")
                    ,sep = '\t'
                    ,row.names = F
                    ,col.names = T)


    # Perform GSEA for NS3 vs Others
                gsea3 <- gsea(rnk = rnk3,
                       gmt = gmt
                       , title = paste0(outdir, "/", name,
                    "_dotplot_NS3vsOthers.pdf"))
         
         gsea3$group1 <- "NS3"
         gsea3$group2 <- "Other"


    # Save GSEA results for NS3 vs Others
        write.table(gsea3
                    , paste0(outdir, "/", name,
                    "_fgsea_NS3vsOthers.txt")
                    ,sep = '\t'
                    ,row.names = F
                    ,col.names = T)

        }

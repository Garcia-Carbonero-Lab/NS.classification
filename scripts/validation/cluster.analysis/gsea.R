###############################################################################
# Script: Gene Set Enrichment Analysis (GSEA)
#
# Description:
#     This script performs Gene Set Enrichment Analysis (GSEA) for subtypes 
#     (NS1, NS2, NS3) and generates ranked gene sets for each comparison.
#
# Inputs:
#     - DEG result files: TopTable_deg_NS1_GroupNS1-GroupOther.txt, 
#       TopTable_deg_NS2_GroupNS2-GroupOther.txt, TopTable_deg_NS3_GroupNS3-GroupOther.txt
#     - Gene sets in GMT format: located in /genesets/ directory
#     - Configuration file: config.tsv
#
# Outputs:
#     - GSEA results for NS1, NS2, and NS3
#     - Dot plots of enriched pathways (PDF format)
#     - fgsea tables for each comparison
#
# Required environment: `gsea`
#
###############################################################################
# Load required libraries for GSEA and data visualization
library(fgsea)
library(GSVA)
library(ggplot2)
library(GSEABase)

# Load custom functions for GSEA and general operations
source("functions/main.functions.R")
source("functions/clustering/gsea.functions.R")
source("functions/clustering/common.functions.R")

# Read the configuration file for directories
config <- read.csv("config/config.tsv",
sep = "\t",
header = T)

keys <- config$key
values <- config$value
names(values) <- keys
# Define working and data directories
wkdir <- values["wkdir"]
datadir <- values["datadir"]

# Create GSEA output directory if it doesn't exist
if (dir.exists(paste0(wkdir, "/validation/clustering.analysis/gsea")) == F) {

    dir.create(paste0(wkdir, "/validation/clustering.analysis/gsea"))

}

outdir <- paste0(wkdir, "/validation/clustering.analysis/gsea")

# Load differential expression results for each subtype
NS1 <- read.table(paste0(wkdir,"/validation/clustering.analysis/deg/TopTable_deg_NS1_GroupNS1-GroupOther.txt"),
sep = "\t",
row.names = 1,
header = T
)
# Create ranked list of genes based on t-statistic for NS1
rnk1 <- NS1$t 
names(rnk1) <- rownames(NS1)
rnk1 <- rnk1[order(rnk1, decreasing = T)]

NS2 <- read.table(paste0(wkdir,"/validation/clustering.analysis/deg/TopTable_deg_NS2_GroupNS2-GroupOther.txt"),
sep = "\t",
row.names = 1,
header = T
)
# Create ranked list of genes based on t-statistic for NS2
rnk2 <- NS2$t 
names(rnk2) <- rownames(NS2)
rnk2 <- rnk2[order(rnk2, decreasing = T)]

NS3 <- read.table(paste0(wkdir,"/validation/clustering.analysis/deg/TopTable_deg_NS3_GroupNS3-GroupOther.txt"),
sep = "\t",
row.names = 1,
header = T
)

# Create ranked list of genes based on t-statistic for NS3
rnk3 <- NS3$t 
names(rnk3) <- rownames(NS3)
rnk3 <- rnk3[order(rnk3, decreasing = T)]

# Load gene sets in GMT format from the genesets directory
gmts <- list.files(paste0(wkdir,"/genesets"))
gmts <- gmts[grep(".gmt", gmts)]
gmts <- paste0(wkdir, "/genesets/", gmts)

# Loop through each gene set file (GMT) and perform GSEA for each subtype
for (gmt in gmts) {

    # Extract the name of the gene set file (without the .gmt extension)
        name <- gsub(".*/|.gmt", "", gmt)
    # Perform GSEA for NS1
        gsea1 <- gsea(rnk = rnk1,
                       gmt = gmt
                       , title = paste0(outdir, "/", name,
                    "_dotplot_NS1vsOthers.pdf"))

         gsea1$group1 <- "NS1"
         gsea1$group2 <- "Other"
    # Save GSEA results for NS1   
        write.table(gsea1
                    , paste0(outdir, "/", name,
                    "_fgsea_G1vsG2.txt")
                    ,sep = '\t'
                    ,row.names = F
                    ,col.names = T)
    # Perform GSEA for NS2
                gsea2 <- gsea(rnk = rnk2,
                       gmt = gmt
                       , title = paste0(outdir, "/", name,
                    "_dotplot_NS2vsOthers.pdf"))

         gsea2$group1 <- "NS2"
         gsea2$group2 <- "Other"
    # Save GSEA results for NS2
        write.table(gsea2
                    , paste0(outdir, "/", name,
                    "_fgsea_NS2vsOthers.txt")
                    ,sep = '\t'
                    ,row.names = F
                    ,col.names = T)
    # Perform GSEA for NS3      
                gsea3 <- gsea(rnk = rnk3,
                       gmt = gmt
                       , title = paste0(outdir, "/", name,
                    "_dotplot_NS3vsOthers.pdf"))
         
         gsea3$group1 <- "NS3"
         gsea3$group2 <- "Other"
        
    # Save GSEA results for NS3        
        write.table(gsea3
                    , paste0(outdir, "/", name,
                    "_fgsea_NS3vsOthers.txt")
                    ,sep = '\t'
                    ,row.names = F
                    ,col.names = T)

        }
####################################################################

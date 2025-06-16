###############################################################################
# Script for comparing two MOFA models (e.g., discovery vs. validation)
#
# Purpose:
# - Compare the weights (feature importance) of expression and methylation layers
#   between two different MOFA models (e.g., trained on different datasets or splits).
#
# Inputs:
# - Two MOFA models: 'mod1' (discovery) and 'mod2' (validation/combination)
# - Methylation annotation file
# - Configuration file with paths
#
# Outputs:
# - Correlation matrices of feature weights (expression and methylation)
# - Correlation plots (PNG/PDF) for expression and methylation comparisons
#
# Required environment: 'MOFA.compare'
###############################################################################
# Load required libraries
library(MOFA2)
library(ggplot2)
library(tidyverse)
library(ggcorrplot)


# Load custom helper functions for MOFA and general analysissource("functions/MOFA/MOFA.functions.R")
source("functions/main.functions.R")


# ------------------------------------------------------------------------------
# Read configuration file containing paths
# ------------------------------------------------------------------------------

config <- read.csv("config/config.tsv",
sep = "\t",
header = T)

keys <- config$key
values <- config$value
names(values) <- keys

# Extract working and data directories from config
wkdir <- values["wkdir"]
datadir <- values["datadir"]


# ------------------------------------------------------------------------------
# Create output folder for storing results
# ------------------------------------------------------------------------------

if (dir.exists(paste0(wkdir, "/validation/MOFA")) == F) {

    dir.create(paste0(wkdir, "/validation/MOFA"))

}


outdir <- paste0(wkdir, "/validation/MOFA")



# ------------------------------------------------------------------------------
# Load MOFA models
# ------------------------------------------------------------------------------

mod1 <- load_model(paste0(wkdir, "/MOFA/model/Mofa.model.hdf5"))
mod2 <- load_model(paste0(wkdir, "/validation/MOFA/model/Mofa.model.combination.hdf5"))


# ------------------------------------------------------------------------------
# Load and process EPIC methylation annotation
# ------------------------------------------------------------------------------
epic <- read.csv(paste0(datadir,
"/methylome/anotation/infinium-methylationepic-v-1-0-b5-manifest-file.csv"),
                       sep = ",",
                       header = T,
                       skip = 7)


# Function to process methylation annotation from EPIC manifest
create.met.anot <- function(epic) {

genes <-   strsplit(as.character(epic$UCSC_RefGene_Name),';')
genes[lengths(genes) == 0] <- NA
group <- strsplit(as.character(epic$UCSC_RefGene_Group),';')
group[lengths(group) == 0] <- NA

met_anot <- data.frame('ID' = rep(epic$IlmnID, sapply(genes, length))
    , "genes" = unlist(genes)
    , "CpG_chrm" = rep(epic$Chromosome_36, sapply(genes, length))
    , "MAPINFO" = rep(epic$MAPINFO, sapply(genes, length))
    , "CpG_Island_Name" = rep(epic$UCSC_CpG_Islands_Name, sapply(genes, length))
    , "CHR" = rep(epic$CHR, sapply(genes, length))
    , "UCSC_RefGene_Group" = unlist(group)
    , "ISLAND" = rep(epic$Relation_to_UCSC_CpG_Island, sapply(genes, length))
    , "FGROUP" = rep(epic$Regulatory_Feature_Group, sapply(genes, length))
    , "DMR" = rep(epic$DMR, sapply(genes, length)))

met_anot <- met_anot[!duplicated(paste(met_anot$ID,met_anot$genes)),]

    met_anot <- met_anot %>% group_by(ID,CpG_chrm,MAPINFO,CpG_Island_Name,ISLAND,FGROUP, DMR) %>% 
        dplyr::summarise(genes = toString(genes),
        UCSC_RefGene_Group = toString(UCSC_RefGene_Group))
    
    met_anot$genes <- gsub(', ','_',met_anot$genes)
        
    
    rownames(met_anot) <- met_anot$ID
    return(met_anot)
}




# Apply annotation processing
met_anot <- create.met.anot(epic)

# Annotate both MOFA models with processed methylation annotation
mod1 <- anotate.metilation(mod1, met_anot)
mod2 <- anotate.metilation(mod2, met_anot)


# ------------------------------------------------------------------------------
# Compare weights for EXPRESSION layer
# ------------------------------------------------------------------------------

exp_mod1 <- get_weights(mod1)$Expression
exp_mod2 <- get_weights(mod2)$Expression

# Align feature sets between models
exp_mod1 <- exp_mod1[rownames(exp_mod1) %in% rownames(exp_mod2),]
exp_mod2 <- exp_mod2[row.names(exp_mod2),]

# Calculate correlation matrix between expression weights of both models
corr.exp <- cor(exp_mod1,exp_mod2)

# Save correlation matrix
write.table(corr.exp,
paste0(outdir,"/Corr.expression.validation.vs.discovery.txt"),
sep = "\t",
row.names = T,
col.names = NA)


# ------------------------------------------------------------------------------
# Compare weights for METHYLATION layer
# ------------------------------------------------------------------------------
met_mod1 <- get_weights(mod1)$Methylation
met_mod2 <- get_weights(mod2)$Methylation

# Align features
met_mod1 <- met_mod1[rownames(met_mod1) %in% rownames(met_mod2),]
met_mod2 <- met_mod2[row.names(met_mod1),]

# Calculate correlation matrix for methylation weights
corr.met <- cor(met_mod1,met_mod2)

# Save correlation matrix
write.table(corr.met,
paste0(outdir,"/Corr.methylation.validation.vs.discovery.txt"),
sep = "\t",
row.names = T,
col.names = NA)

# ------------------------------------------------------------------------------
# Plot correlation matrix for EXPRESSION weights
# ------------------------------------------------------------------------------

  p <- ggcorrplot(corr.exp, method = "circle", lab = T, pch.cex = 5)
        p + xlab("Model1") +
            ylab("Model2") + theme_classic() +
            theme(
            text = element_text(size = 15),
            legend.key.size = unit(0.5, "cm"),
            legend.text = element_text(size = 12),
            axis.text.x = element_text(angle = 30, hjust = 1)
        )
        

# Save expression correlation plot

        ggsave(paste0(outdir, "/Correlation_plot.expression.png"),
            device = "png",
            dpi = 500, width = 8,
            height = 8
        )

        ggsave(paste0(outdir, "/Correlation_plot.expression.pdf"),
            device = "pdf",
            dpi = 500, width = 8,
            height = 8
        )

# ------------------------------------------------------------------------------
# Plot correlation matrix for METHYLATION weights
# ------------------------------------------------------------------------------


  p <- ggcorrplot(corr.met, method = "circle", lab = T, pch.cex = 5)
        p + xlab("Model1") +
            ylab("Model2") + theme_classic() +
            theme(
            text = element_text(size = 15),
            legend.key.size = unit(0.5, "cm"),
            legend.text = element_text(size = 12),
            axis.text.x = element_text(angle = 30, hjust = 1)
        )


# Save methylation correlation plot

        ggsave(paste0(outdir, "/Correlation_met.expression.png"),
            device = "png",
            dpi = 500, width = 8,
            height = 8
        )
        
        ggsave(paste0(outdir, "/Correlation_met.expression.pdf"),
            device = "pdf",
            dpi = 500, width = 8,
            height = 8
        )
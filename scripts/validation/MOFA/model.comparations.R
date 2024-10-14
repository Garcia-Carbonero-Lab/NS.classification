###############################################################################
# Script: Model Comparison for MOFA Models
#
# Description:
#     This script compares the latent factor structures of multiple MOFA (Multi-Omics 
#     Factor Analysis) models by calculating correlation matrices between corresponding 
#     factors from different models. It generates correlation plots and outputs text 
#     files with correlation values and comparative factor information.
#
# Requirements:
#     - Pre-trained MOFA models (discovery, validation, combination).
#     - Required libraries: MOFA2, ggplot2, tidyverse, ggcorrplot.
#
# Inputs:
#     - MOFA model files in HDF5 format.
#     - Configuration file for directory paths.
#
# Outputs:
#     - Correlation matrices in `.txt` format.
#     - Comparative factor assignments in `.txt` format.
#     - Correlation plots in `.pdf` format.
#
# Required environment: `MOFA.compare`
###############################################################################



###############################################################################
#This scripts need environment MOFA.compare 


############### Load Required Libraries ###############
library(MOFA2)
library(ggplot2)
library(tidyverse)
library(ggcorrplot)


############### Function: Compare MOFA Models ###############
# This function compares the latent factors of two MOFA models by calculating
# a Pearson correlation matrix and generating comparative factor information
Compare_models <- function(models, outdir, name, th = 0.80) {
        m1 <- models[[1]] # First model
        m2 <- models[[2]] # Second model
        z1 <- Reduce(rbind, m1@expectations$Z)
        z2 <- Reduce(rbind, m2@expectations$Z)

        # Find intersecting samples and factors between models
        p <- intersect(rownames(z1), rownames(z2))
        f <- intersect(colnames(z1), colnames(z2))
        
        # Calculate Pearson correlation between common samples and factors
        corr <- cor(z1[p, ], z2[p, ], method = "pearson")
        corr <- cor(z1[p, f], z2[p, f], method = "pearson")

        # Write correlation matrix to a file
        write.table(corr,
        paste0(outdir, "/", name, "_Correlation.models.txt"),
        sep = "\t",
        row.names = T,
        col.names = NA
        )

        # Assign factors based on maximum absolute correlation        
        corr.df <- as.data.frame(corr)
        fz1 <- apply(corr.df, 1, function(x){
        factor <- names(x)[abs(x) == max(abs(x))]
        res <- ifelse(abs(x)[factor] >= th, factor, "new")
        return(res)
        })
        fz2 <- apply(corr.df, 2, function(x){
        factor <- names(x)[abs(x) == max(abs(x))]
        res <- ifelse(abs(x)[factor] >= th, factor, "new")
        return(res)
        })

        # Create a data frame for comparative factor assignment between models        
        corr.df <- data.frame(nm1 = fz1, nm2 = fz2, row.names = rownames(corr))
        colnames(corr.df) <- c("Model1", "Model2")

        # Write comparative factor assignment to a file
        write.table(corr.df,
        paste0(outdir, "/", name, "_Comparative.models.txt"),
        sep = "\t",
        row.names = T,
        col.names = NA
        )

        # Generate correlation plot        
        p <- ggcorrplot(corr, method = "circle", lab = T, pch.cex = 5)
        p + xlab("Model1") +
            ylab("Model2") + theme_classic() +
            theme(
            text = element_text(size = 15),
            legend.key.size = unit(0.5, "cm"),
            legend.text = element_text(size = 12),
            axis.text.x = element_text(angle = 30, hjust = 1)
        )

        # Save the correlation plot as a PDF        
        ggsave(paste0(outdir, "/", name,"_correlation_plot.pdf"),
            device = "pdf",
            dpi = 500, width = 8,
            height = 8
        )
return(corr.df)
}



#load functions
# we read config file
config <- read.csv("config/config.tsv",
sep = "\t",
header = T)

keys <- config$key
values <- config$value
names(values) <- keys
###############
wkdir <- values["wkdir"]
datadir <- values["datadir"]

#We create output folder
if (dir.exists(paste0(wkdir, "/validation/MOFA/")) == F) {

    dir.create(paste0(wkdir, "/validation/MOFA/"))

}

if (dir.exists(paste0(wkdir, "/validation/MOFA/comparations")) == F) {

    dir.create(paste0(wkdir, "/validation/MOFA/comparations"))

}



# We indicate the out folder
outdir <- paste0(wkdir, "/validation/MOFA/comparations")
# Introduce clinical data and Factor value in samples_metadata
mod.discovery <-  load_model(paste0(wkdir, "/MOFA/model/Mofa.model.hdf5"))
mod.validation <- load_model(paste0(wkdir, "/validation/MOFA/model/Mofa.model.validation.hdf5"))
mod.combination <- load_model(paste0(wkdir, "/validation/MOFA/model/Mofa.model.combination.hdf5"))

models <-c(mod.discovery,mod.combination)

Compare_models(models,outdir, "Discovery.vs.Combination")

models <-c(mod.combination,mod.validation)
Compare_models(models,outdir, "Combination.vs.Discovery")

###############################################################################
# Script: Survival Analysis and Kaplan-Meier Estimation
#
# Description:
#     This script performs survival analysis using the Cox proportional hazards 
#     model and generates Kaplan-Meier survival curves for subtypes.
#
# Inputs:
#     - Clinical data: Base_validation.subtype.tsv
#     - Configuration file: config.tsv
#
# Outputs:
#     - Univariate and multivariate Cox model results
#     - Kaplan-Meier survival curves
#
# Required environment: `survival`
#
###############################################################################


#This script need the environment survival
library(survival)
library(survminer)
library(forestplot)
library(tidyverse)
library(aod)
library(forestmodel)

# Load custom functions for survival analysis
source("functions/clustering/survival.functions.R")

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


# Create output directory for survival analysis if it doesn't exist
if (dir.exists(paste0(wkdir, "/validation/clustering.analysis/survival")) == F) {

    dir.create(paste0(wkdir, "/validation/clustering.analysis/survival"))

}



outdir <- paste0(wkdir, "/validation/clustering.analysis/survival")

# Load clinical data
base <- read.table(paste0(wkdir, "/validation/base/Base_validation.subtype.tsv"),
                              sep = '\t', header = TRUE, row.names = 1)

# Modify the Stage column to create a binary Stage.IV column
base$Stage.IV <- base$Stage
base$Stage.IV[base$Stage %in% c("I","II","III")] <- "No"
base$Stage.IV[base$Stage == "IV"] <- "Yes"

# Create univariate Cox models for features known to be related to Overall Survival (OS)
for (feature in c("Subtype.validation",
"Grade",
"Stage.IV",
"Age",
"primary.tumor",
"Gender")) {

Hazard_Cox(base,
class = feature,
OS = "Overall.survival",
OS_event = "Status",
covariates = NULL,
outdir = outdir,
flag = paste0("Univariate_",feature))
}

# Create a multivariate Cox model including all relevant covariates
Hazard_Cox(base,
class = "Subtype.validation",
OS = "Overall.survival",
OS_event = "Status",
covariates = c("Grade","Stage.IV", "Age", "primary.tumor", "Gender"),
outdir = outdir,
flag = "Multivariate")

# Create a multivariate Cox model excluding Grade
Hazard_Cox(base,
class = "Subtype.validation",
OS = "Overall.survival",
OS_event = "Status",
covariates = c("Stage.IV", "Age", "primary.tumor", "Gender"),
outdir = outdir,
flag = "Multivariate_sn_Grade")


# Create Kaplan-Meier survival curves for subtypes
Kapplan_Meyer(cut_survival(base, 150,"Overall.survival", "Status"),
class = "Subtype.validation",
OS = "Overall.survival",
OS_event = "Status",
outdir = outdir,
flag = "cut_150",
custcolor = c("gold", "cornflowerblue", "tomato1")
)

###################################################
# Script: Survival Analysis for Cancer Subtypes
#
# Description:
#      This script performs survival analysis on clinical data 
#      from cancer patients using Kaplan-Meier survival curves 
#      and Cox proportional hazards models. It generates forest plots
#      to visualize the effect of clinical covariates and subtypes 
#      on survival outcomes. 
#
# Requirements:
#      - Configuration file (`config.tsv`) with directory paths.
#      - Clinical data file (`clinical_data.txt`).
###################################################
# Required environment: survival


############### Library Setup ###############
# Load necessary libraries for survival analysis, plotting, and model fitting


library(survival)
library(survminer)
library(forestplot)
library(tidyverse)
library(aod)
library(forestmodel)

# Load custom survival functions
source("functions/clustering/survival.functions.R")

############### Configuration File Loading ###############
# Read the configuration file (config.tsv) to obtain directory paths and other settings

config <- read.csv("config/config.tsv",
sep = "\t",
header = T)

# Extract key-value pairs from the config file for directory setup
keys <- config$key
values <- config$value
names(values) <- keys
wkdir <- values["wkdir"]
datadir <- values["datadir"]


############### Directory Setup ###############
# Check if the survival analysis output directory exists, if not create it

if (dir.exists(paste0(wkdir, "/clustering/survival")) == F) {

    dir.create(paste0(wkdir, "/clustering/survival"))

}

# Set the output directory for survival analysis results
outdir <- paste0(wkdir, "/clustering/survival")



############### Data Loading ###############
# Load clinical data from the specified file (clinical_data.txt)
base <- read.csv(paste0(datadir, "/clinical_data.txt"),
                 sep = "\t",
                 row.names = 1,
                 header = T)

# Set the factor levels for the primary tumor column
base$PRIMARY_TUMOR <- factor(base$PRIMARY_TUMOR,
levels = c("GASTRIC","COLORECTAL", "LUNG", "PANCREAS", "SMALL INTESTINE"))
############################################


############################################
############### Forest Plot Generation ###############
# Prepare data for the forest plot by selecting relevant columns (time, status, and subtype)

df <- base %>%
    transmute(
        time = OS.time,
        status = EXITUS,
        SUBTYPE = Subtype
    )

# Create a univariate forest plot of the Cox proportional hazards model for survival

p <- forest_model(coxph(Surv(time, status) ~ ., df))


# Save the forest plot as a PDF file
ggsave(p,
file =paste0(outdir, "/Forest_plot_Univariate_.pdf"),
device = "pdf", width =8 , height = 5)


############### Kaplan-Meier Curve Generation ###############
# Create Kaplan-Meier curves for the subtypes, cutting survival time at 150 months

Kapplan_Meyer(cut_survival(base, 150),
class = "Subtype",
OS = "OS.time",
OS_event = "EXITUS",
outdir = outdir,
flag = "cut_150",
custcolor = c("gold", "cornflowerblue", "tomato1")
)



############### Univariate Cox Model Analysis ###############
# Run univariate Cox models for multiple clinical features related to OS

for (feature in c("Subtype",
"GRADE",
"STAGE_IV",
"AGE",
"PRIMARY_TUMOR",
"GENDER")) {

Hazard_Cox(base,
class = feature,
OS = "OS.time",
OS_event = "EXITUS",
covariates = NULL,
outdir = outdir,
flag = paste0("Univariate_",feature))
}


############### Multivariate Cox Model Analysis ###############
# Run a multivariate Cox model with additional covariates
Hazard_Cox(base,
class = "Subtype",
OS = "OS.time",
OS_event = "EXITUS",
covariates = c("GRADE", "STAGE_IV", "AGE", "PRIMARY_TUMOR", "GENDER"),
outdir = outdir,
flag = "Multivariate_stage_IV")


############### Additional Forest Plot Generation ###############
# Create forest plots adjusting for multiple covariates

df <- base %>%
    transmute(
        time = OS.time,
        status = EXITUS,
        GENDER,
        GRADE,
        AGE,
        STAGE.IV = STAGE_IV,
        PRIMARY.TUMOR = PRIMARY_TUMOR,
        SUBTYPE = Subtype
    )

# Create the multivariate forest plot
p <- forest_model(coxph(Surv(time, status) ~ ., df))

# Save the forest plot as a PDF file
ggsave(p,
file =paste0(outdir, "/Forest_plot_covariates_.pdf"),
device = "pdf", width =12 , height = 10)

############### Stratified Kaplan-Meier Curves by Stage and Grade ###############
# Generate Kaplan-Meier curves and forest plots for stage IV and other stages

for (meta in unique(base$STAGE_IV)){
    if (is.na(meta) == F){
    df <- base[which(base$STAGE_IV == meta),]

        # Kaplan-Meier curve for the specific stage
    Kapplan_Meyer(cut_survival(df, 150),
    class = "Subtype",
    OS = "OS.time",
    OS_event = "EXITUS",
    outdir = outdir,
    flag = paste0("cut_150_Metastasis", meta),
    custcolor = c("gold", "cornflowerblue", "tomato1")
    )
    
    
         # Generate forest plot for this stage   
df <- df %>%
    transmute(
        time = OS.time,
        status = EXITUS,
        GENDER,
        GRADE,
        AGE,
        PRIMARY.TUMOR = PRIMARY_TUMOR,
        SUBTYPE = Subtype
    )

p <- forest_model(coxph(Surv(time, status) ~ ., df))

# Save the plot
ggsave(p,
file =paste0(outdir, "/Forest_plot_all_",meta,".pdf"),
device = "pdf", width =12 , height = 10)
        
    }

}

# Generate Kaplan-Meier curves and forest plots for each tumor grade
for (meta in unique(base$GRADE)){

    df <- base[base$GRADE == meta,]
    # Kaplan-Meier curve for the specific grade
    Kapplan_Meyer(cut_survival(df, 150),
    class = "Subtype",
    OS = "OS.time",
    OS_event = "EXITUS",
    outdir = outdir,
    flag = paste0("cut_150_Grade", meta),
    custcolor = c("gold", "cornflowerblue", "tomato1")
    )

    # Generate forest plot for this grade
df <- df %>%
    transmute(
        time = OS.time,
        status = EXITUS,
        GENDER,
        STAGE.IV = STAGE_IV,
        AGE,
        PRIMARY.TUMOR = PRIMARY_TUMOR,
        SUBTYPE = Subtype
    )

p <- forest_model(coxph(Surv(time, status) ~ ., df))

    # Save the plot
ggsave(p,
file =paste0(outdir, "/Forest_plot_all_",meta,".pdf"),
device = "pdf", width =12 , height = 10)



}

###############################################################
# Run multivariate Cox model without GRADE as covariate

df <- base %>%
    transmute(
        time = OS.time,
        status = EXITUS,
        GENDER,
        AGE,
        STAGE.IV = STAGE_IV,
        PRIMARY.TUMOR = PRIMARY_TUMOR,
        SUBTYPE = Subtype
    )

p <- forest_model(coxph(Surv(time, status) ~ ., df))

ggsave(p,
file =paste0(outdir, "/Forest_plot_covariates_not_GRADE_.pdf"),
device = "pdf", width =12 , height = 10)



Hazard_Cox(base,
class = "Subtype",
OS = "OS.time",
OS_event = "EXITUS",
covariates = c("STAGE_IV", "AGE", "PRIMARY_TUMOR", "GENDER"),
outdir = outdir,
flag = "Multivariate_not_Grade")



############### Kaplan-Meier by Tissue Type ###############
# Generate Kaplan-Meier curves and Cox models by tissue type

for (tissue in unique(base$PRIMARY_TUMOR)){

    df <- base[base$PRIMARY_TUMOR == tissue,]
    # Kaplan-Meier curve for the specific tissue
    Kapplan_Meyer(cut_survival(df, 150),
    class = "Subtype",
    OS = "OS.time",
    OS_event = "EXITUS",
    outdir = outdir,
    flag = paste0("cut_150_", tissue),
    custcolor = c("gold", "cornflowerblue", "tomato1")
    )

    # Generate Cox model for this tissue
Hazard_Cox(df,
class = "Subtype",
OS = "OS.time",
OS_event = "EXITUS",
covariates = c("GRADE", "STAGE_IV", "AGE", "GENDER"),
outdir = outdir,
flag = paste0("Multivariate_stage_IV_",tissue))


df <- df %>%
    transmute(
        time = OS.time,
        status = EXITUS,
        GENDER,
        GRADE,
        AGE,
        STAGE.IV = STAGE_IV,
        PRIMARY.TUMOR = PRIMARY_TUMOR,
        SUBTYPE = Subtype
    )

p <- forest_model(coxph(Surv(time, status) ~ ., df))

ggsave(p,
file =paste0(outdir, "/Forest_plot_all_stage_IV_",tissue,".pdf"),
device = "pdf", width =12 , height = 10)


}

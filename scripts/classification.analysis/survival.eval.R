###################################################
# Script: Survival Analysis with Cox Models
#
# Description:
#      This script evaluates the performance of Cox proportional hazards models 
#      for predicting overall survival (OS) in cancer patients. The models are 
#      built with and without including the "Subtype" variable, and the C-index 
#      of both models is compared to assess the added value of the subtype 
#      information. The script also generates a box plot comparing the C-indices 
#      and performs statistical testing using paired t-tests.
#
# Requirements:
#      - Clinical data file (`clinical_data.txt`) containing patient survival 
#        information and covariates like tumor type, stage, and gender.
###################################################
# Required environment: surveval

############### Library Setup ###############
# Load necessary libraries for survival analysis, plotting, and statistical testing

library(survival)  
library(ggplot2)
library(ggpubr)
library(rstatix)
library(tidyverse)

# Set a random seed for reproducibility
set.seed(123)

############### Configuration File Loading ###############
# Read the configuration file (config.tsv) to obtain directory paths and other settings

config <- read.csv("config/config.tsv",
sep = "\t",
header = T)

# Extract key-value pairs from the config file for directory setup
keys <- config$key
values <- config$value
names(values) <- keys

# Set working directory and data directory paths based on the config file

wkdir <- values["wkdir"]
datadir <- values["datadir"]

############### Directory Setup ###############
# Check if the directory for survival evaluation results exists, if not create it

if (dir.exists(paste0(wkdir, "/clustering/survival.eval")) == F) {

    dir.create(paste0(wkdir, "/clustering/survival.eval"))

}

# Set the output directory for survival evaluation results
outdir <- paste0(wkdir, "/clustering/survival.eval")



############### Data Loading and Preparation ###############
# Load clinical data from the specified file (clinical_data.txt)

base <- read.csv(paste0(datadir, "/clinical_data.txt"),
                 sep = "\t",
                 row.names = 1,
                 header = T)


# Set PRIMARY_TUMOR as a factor with specific levels (types of tumors)
base$PRIMARY_TUMOR <- factor(base$PRIMARY_TUMOR,
levels = c("GASTRIC","COLORECTAL", "LUNG", "PANCREAS", "SMALL INTESTINE"))


# Ensure survival time and event status (OS.time and EXITUS) are numeric  
base$OS.time <- as.numeric(base$OS.time)
base$OS.event <- as.numeric(base$EXITUS)

# Generate the survival object (Surv) using overall survival time (OS.time) and event status (OS.event)
base$Surv_OS <- with(base, Surv(base$OS.time, base$OS.event == 1))

############### Cox Proportional Hazards Models and C-index Calculation ###############
# Define the features to include in the Cox models
features <- c("GRADE", "AGE", "PRIMARY_TUMOR", "GENDER","STAGE_IV")

# Initialize empty lists to store C-index values for the models
cindexm1 <- list()
cindexm2 <- list()

# Run the Cox models 100 times with random sampling
for (i in 1:100){

    # Randomly sample 80% of the data
index_data = sample(1:nrow(base), nrow(base)*0.8)

# Fit Cox model without Subtype (using only the predefined features)
m1 = coxph( as.formula(paste('Surv_OS ~ ',paste(features, collapse = ' + '))), data = base[index_data,], x = T)

# Fit Cox model with Subtype added to the predefined features
m2 = coxph( as.formula(paste('Surv_OS ~ Subtype +',paste(features, collapse = ' + '))), data = base[index_data,], x = T)

# Obtain the concordance index (C-index) for both models
cindex <- concordance(m1,m2)

# Store the C-index for each model in the lists
cindexm1 <- append(cindexm1,cindex[[1]][1])
cindexm2 <- append(cindexm2,cindex[[1]][2])
}


############### C-index Comparison and Statistical Test ###############
# Create a dataframe containing the C-index values for both models
df <- data.frame("Cox.model" = unlist(cindexm1), "Cox.model.NS" = unlist(cindexm2))

# Reshape the data to long format for easier plotting and statistical analysis
df <- df %>% pivot_longer(c(Cox.model, Cox.model.NS),names_to = "Model", values_to = "C.index")

# Perform a paired t-test to compare the C-index values between the two models
stat.test <- df  %>% 
  t_test(C.index ~ Model, paired = T) %>%
  add_significance() %>% add_xy_position(x = "Model")

# Clear the "Model" column for plotting (no need to display the model names in test results)
stat.test[["Model"]] <- NA


############### C-index Comparison Box Plot ###############
# Create a box plot comparing the C-index values for the two models

# Add a box plot with jittered points to show individual C-index values

plt <- ggplot( df , aes(y = .data[["C.index"]],
    x = .data[["Model"]], fill = .data[["Model"]])) 
        

        plt <- plt + geom_boxplot(
            width = 0.6, outlier.shape = NA
        ) + 
        geom_point(position = position_jitterdodge()) 

# Customize the plot with colors, theme, and labels         
plt <- plt + scale_fill_manual(values = list("Cox.model" = "darkslategray3",
    "Cox.model.NS" = "slateblue3")) +

        theme_classic() +
        stat_pvalue_manual(stat.test,
        label = "p.signif",  hide.ns = T,
        step.increase = 0.1)  +

        theme(text = element_text(size = 40, colour = "black" ),
        axis.line = element_line(size = 1, colour = "black"),
        axis.text.y = element_text(margin = margin(20, 10, 20, 20),
        colour = "black"),
        axis.text.x = element_text(colour = "black"),
        legend.key.size = unit(3,"cm"),
        axis.ticks.length.y = unit(10,"pt"),
    plot.title = element_text(hjust = 0.5, vjust= 4)
) +
        xlab("") 

# Customize the size of significance labels
    plt$layers[[3]]$aes_params$size <- 1
    plt$layers[[3]]$aes_params$label.size <- 20
    plt$layers[[3]]$aes_params$vjust <- 0.5

# Save the plot as a PDF
    ggsave(
        plot = plt, filename = paste0(
            outdir, "/Cindex.pdf"
        ),
        device = "pdf",
        dpi = 500, width = 20,
        height = 15
    )


############### Save the Statistical Test Results ###############
# Write the t-test results to a text file

write.table(stat.test[,c(1:9)],
paste0(outdir,"/T_test_results,txt"),
sep = "\t",
row.names = F,
col.names = T)

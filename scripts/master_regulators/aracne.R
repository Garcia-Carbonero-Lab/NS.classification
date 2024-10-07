# This script processes gene expression data to identify master regulators using the ARACNe algorithm.
# The script prepares expression data, runs ARACNe for network inference using bootstrapping, consolidates the results, 
# and generates a consolidated regulatory network for further analysis.
#
# Requirements:
# - ARACNe-AP must be installed and available in the working directory.
# - A quality-adjusted expression matrix must be available.
# - Master regulator gene lists in HUGO format must be present.
# - Configurations for directories are specified in "config/config.tsv".
#
# Inputs:
# - expression.qual.txt: Quality-adjusted gene expression matrix.
# - Master regulator (MR) HUGO gene list(s) in master_regulators directory.
#
# Outputs:
# - Network files for individual Master regulator genes.
# - A consolidated network file combining results from all Master regulator gene lists.
#
# Create directories if necessary, install ARACNe, run ARACNe analysis, and consolidate the results.

# Required environment: 'viper.yml' (containing prerequiste libraries).



####################################################
# Load configuration and set directories
####################################################
# Open config file
config <- read.csv("config/config.tsv",
  sep = "\t",
  header = T
)

keys <- config$key
values <- config$value
names(values) <- keys
# Set working and data directories from the config file
wkdir <- values["wkdir"]
datadir <- values["datadir"]


####################################################
# Create Directories for Results
####################################################

# Create directories for storing master regulator results, if not already present
if (dir.exists(paste0(wkdir, "/master_regulators")) == F) {
  dir.create(paste0(wkdir, "/master_regulators"))
  dir.create(paste0(wkdir, "/master_regulators/aracne"))
}


####################################################
# Install ARACNe (commented out if already installed)
####################################################

setwd(paste0(datadir, "/master_regulators"))


# Clone ARACNe-AP repository and compile it (uncomment if running for the first time)
# system("git clone https://github.com/califano-lab/ARACNe-AP.git")


setwd("ARACNe-AP")

#system("ant main")


####################################################
# Prepare Expression Data for ARACNe
####################################################

# Load quality-adjusted gene expression matrix
expression <- read.table(paste0(wkdir,
"/preprocess/transcriptome/expression.qual.txt"),
  sep = "\t",
  header = T,
  check.names = F
)


# Prepare the expression matrix by adding a "gene" column as required by ARACNe
colnames(expression) <- c("gene", colnames(expression[2:ncol(expression)]))

# Write the formatted expression matrix to a file for ARACNe input
write.table(expression,
  paste0(datadir, "/master_regulators/expression_aracne.tsv"),
  sep = "\t",
  row.names = F,
  col.names = T,
  quote = F
)


####################################################
# Run ARACNe for Each Master Regulator Gene List
####################################################

# Identify all MR gene lists in the master_regulators directory (files containing 'hugo' in the name)
reg <- list.files((paste0(datadir, "/master_regulators")))
reg <- reg[grep("hugo", reg)]


# Loop over each MR list to calculate networks
for (r in reg) {
  name <- gsub("\\.txt", "", r)
  # Step 1: Calculate threshold for mutual information
   
    system(paste0(
    "java -Xmx5G ",
    "-jar dist/aracne.jar ",
    "-e ", datadir, "/master_regulators/expression_aracne.tsv ",
    "-o ", wkdir, "/master_regulators/aracne/", name,
    " --tfs ", datadir, "/master_regulators/", r,
    " --pvalue 1E-8 --seed 1 --calculateThreshold"
  ))

  # Step 2: Perform bootstrapping to infer networks (100 bootstrap iterations)
  
    system(paste0(
    "for i in $(seq 1 100) ; do ",
    "java -Xmx5G ",
    "-jar dist/aracne.jar ",
    "-e ", datadir, "/master_regulators/expression_aracne.tsv ",
    "-o ", wkdir, "/master_regulators/aracne/", name,
    " --tfs ", datadir, "/master_regulators/", r,
    " --pvalue 1E-8 --seed $i; done"
  ))
  
  # Step 3: Consolidate the bootstrapped networks into a single network

  system(paste0(
    "java -Xmx10G -jar dist/aracne.jar ",
    "-o ", wkdir, "/master_regulators/aracne/", name,
    " --consolidate"
  ))


  ####################################################
  # Save ARACNe Result as TSV File
  ####################################################

  # Load the consolidated network and remove rows with missing values
  aracne <- read.table(paste0(wkdir, "/master_regulators/aracne/",
  name, "/network.txt")
  , sep = "\t"
  , header = T)

  aracne <- aracne[complete.cases(aracne[, 3]),]
  
  # Write the cleaned network to a .tsv file
  write.table(aracne, paste0(wkdir, "/master_regulators/aracne/",
  name, "/network.tsv"),
  sep = "\t",
  row.names = F,
  col.names = F,
  quote = F)

}


####################################################
# Concatenate All Networks
####################################################

# Concatenate all individual network files into one combined network file

system(paste0("cat ", wkdir, "/master_regulators/aracne/*/network.tsv > ",
wkdir, "/master_regulators/aracne/concatenated_networks.tsv"))

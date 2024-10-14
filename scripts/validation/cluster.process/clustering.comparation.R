###############################################################################
# Script: Rename Subtypes in Validation Dataset
#
# Description:
#     This script renames subtypes in the validation dataset based on the most
#     frequent subtype matches found in the combination dataset. The resulting
#     validation dataset is saved with updated subtype names.
#
# Inputs:
#     - Combination dataset: Base_comparation.subtype.tsv
#     - Validation dataset: Base_validation.subtype.tsv
#     - Configuration file: config.tsv
#
# Outputs:
#     - Updated validation dataset with renamed subtypes
###############################################################################
# Load configuration file to get working and data directories
config <- read.csv("config/config.tsv",
sep = "\t",
header = T)
# Extract key-value pairs for configuration settings
keys <- config$key
values <- config$value
names(values) <- keys
# Define working directory and data directory
wkdir <- values["wkdir"]
datadir <- values["datadir"]

# Load combination and validation datasets
base_combination <- read.csv(paste0(wkdir, "/validation/base/Base_comparation.subtype.tsv"),
                 sep = "\t",
                 row.names = 1,
                 header = T)


base.validation <- read.table(paste0(wkdir, "/validation/base/Base_validation.subtype.tsv"),
sep = '\t',
header = T,
row.names = 1
)
##################################################################
# Rename subtypes in the combination dataset
# Create a frequency table comparing Subtype.combination with Subtype
df.comp <- table(base_combination$Subtype.combination,base_combination$Subtype)

# Get the most common subtype matches for each combination subtype
levels.comb <- apply(df.comp, 2, function(x) {
names(x)[x == max(x)]
})

# Convert Subtype.comb to character for renaming
base_combination$Subtype.comb <- as.character(base_combination$Subtype.comb)

# Rename subtypes in the combination dataset based on frequency matches
for (i in 1:length(levels.comb)){
base_combination$Subtype.comb <- gsub(paste0("^",levels.comb[i],"$")
,names(levels.comb)[i],
base_combination$Subtype.comb)
}
# Create a frequency table comparing the renamed subtypes with validation subtypes
df.val <- table(base_combination$Subtype.comb,base_combination$Subtype.validation)

# Get the most common validation subtype matches
levels.val <- apply(df.val, 2, function(x) {
names(x)[x == max(x)]
})

# Rename subtypes in the validation dataset based on frequency matches
for (i in 1:length(levels.val)){
base.validation$Subtype.validation <- gsub(paste0("^", names(levels.val[i]),"$")
,levels.val[i],
base.validation$Subtype.validation)

}
# Save the updated validation dataset
write.table(base.validation,
paste0(wkdir, "/validation/base/Base_validation.subtype.tsv"),
sep = "\t",
row.names = T,
col.names = T)

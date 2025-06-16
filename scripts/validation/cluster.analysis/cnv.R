#clinical
#environment is figuras
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(reshape2)

source("functions/main.functions.R")
source("functions/cluster.analysis/additional.functions.R")
source("functions/cluster.analysis/clinical.functions.R")

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



if (dir.exists(paste0(wkdir, "/validation/clustering.analysis/cnv")) == F) {

    dir.create(paste0(wkdir, "/validation/clustering.analysis/cnv"))

}


outdir <- paste0(wkdir, "/validation/clustering.analysis/cnv")




base <- read.csv(paste0(wkdir, "/validation/clustering.rgc.all/Base_validation.subtype.tsv"),
                 sep = "\t",
                 row.names = 1,
                 header = T)

#study differences in CNV among subtypes
table <- clinical.table(data = base, 
features = c("chr18del", "chr4gain", "chr5gain", "chr7gain", "chr10gain", "chr14gain", "chr20gain"),
group = "Subtype.validation")
write.table(table, paste0(outdir, "/CNV.txt"),
sep = "\t", row.names = T, col.names = T, quote = F)


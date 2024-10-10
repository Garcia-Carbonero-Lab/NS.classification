#this script needs met.arrays environment

library(minfi)
library(ENmix)
library(wateRmelon)
library(MASS)
library(broom)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)



if (dir.exists(paste0(wkdir, "/validation")) == F) {

    dir.create(paste0(wkdir, "/validation"))

}

if (dir.exists(paste0(wkdir, "/validation/Laddha")) == F) {

    dir.create(paste0(wkdir, "/validation/Laddha"))

}

outdir <- paste0(wkdir, "/validation/Laddha")


config <- read.csv("config/config.tsv",
sep = "\t",
header = T)

keys <- config$key
values <- config$value
names(values) <- keys
###############
wkdir <- values["wkdir"]
datadir <- values["datadir"]
source("functions/validation/validation.functions.R")

path <- paste0(datadir, "/external_data/Laddha_LUNET/IDAT")

generate.sample.sheet(path)

methylation_matrix(paste0(datadir, "/external_data/Laddha_LUNET"),
arraytype = "450K",
outdir)

# This script needs preprocess environment
library(tidyverse)

config <- read.csv("config/config.tsv",
    sep = "\t",
    header = T
)

keys <- config$key
values <- config$value
names(values) <- keys
###############
wkdir <- values["wkdir"]
datadir <- values["datadir"]

if (dir.exists(paste0(wkdir, "/validation/whole.matrix")) == F) {

    dir.create(paste0(wkdir, "/validation/whole.matrix"))

}

outdir <- paste0(wkdir, "/validation/whole.matrix")

#Open base

base <- read.table(paste0(wkdir, "/validation/base/Base_validation.tsv"),
sep = "\t",
row.names = 1,
header = T)

#Alcala & Laddha
#expression data from Laddha is in the matrix of Alcala
exp_Alcala <- read.table(paste0(wkdir,"/validation/Alcala/Alcala.expression.tsv"),
sep = "\t",
row.names = 1,
header = T)

#Select samples in base with both omics
exp_Alcala <- exp_Alcala[,colnames(exp_Alcala) %in% rownames(base)]

#Lakis
exp_Lakis<- read.table(paste0(wkdir,"/validation/Lakis/matrix_Lakis_expression.txt"),
sep = "\t",
row.names = 1,
header = T)
# log to Lakis
#check colnames in base
exp_Lakis <- exp_Lakis[,colnames(exp_Lakis) %in% rownames(base)]
#################################
#Karpathakis
exp_Karpathakis<- read.table(paste0(wkdir,"/validation/Karpathakis/Karpathakis.expression.tsv"),
sep = "\t",
row.names = 1,
header = T,
check.names = F)
exp_Karpathakis <- log(exp_Karpathakis +1)
exp_Karpathakis <- exp_Karpathakis[,colnames(exp_Karpathakis) %in% rownames(base)]
############################################################################################
#Select expression fron datasets with both omics
#Lakis; Alcala; Ladha ;  Karpathakis
combination <- merge(exp_Lakis, exp_Alcala,by.x = "row.names",by.y = "row.names")
rownames(combination) <- combination[,1]
combination <- combination[,-1]
combination <- merge(combination, exp_Karpathakis,by.x = "row.names",by.y = "row.names")
rownames(combination) <- combination[,1]
combination <- combination[,-1]
write.table(combination, paste0(outdir,"/expression.combination.all.tsv"),
sep = "\t",
row.names = T,
col.names= NA)

###########################################################################################
#selection of genes used in MOFA model


exp <- read.table(paste0(wkdir,
"/preprocess/transcriptome/expression_fs.txt"),
    sep = "\t",
    header = T,
    row.names = 1,
    check.names = F
)

#filter genes used in MOFA

combination <- combination[rownames(exp),]




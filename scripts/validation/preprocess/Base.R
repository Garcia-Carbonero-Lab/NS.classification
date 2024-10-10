# this  script uses geo environment
library(tidyverse)
library(GEOquery)

config <- read.csv("config/config.tsv",
sep = "\t",
header = T)

keys <- config$key
values <- config$value
names(values) <- keys
wkdir <- values["wkdir"]
datadir <- values["datadir"]


if (dir.exists(paste0(wkdir, "/validation")) == F) {

    dir.create(paste0(wkdir, "/validation"))

}

if (dir.exists(paste0(wkdir, "/validation/base")) == F) {

    dir.create(paste0(wkdir, "/validation/base"))

}

outdir <- paste0(wkdir, "/validation/base")


##################################################################
#Lakis

base_Lakis <- read.table(paste0(wkdir, "/validation/Lakis/base_Lakis.tsv"),
sep = "\t", row.names = 1, header = T)

rownames(base_Lakis) <- base_Lakis[,1]

base_Lakis <- base_Lakis[,-1]
head(base_Lakis)

base_Lakis$Gender <- as.factor(base_Lakis$Gender)
levels(base_Lakis$Gender) <- c("Female", "Male")

base_Lakis$tissue <- as.factor(base_Lakis$primary.tumor)
levels(base_Lakis$tissue) <- c("Primary.tumor", "Primary.tumor")


base_Lakis$primary.tumor <- as.factor(base_Lakis$primary.tumor)
levels(base_Lakis$primary.tumor) <- c("pancreas", "pancreas")

base_Lakis$Differenciation <- as.factor(base_Lakis$Differenciation)
levels(base_Lakis$Differenciation) <- c("NET", "NET", "NET","NET")

base_Lakis$Grade <- as.factor(base_Lakis$Grade)
levels(base_Lakis$Grade) <- c(NA, "G1", "G2")

base_Lakis$Stage <- as.factor(base_Lakis$Stage)
levels(base_Lakis$Stage) <- c(NA, "I", "I", "I","II","II","III","III",
"III","II","II","II","II", "IV","III","III","III")

base_Lakis$Overall.survival <- base_Lakis$Overall.survival / 30.41

base_Lakis$Status <- ifelse(base_Lakis$Status == "alive",0,1)

base_Lakis$paper <- rep("Lakis_PANET", nrow(base_Lakis))

columns_prev <- colnames(base)[!colnames(base) %in% colnames(base_Lakis)]

base <- base_Lakis
#############################################################################
# Alcala

base_Alcala <- read.table(paste0(datadir, "/external_data/Alcala_LUNET/Attributes.txt"),
sep = "\t", row.names = 1, header = T)

columns <- c("Age", "Sex", "Smoking_status", "Survival_months", "Survival_censor",
"Histopathology", "Stage_UICC", "cluster_LNET")

base_Alcala <- base_Alcala[,columns]

colnames(base_Alcala) <- c("Age", "Gender", "Smoking",
"Overall.survival","Status", "Differenciation", "Stage","Clusters")

base_Alcala$Gender <- as.factor(base_Alcala$Gender)
levels(base_Alcala$Gender) <- c("Female", "Male")

base_Alcala$Grade <- as.factor(base_Alcala$Differenciation)
levels(base_Alcala$Grade) <- c("G2",NA,"G3","G3", "G3", "G3","G3","G1")

base_Alcala$Differenciation <- as.factor(base_Alcala$Differenciation)
levels(base_Alcala$Differenciation) <- c("NET.atypical", "NET", "Carcinoma","Carcinoma",
"Carcinoma", "Carcinoma","Carcinoma","NET.typical")

base_Alcala <- base_Alcala[base_Alcala$Differenciation != "Carcinoma",]

base_Alcala$primary.tumor <- rep("lung", nrow(base_Alcala))

base_Alcala$tissue <- rep("Primary.tumor", nrow(base_Alcala))

base_Alcala$paper <- rep("Alcala_LUNET", nrow(base_Alcala))

base_Alcala$Status <- as.factor(base_Alcala$Status)
levels(base_Alcala$Status) <- c(0,1)

base_Alcala$Stage <- as.factor(base_Alcala$Stage)
levels(base_Alcala$Stage) <- c("I","I","I","I","I","I","II","II","III","III","IV","IV")

base_Alcala$Clusters <- as.factor(base_Alcala$Clusters)
levels(base_Alcala$Clusters) <- c("Carcinoid-A1", "Carcinoid-A2","Carcinoid-B")

columns_prev <- colnames(base)[!colnames(base) %in% colnames(base_Alcala)]

for (column in columns_prev){

base_Alcala[,column] <- rep(NA, nrow(base_Alcala))


}

columns_next <- colnames(base_Alcala)[!colnames(base_Alcala) %in% colnames(base)]

for (column in columns_next){

base[,column] <- rep(NA, nrow(base))


}

base_Alcala <- base_Alcala[,colnames(base)]

base <- rbind(base,base_Alcala)
############################################################################################
#Laddha
#SAMPLES IN LADDHA ARE IN ALCALA

GSE118131 <- getGEO("GSE118131", GSEMatrix =TRUE)
GSE118131 <- GSE118131[[1]]
base_Laddha <- pData(phenoData(GSE118131))

srr <- c("SRR7646229" ,"SRR7646230" ,"SRR7646231" ,"SRR7646232"
 ,"SRR7646233" ,"SRR7646234" ,"SRR7646235" ,"SRR7646236" ,"SRR7646237"
  ,"SRR7646238" ,"SRR7646239" ,"SRR7646240" ,"SRR7646241" ,"SRR7646242" 
  ,"SRR7646243" ,"SRR7646244" ,"SRR7646245" ,"SRR7646246" ,"SRR7646247" 
  ,"SRR7646248" ,"SRR7646249" ,"SRR7646250" ,"SRR7646251" ,"SRR7646252" ,"SRR7646253"
 ,"SRR7646254" ,"SRR7646255" ,"SRR7646256" ,"SRR7646257" ,"SRR7646258")

base_Laddha$srr <- srr 

columns <- c("title","srr","age (years):ch1", "cell type:ch1", "diagnosis:ch1"
, "gender:ch1", "lcs subtype:ch1", "phenotype:ch1", "tumor location:ch1")


base_Laddha <- base_Laddha[,columns]

rownames(base_Laddha) <- base_Laddha$srr

base_Laddha <- base_Laddha[,3:ncol(base_Laddha)]

colnames(base_Laddha) <- c("Age", "tissue", "Differenciation", "Gender", "Clusters", "primary.tumor", "tumor.location") 

base_Laddha$tissue <- rep("Primary.tumor",nrow(base_Laddha))
base_Laddha$Differenciation <- as.factor(base_Laddha$Differenciation)
levels(base_Laddha$Differenciation) <- c("NET.atypical", "NET.typical")

base_Laddha$Gender <- as.factor(base_Laddha$Gender)
levels(base_Laddha$Gender) <- c("Female", "Male")


base_Laddha$Clusters <- as.factor(base_Laddha$Clusters)
levels(base_Laddha$Clusters) <- c("LC1", "LC2", "LC3")

base_Laddha$primary.tumor <- rep("lung", nrow(base_Laddha))
base_Laddha$paper <- rep("Laddha_LUNET",nrow(base_Laddha))

columns_prev <- colnames(base)[!colnames(base) %in% colnames(base_Laddha)]

for (column in columns_prev){

base_Laddha[,column] <- rep(NA, nrow(base_Laddha))


}

columns_next <- colnames(base_Laddha)[!colnames(base_Laddha) %in% colnames(base)]

for (column in columns_next){


base[,column] <- rep(NA, nrow(base))


}

base_Laddha <- base_Laddha[,colnames(base)]

base <- rbind(base,base_Laddha)
######################################################################################
# Karpathakis

base_Karpathakis <- read.table(paste0(datadir,"/external_data/Karpathakis_SINET/Base_Karpathakis.tsv"),
sep = "\t", row.names = 1, header = T)

base_Karpathakis <- base_Karpathakis[which(base_Karpathakis$Gender != "cell line" | is.na(base_Karpathakis$Gender)),]

base_Karpathakis$Gender <- as.factor(base_Karpathakis$Gender)
levels(base_Karpathakis$Gender) <- c("Female", "Male")

base_Karpathakis$Grade <- as.factor(base_Karpathakis$Grade)
levels(base_Karpathakis$Grade) <- c("G1", "G2", "G3", NA, NA)

base_Karpathakis$primary.tumor <- as.factor(base_Karpathakis$primary.tumor)
levels(base_Karpathakis$primary.tumor) <- c("appendix", "colorrectal","stomach","small.intestine")

base_Karpathakis <- base_Karpathakis[which(base_Karpathakis$primary.tumor == "small.intestine"),]

#N are normals

base_Karpathakis$tissue[grep("_N",base_Karpathakis$Patient)] <- "Normal"
base_Karpathakis$tissue <- as.factor(base_Karpathakis$tissue)
levels(base_Karpathakis$tissue) <- c("Liver", "Mesenteric","Lymph.node", "Primary.tumor", "Mesenteric")

base_Karpathakis$paper <- rep("Karpathakis_SINET", nrow(base_Karpathakis))

columns_prev <- colnames(base)[!colnames(base) %in% colnames(base_Karpathakis)]

for (column in columns_prev){

base_Karpathakis[,column] <- rep(NA, nrow(base_Karpathakis))


}

columns_next <- colnames(base_Karpathakis)[!colnames(base_Karpathakis) %in% colnames(base)]

for (column in columns_next){

base[,column] <- rep(NA, nrow(base))


}

base_Karpathakis <- base_Karpathakis[,colnames(base)]

base <- rbind(base,base_Karpathakis)

base <- base[complete.cases(base[,1]),]

rownames(base)[base$paper == "Laddha_LUNET"] <- rownames(base_Laddha)


write.table(base, paste0(outdir, "/Base_validation.tsv"),
sep = "\t",
row.names = T,
col.names = NA)

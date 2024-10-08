#this script needs geo env
library(tidyverse)

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

if (dir.exists(paste0(wkdir, "/validation/Lakis")) == F) {

    dir.create(paste0(wkdir, "/validation/Lakis"))

}

outdir <- paste0(wkdir, "/validation/Lakis")

base_icgc_donor <- read.table(paste0(datadir, "/external_data/Lakis_PANET/clinical_data_donor.tsv"),
sep = "\t",
header= T)


base_icgc_specimen <- read.table(paste0(datadir, "/external_data/Lakis_PANET/specimen.PAEN-AU.tsv"),
sep = "\t",
header= T)

base_icgc_sample<- read.table(paste0(datadir, "/external_data/Lakis_PANET/clinical.data.tsv"),
sep = "\t",
header= T)


base <- merge(base_icgc_donor,base_icgc_specimen,by.x = "icgc_donor_id", by.y = "icgc_donor_id")
base <- merge(base,base_icgc_sample,by.x = "icgc_specimen_id", by.y = "icgc_specimen_id")



columns <- c("icgc_sample_id","icgc_donor_id.x", "donor_sex", 
"donor_vital_status", "donor_age_at_diagnosis", "donor_tumour_stage_at_diagnosis",
"donor_survival_time", "percentage_cellularity.x", "disease_status_last_followup",
"specimen_type", "tumour_histological_type", "tumour_grade", "tumour_stage")

base <- base[,columns]

colnames(base) <- c("id_sample", "id_patient", "Gender", "Status", "Age", "Stage"
, "Overall.survival", "Celularity", "dissease", "primary.tumor",
"Differenciation","Grade", "Stage")

write.table(base,
paste0(outdir, "/base_Lakis_prev.tsv"),
sep = "\t",
row.names = F,
col.names = T)


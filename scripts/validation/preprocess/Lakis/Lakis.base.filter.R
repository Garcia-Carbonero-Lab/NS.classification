base <- read.table(paste0(outdir, "/base_Lakis_prev.tsv"),
sep = "\t",
header = T)

########################################################################
#Select base there is in data
expression <- read.table(paste0(outdir, "/matrix_Lakis_expression.txt"),
sep = "\t",
row.names = 1,
header = T)

methylation <- read.table(paste0(outdir, "/Methylation_M_Lakis.tsv"),
sep = "\t",
row.names = 1,
header = T)

patients <- unique(c(colnames(expression), colnames(methylation)))
base <- base[rownames(base)[base$id_patient %in% patients],]

base <- base[order(base$id_patient),]
base <- base[!base$primary.tumor %in% c(
    "Normal - tissue adjacent to primary", "Normal - blood derived",
    "Normal - solid tissue"

),]
base <- base[order(base$Celularity),]


base <- base[!duplicated(base$id_patient),]
############################################################

###############################################################33
#mutations

mutations <- read.table(paste0(datadir, "/external_data/Lakis_PANET/simple_somatic_mutation.open.PAEN-AU.tsv"),
sep = "\t",
header = T)

mutations <- mutations[!duplicated(mutations$icgc_mutation_id),]
genes <- unique((mutations$gene_affected))
genesID <-mygene::queryMany(genes, scopes="ensembl.gene", fields="symbol", species="human")
genesID <- as.data.frame(genesID[,c("query","symbol")])

mutations <- merge(mutations,genesID, by.x = "gene_affected", by.y = "query", all.x = T)

interest <- c("MEN1", "ATRX", "DAXX", "PIK3CA", "PTEN", "TSC2", "VHL", "TP53", "RB1",
"CDKN1B", "APC", "KRAS", "BRAF", "PSIP1", "ARID1A", "SMARCA1", " SMARCA2", "SMARCA4",
"ARID2", "SMARCC2", "SMARCB1", "BCL11A", "SETD1B", "SETDB1", "NSD1", "KDM4A", "PHF8",
"JMJD1", "CBX6", "EZH1", "YY1", "BRWD3", "HDAC4", "STAG1", "NIPBL", "DICER", "ERCC6L",
"KMT2A", "KMT2C", "KMT2D", "ARID1A", "ARID1B", "ARID2", "NOTCH1", "NOTCH2", "NOTCH3",
"NOTCH4", "NF1", "INSR", "TSC2", "RICTOR", "PTPRN", "BEX1", "NOVA1")

mutations <- mutations[mutations$symbol %in% interest,]
mutations$mutated <- rep("mutated",nrow(mutations))
mutations <- mutations[,c("icgc_sample_id", "symbol", "mutated")]
mutations <- mutations[!duplicated(mutations),]


mutations_wider <- mutations %>% pivot_wider(

    names_from = icgc_sample_id,
    values_from = mutated
)

mutations_wider <- t(mutations_wider)
colnames(mutations_wider) <- mutations_wider[1,]
mutations_wider <- mutations_wider[-1,]
mutations_wider[is.na(mutations_wider)] <- "WT"

base <- merge(base,mutations_wider, by.x = "id_sample", 
by.y = "row.names", all.x = TRUE)

base <- base[,-6]

write.table(base,
paste0(outdir, "/base_Lakis.tsv"),
sep = "\t",
row.names = F,
col.names = T)


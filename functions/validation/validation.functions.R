##############################
# Function: normalization
# Input:
#    - matrix: Data.frame with samples as columns and genes as rows
#    - method: Normalization method, either "uq" (upper quartile) or "cpm" (counts per million)
#
# Output:
#    - Normalized matrix with samples as columns and genes as rows
#
# Description:
#    This function normalizes RNA-seq data either using upper quartile normalization or CPM (Counts Per Million) normalization.
##############################


normalization <- function(matrix, method = "uq"){
   # Upper quartile (UQ) normalization
    if (method == "uq") {
        # Apply log2 transformation to the matrix
        matrix <- log2(matrix + 1)
        
        # Calculate the 75th percentile of the expression values for each sample
        data.quantileExpressed <- apply(matrix, 2, function(x) { quantile(x[x > 0], 0.75, na.rm = TRUE) })
        
        # Normalize the matrix using the upper quartile
        matrix <- as.data.frame(t(t(matrix) / data.quantileExpressed)) * 1000
    }
    
    # CPM (Counts Per Million) normalization
    if (method == "cpm") {
        # Create DGEList object for normalization
        dge <- DGEList(matrix)
        
        # Apply TMM normalization (Trimmed Mean of M-values)
        y <- calcNormFactors(dge, method = "TMM")
        
        # Calculate CPM and apply log2 transformation
        cpm <- cpm(y, log = FALSE, normalized.lib.sizes = TRUE)
        matrix <- log2(cpm + 1)
    }
    
    # Return the normalized matrix
    return(matrix)
}


##############################
# Function: generate.sample.sheet
# Input:
#    - path: Directory path containing the IDAT files
#
# Output:
#    - Generates a sample sheet CSV file for methylation analysis
# Description:
#    This function processes the IDAT files in the specified path, organizes them into folders by Sentrix ID, 
#    and generates a sample sheet required for for methylation analysis.
##############################

generate.sample.sheet <- function(path) {

# List all files in the specified directory
files <- list.files(path)

# Filter for files that contain "GSM" in the name
files <- files[grep("GSM",files)]


# Split the filenames into components to extract relevant information 
flist <- strsplit(files, "_")
names <- unlist(lapply(flist,function(x){x[[1]]}))
s_entrix_id <- unlist(lapply(flist,function(x){x[[2]]}))
s_entrix <- unlist(lapply(flist,function(x){x[[3]]}))


# Create directories for each unique Sentrix ID and move files accordingly
s_entrix_id <- unlist(s_entrix_id)
for(i in unique(s_entrix_id)) {

# Create folder if it doesn't already exist
 if (dir.exists(paste0(path,"/",i)) == F) {

    dir.create(paste0(path,"/",i))

}

# Move files to the corresponding folder and rename them
    files.i <- files[grep(i,files)]
    lapply(files.i, function(x){
                
        system(paste0("mv ", path,"/",x, " ",
        path,"/",i),
        )
    # Rename the files by removing the GSM prefix    
        nm <- strsplit(x,"_")
        nm <- nm[[1]]
        nml <- nchar(gsub("GSM","",nm))
        system(paste0("mv ", path,"/",i,"/",x, " ",
        path,"/",i,"/",gsub(paste0("GSM.{",nml,"}_"),"",x)),
        )

    })
}


# Create a sample sheet with required columns for the ChAMP pipeline

sample_sheet <- data.frame("Sample_Name" = names,
"Sample_Well" = rep("",length(names)),
"Sample_Plate" = rep("",length(names)),
"Sample_Group" = rep("",length(names)),
"Pool_ID" = rep("",length(names)),
"Sentrix_ID"   = s_entrix_id,
"Sentrix_Position" = s_entrix
)

# Remove duplicate rows and write the sample sheet to a CSV file
sample_sheet <- sample_sheet[!duplicated(sample_sheet),]

if(file.exists(paste0(gsub("/IDAT","",path),"/Sample_sheet.csv")) == F){

write.table(sample_sheet,
paste0(gsub("/IDAT","",path),"/Sample_sheet.csv"),
sep = ",",
row.names = F,
col.names = T
)
}
}


##############################
# Function: methylation_matrix
# Input:
#    - datadir: Directory containing IDAT files and sample sheet
#    - arraytype: Type of array used (450K or EPIC)
#    - outdir: Output directory for saving processed methylation data
#
# Output:
#    - Writes the normalized methylation data (beta and M values) to CSV files
#
# Description:
#    This function processes methylation data from IDAT files, performs quality control and normalization, 
#    and generates methylation beta and M matrices for further analysis.
# 
# Note: The code in this function is obtained from https://github.com/IARCbioinfo/Methylation_analysis_scripts
##############################

methylation_matrix <- function(datadir, arraytype, outdir){

# Load the appropriate annotation based on the array type
if (arraytype == "450K"){
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
}else{
ann <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

}

# Read the sample sheet and IDAT files
targets <- read.metharray.sheet(datadir)
RGSet <- read.metharray.exp(targets = targets)
targets$ID = paste(targets$Sample_Name)
sampleNames(RGSet) = targets$ID

# Preprocess the raw data and convert to a genomic ratio se
MSet <- preprocessRaw(RGSet)
ratioSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
gset <- mapToGenome(ratioSet, mergeManifest=T)
# Perform quality control and filter samples with low intensity
qc <- getQC(MSet)
# If required: remove any samples with meth or unmeth intensity < 10.5, then remove from all previous files
keep.samples <- apply(as.data.frame(qc),1,mean) > 10.5
RGSet <- RGSet[,keep.samples]
MSet <- MSet[,keep.samples]
ratioSet <- ratioSet[,keep.samples]
gset <- gset[,keep.samples]
targets <- targets[keep.samples,]
# Calculate detection p-values and filter out low-quality samples 
# Here the samples are coloured by sentrix ID to check for poor performing chips
detP <- detectionP(RGSet)
if(sum(colMeans(detP)>0.01)>0 ) print(paste("Warning: sample", names(colMeans(detP))[colMeans(detP)>0.01], "has >0.01 mean detection p-value. Remove sample"))
# If required: remove any samples with detection p value > 0.01, then remove from all previous files
keep.samples <- apply(detP,2,mean) < 0.01 
RGSet <- RGSet[,keep.samples]
MSet <- MSet[,keep.samples]
ratioSet <- ratioSet[,keep.samples]
gset <- gset[,keep.samples]
targets <- targets[keep.samples,]
#Normalization
# includes NOOB background/dye correction
fun <- preprocessFunnorm(RGSet) # If clinical sex is available can use preprocessFunnorm(RGSet, sex=targets$sex), default of sex=NULL uses function getSex to guess the sex
size.fun <- as.data.frame(t(dim(fun)))
# Remove poor performing probes
detP2 = detP[rownames(fun),colnames(fun)]
failed = which(rowSums(detP2 < 0.01) != ncol(fun))
fun1 <- fun[ -failed, ]
print(fun1)
# Check that the correct probes have been removed - so no probe names in failed should be in fun1
failedDF <- as.data.frame(failed)
check <- fun1[featureNames(fun1) %in% rownames(failedDF), ]
if(sum(rownames(check)>0)>0) print(paste("Warning: poor performing probes remain in normalised data"))
size.fun1 <- as.data.frame(t(dim(fun1)))

#Remove XY probes
autosomes = !(rownames(fun1) %in% ann$Name[ann$chr %in% c("chrX","chrY")])
fun2 = fun1[autosomes,]
print(fun2)

size.fun2 <- as.data.frame(t(dim(fun2)))

#Remove SNP-containing probes (mandatory)
print("Remove SNP-associated probes")
fun3 <- dropLociWithSnps(fun2, snps=c("SBE", "CpG"), maf = 0.05) 
size.fun3 <- as.data.frame(t(dim(fun3)))

#Create analysis tables 
betaRaw <- getBeta(gset)
betaNorm <- getBeta(fun3) 
mNorm <- getM(fun3) 

size.betaNorm <- as.data.frame(t(dim(betaNorm)))
size.mNorm <- as.data.frame(t(dim(mNorm)))

# Check for NA and -Inf values
betaRaw.na <- betaRaw[!complete.cases(betaRaw),]
dim(betaRaw.na) 
NAbetas <- rownames(betaRaw.na)
betaNAsNorm <- betaNorm[rownames(betaNorm) %in% NAbetas,] # Should be [1] 0 x ncol(betaNorm)
# Check for infinite values in m table, replace in the no infinite values m table
TestInf <- which(apply(mNorm,1,function(i) sum(is.infinite(i)))>0)
mNoInf <- mNorm
mNoInf[!is.finite(mNoInf)] <- min(mNoInf[is.finite(mNoInf)])
TestInf2 <- which(apply(mNoInf,1,function(i) sum(is.infinite(i)))>0)
TestInf2 # Should be named integer(0)
size.mNoInf <- as.data.frame(t(dim(mNoInf)))
# Write tables
write.table(betaNorm, file=paste0(outdir,
"/Beta_matrix.csv"), sep=",", row.names = T, col.names=NA)
write.table(mNorm, file=paste0(outdir,
"/M_matrix.csv"), sep=",", row.names = T, col.names=NA)
write.table(mNoInf, file=paste0(outdir,
"/M_matrix_noinf.csv"), sep=",", row.names = T, col.names=NA)
}
#Function: qual.met
#
# Input:
#   import: A list generated from the 'champ.import' function, 
#           containing methylation data and detection p-values.
#
# Output:
#   A data frame with sample names as rows and the -log of the mean p-value 
#   for each sample as the only column, named "qual.met". This serves as a 
#   measure of methylation quality.
#
# Description:
#   This function computes the methylation quality for each sample by calculating 
#   the -log of the mean detection p-value. Lower p-values indicate higher quality 
#   methylation data.
qual.met <- function(import) {


    # Extract the detection p-values from the 'import' object
    pval <- import$detP
    # Calculate the -log of the mean p-value for each sample
    qual <- as.data.frame(-log(apply(pval, 2, mean)))
    # Set column name to "qual.met"
    colnames(qual) <- "qual.met"
    
    # Return the quality metrics    
    return(qual)
}



# Function: addCIMPstatus
#
# Input:
#   values: A matrix of methylation beta or M values, where rows represent probes 
#           and columns represent samples.
#
# Output:
#   A data frame with sample names as rows and a "CIMP.status" column, which 
#   classifies samples as "High" or "Low" based on their median methylation 
#   values relative to the global median.
#
# Description:
#   This function assigns a CIMP (CpG Island Methylator Phenotype) status to 
#   each sample by comparing its median methylation value to the overall median 
#   across all samples. Samples with a median greater than or equal to the global 
#   median are classified as "High", while those below are classified as "Low".


addCIMPstatus <- function(values){
    
    # Calculate the global median methylation value across all samples
        median_score <- median(values)
    
    # Classify each sample as "High" or "Low" based on its median value

    f <- function(x){
        m <- median(x,na.rm = T) # Compute the sample's median methylation value
        if (m >= median_score){l <- 'High' }
        if (m <= median_score){l <- 'Low' }
        return(l)
    }
    
    # Apply the classification function across all samples
    CIMP <- apply(values,2,f)

    # Convert to a data frame and set column name to "CIMP.status"
    CIMP <- as.data.frame(CIMP)
    colnames(CIMP) <- "CIMP.status"
    
    # Return the CIMP status data frame    
    return(CIMP)

}


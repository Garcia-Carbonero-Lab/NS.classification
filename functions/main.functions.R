# Function: zscore.rows2
#
# Input:
# - x: A matrix with features as rows and samples as columns, where each 
#       row represents a different feature (e.g., gene expression levels) 
#       and each column corresponds to a sample.
#
# Output:
# A matrix of the same dimensions as the input, where each row has been 
# standardized to have a mean of 0 and a standard deviation of 1. 
# This is achieved by subtracting the mean and dividing by the standard 
# deviation for each feature across samples.
#
# Description:
# This function standardizes each feature (row) in the input matrix by 
# removing the mean and scaling by the standard deviation. The output is 
# useful for subsequent analyses where feature scaling is necessary, such 
# as clustering or principal component analysis (PCA).

zscore.rows2 <- function(x) {
  return(t(apply(x, 1, function(x) (x - mean(na.omit(x))) / sd(na.omit(x)))))
}

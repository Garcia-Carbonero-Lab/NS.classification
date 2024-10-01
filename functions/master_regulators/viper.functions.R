
# Function: RegProcess
# Input:
#      a.file: Character. Path to the ARACNe network file (3-column format: master regulator, target gene, confidence score).
#      exp.mat: Matrix. Expression matrix where genes are in rows and samples are in columns.
#      out.dir: Character. Directory where the output regulon objects (unpruned and pruned) will be saved.
#      out.name: Character. Optional. Prefix for the output filenames (default is ".").
#
# Output:
#      Two RDS files saved in the specified output directory:
#      1. Unpruned regulon object (.rds) - this is the raw network obtained from ARACNe.
#      2. Pruned regulon object (.rds) - network pruned to retain only the strongest interactions (top 50 interactions per master regulator).
#
# Description:
#      This function processes an ARACNe network file to generate a regulatory network (regulon) object using the viper package. 
#      The regulon is then pruned to keep the most significant interactions, and both the unpruned and pruned regulons are saved as RDS files.
#      The pruning process eliminates weak interactions by retaining only the top 50 targets per master regulator.
#
RegProcess <- function(a.file, exp.mat, out.dir, out.name = '.') {
  processed.reg <- aracne2regulon(afile = a.file, eset = exp.mat, format = '3col')
  saveRDS(processed.reg, file = paste(out.dir, out.name, 'unPruned.rds', sep = ''))
  pruned.reg <- pruneRegulon(processed.reg, 50, adaptive = FALSE, eliminate = TRUE)
  saveRDS(pruned.reg, file = paste(out.dir, out.name, 'pruned.rds', sep = ''))
}

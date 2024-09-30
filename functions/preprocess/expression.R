#Function: norm_exp
# 
# Input:
#   targetfile: Path to a tab-separated values (TSV) file containing a column 
#               "FileName" with each row specifying the path to individual .CEL files.
#   inpath:     Path to the directory containing sample annotation files 
#               "sample_anotation.txt" and "Anotacion_Clariom_Human_S.txt".
#   outputpath: Path to the directory where output files will be saved.
#   norm:       Normalization method, one of "rma", "vsn", or "none".
#   flag:       A custom string used as a prefix for the output filenames.
# 
# Output:
#   A normalized expression matrix, with genes as rows and samples as columns, 
#   written to the specified output path.
#
# Description:
#   This function processes Affymetrix .CEL files to generate a normalized 
#   gene expression matrix. Depending on the chosen normalization method ("rma", 
#   "vsn", or none), the data undergoes background correction, normalization, 
#   and summarization. The function also integrates sample annotations, 
#   performs probe-to-gene mapping, and resolves multiple probes mapping to the 
#   same gene by selecting the probe with the highest variance.
#   Finally, the expression matrix is saved in the specified output path.



norm_exp <- function(targetfile, inpath, outputpath, norm = "rma", flag) {
    
    browser()
    # Load target file to obtain .CEL file paths    
    
    Target <- readTargets(targetfile, row.names = "FileName")
    
    # Read in the raw expression data from .CEL files
    data <- ReadAffy(filenames = Target$FileName, cdfname = "clariomshumancdf")


    # Perform normalization based on the selected method    
    if (norm == "vsn") {
        # Apply variance-stabilizing normalization (VSN)
        ndata <- vsnrma(data)
    }

    if (norm == "rma") {

        # Apply Robust Multi-array Average (RMA) normalization
        ndata <- expresso(data,
            bg.correct = TRUE,
            bgcorrect.method = "rma",
            normalize = TRUE,
            normalize.method = "quantiles",
            pmcorrect.method = "pmonly",
            summary.method = "medianpolish",
            verbose = TRUE
        )
    }


    if (norm == "none") {
        # Only background correction, without normalization
        ndata <- expresso(data,
            bg.correct = TRUE,
            bgcorrect.method = "rma",
            normalize = FALSE,
            pmcorrect.method = "pmonly",
            summary.method = "medianpolish",
            verbose = TRUE
        )
    }


    # Harmonize sample IDs by removing unnecessary parts of file names

    data_rma_tab <- as.data.frame(exprs(ndata))
    colnames(data_rma_tab) <- gsub(
        "_\\(Clariom_S_Human\\).CEL", "",
        colnames(data_rma_tab)
    )
     # Clean up NENs sample identifiers
    colnames(data_rma_tab)[grep("NENs", colnames(data_rma_tab))] <- gsub(
        "^[0-9]*\\.", "",
        colnames(data_rma_tab)[grep("NENs", colnames(data_rma_tab))]
    )

    # Load sample annotations
    sample_anot <- read.table(paste0(inpath, "/sample_anotation.tsv"),
        sep = "\t",
        header = T
    )

     # Transpose and merge the expression data with sample annotations
    data_rma_tab <- t(data_rma_tab)



    data_rma_tab <- merge(sample_anot, data_rma_tab,
        by.x = "Samples", by.y = "row.names"
    )

    # Use 'Code' as the row names and clean up the data frame
    rownames(data_rma_tab) <- data_rma_tab$Code
    data_rma_tab <- data_rma_tab[, 3:ncol(data_rma_tab)]
    data_rma_tab <- as.data.frame(t(data_rma_tab))
 
     # Save the raw expression data matrix (after background correction)
    write.table(data_rma_tab,
    paste0(outputpath, "/RMA_probes.txt"),
    sep = "\t",
    row.names = T,
    col.names = NA)
    
    # Load gene annotation file
        anotation <- read.table(paste0(inpath, "/Anotacion_Clariom_Human_S.txt"),
        sep = "\t",
        header = T
    )

    # Merge expression data with gene annotations based on probes
    data_rma_tab <- merge(anotation, data_rma_tab,
        by.x = "Probe",
        by.y = "row.names"
    )

    # Remove probes and select the highest variance probe for each gene
    data_rma_tab <- data_rma_tab[, -1]
    data_rma_tab <- data_rma_tab[complete.cases(data_rma_tab$gene), ]
    num <- apply(data_rma_tab[, 2:ncol(data_rma_tab)], c(1, 2), as.numeric)
    vari <- apply(num, 1, var)
    data_rma_tab <- data_rma_tab[order(vari), ]
    data_rma_tab <- data_rma_tab[!duplicated(data_rma_tab$gene), ]

    # Set gene names as row names and remove 'gene' column
    rownames(data_rma_tab) <- data_rma_tab$gene

    data_rma_tab <- data_rma_tab[, -1]

    # Write the final expression matrix to file
    write.table(data_rma_tab,
        paste0(outputpath, "/", flag, "_Expression_data_genes.tsv"),
        sep = "\t",
        quote = F,
        col.names = NA,
        row.names = T
    )
}
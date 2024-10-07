####################################################
# Function: gsea
#
# Description:
#      This function performs Gene Set Enrichment Analysis (GSEA) using ranked gene 
#      expression data and gene set information in GMT format. It generates a dot plot 
#      for the top enriched pathways with a significance threshold (adjusted p-value < 0.05).
#
# Input:
#      rnk: A named vector of ranked statistics (e.g., log fold-change or t-statistics).
#      gmt: Path to the GMT file containing gene sets (e.g., pathways).
#      title: The title for the output plot file (PDF).
#
# Output:
#      gsea: A data frame with the GSEA results, including pathway enrichment statistics, 
#      normalized enrichment score (NES), and adjusted p-values.
####################################################

gsea <- function(rnk, gmt, title) {
    
    # Load the gene set pathways from the GMT file
    gmt <- gmtPathways(gmt)
    
    # Perform GSEA using the fgsea function
    gsea <- fgsea(
        pathways = gmt,  # Provide the pathways loaded from GMT
        stats = rnk,     # Use the ranked statistics for GSEA
        minSize = 15,    # Minimum gene set size to consider
        maxSize = 500    # Maximum gene set size to consider
    )

    # Concatenate genes in the leading edge into a single string for easier reading
    gsea$leadingEdge <- unlist(lapply(gsea$leadingEdge, function(x) {
        paste(x, collapse = "/")
    }))
    
    # Check if there are significant pathways (adjusted p-value < 0.05)
    if (nrow(gsea[gsea$padj < 0.05, ]) > 0) {
        
        # Get the number of significant pathways (up to 10)
        n <- nrow(gsea[gsea$padj < 0.05, ])
        if (n > 10) {
            n <- 10
        }

        # Filter the results for significant pathways (padj < 0.05)
        gsea <- gsea[which(gsea$padj < 0.05), ]
        
        # Sort the significant pathways by adjusted p-value
        gsea <- gsea[order(abs(gsea$padj)), ]
        
        # Set the pathway names as factors for the plot
        gsea$pathway <- factor(gsea$pathway, levels = gsea$pathway)

        # Create a dot plot for the top significant pathways
        ggplot(data = gsea[1:n, ], aes(x = NES, y = pathway, size = size, color = -1 * log10(padj))) +
            geom_point() +
            scale_size_area(max_size = 10) +   # Scale the dot size by the gene set size
            scale_colour_gradient(low = "green", high = "red") +   # Color by significance
            theme(text = element_text(size = 26))  # Set text size for readability

        # Save the plot to a PDF file
        ggsave(
            filename = title,  # Use the provided title as the filename
            device = "pdf", dpi = 500, width = 22, height = 12
        )

    }

    # Return the GSEA results
    return(gsea)
}

####################################################
# Function: anot.dep
#
# Description:
#      This function annotates a differential methylation (DE) results table with probe 
#      annotations, ensuring that only the most significant unique genes are retained.
#
# Input:
#      dep: Data frame of differential methylation results with gene identifiers and adjusted p-values.
#      anot: Data frame containing probe annotations with probes IDs and gene names.
#
# Output:
#      dep: Annotated and filtered differentially methylated results with genes as row names.
####################################################

anot.dep <- function(dep, anot) {
    
    # Merge DE results with gene annotations using gene IDs
    dep <- merge(dep, anot, by.x = "row.names", by.y = "ID")
    
    # Sort the merged table by adjusted p-value (adj.P.Val)
    dep <- dep[order(dep$adj.P.Val),]
    
    # Remove duplicate genes, keeping only the first occurrence
    dep <- dep[!duplicated(dep$genes),]
    
    # Set the gene names as row names for easier reference
    rownames(dep) <- dep$genes
    
    # Return the annotated and filtered DE results
    return(dep)
}
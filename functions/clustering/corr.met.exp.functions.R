# Function: corr.function
# Input:
#    - exp: Gene expression matrix with genes as rows and samples as columns
#    - met: Methylation matrix with probes as rows and samples as columns
#    - id: ID of the methylation probe
#    - gene: Gene symbol
#
# Output:
#    - A list containing correlation coefficient (estimate), r-squared, and p-value
#      for the correlation between the gene expression and methylation probe.
#
# Description:
#    This function calculates the correlation between a gene's expression and 
#    its corresponding methylation probe using Pearson correlation and a linear model.
##############################
corr.function <- function(exp, met, id, gene){

    # Combine the expression and methylation data for the specified gene and probe
    df <- cbind(
        as.numeric(exp[gene,]),  # Expression data for the gene
        as.numeric(met[id,])     # Methylation data for the probe
    )
    
    # Convert to a dataframe
    df <- as.data.frame(df)
    
    # Assign column names
    colnames(df) <- c("exp", "met")
    
    # Perform correlation test (Pearson)
    cor <- cor.test(as.numeric(df[,"exp"]), as.numeric(df[,"met"]))
    
    # Perform linear regression and extract the adjusted R-squared value
    lmx <- summary(lm(exp ~ met, df))
    
    # Get adjusted r-squared, p-value, and correlation estimate
    r2 <- lmx$adj.r.squared[[1]]
    p.val <- cor$p.value[[1]]
    est <- cor$estimate[[1]]
    
    # Return the result as a list
    return(list("gene" = gene,
                "CpG" = id,
                "corr" = est,
                "r2" = r2,
                "p.val" = p.val))
}

##############################
# Function: corr.met.exp.unique
# Input:
#    - exp: Gene expression matrix with genes as rows and samples as columns
#    - met: Methylation matrix with probes as rows and samples as columns
#    - met_anot: Annotation matrix of methylation probes
#
# Output:
#    - A dataframe with the correlation coefficient, r-squared, and p-value for each 
#      probe-gene pair from the correlation analysis.
#
# Description:
#    This function correlates gene expression with methylation probes, filtering the data
#    by overlapping genes, removing duplicates, and then calculating correlations for unique gene-probe pairs.
##############################
corr.met.exp.unique <- function(exp, met, met_anot){

    # Separate genes associated with probes if multiple genes are listed in the same row
    met_anot <- met_anot %>%
        separate_rows(genes, sep = "_")
    
    # Convert gene IDs to gene symbols
    symbol_genes <- alias2SymbolTable(rownames(exp))
    
    # Create a dataframe with gene symbols, aliases, and variances
    anot <- data.frame(
        "symbol" = symbol_genes,
        "alias" = rownames(exp),
        "variance" = apply(exp, 1, var)
    )
    
    # Sort data based on variance for filtering duplicates
    exp <- exp[order(anot$variance),]
    anot <- anot[order(anot$variance),]
    
    # Remove duplicate gene symbols
    exp <- exp[!duplicated(anot$symbol),]
    anot <- anot[!duplicated(anot$symbol),]
    
    # Replace rownames with gene symbols where available
    rownames(exp)[!is.na(anot$symbol)] <- anot$symbol[!is.na(anot$symbol)]
    
    # Convert probe annotation to gene symbols
    symbol_met <- alias2SymbolTable(met_anot$genes)
    met_anot$genes[!is.na(symbol_met)] <- symbol_met[!is.na(symbol_met)]
    
    # Filter annotation data to keep only overlapping genes and probes
    met_anot <- met_anot[met_anot$genes %in% rownames(exp),]
    met_anot <- met_anot[met_anot$ID %in% rownames(met),]
    
    # Filter methylation matrix for probes in the annotation
    met <- met[met_anot$ID,]
    
    # Ensure expression and methylation matrices have matching samples
    exp <- exp[, colnames(met)]
    
    # Perform correlation analysis for each probe
    corr.df <- lapply(met_anot$ID, function(x){
        gene <- met_anot$genes[met_anot$ID == x]
        id <- x
        
        # Calculate correlation for each gene-probe pair
        result <- lapply(gene, function(y){
            return(corr.function(exp, met, id, y))
        })
        
        # Convert list of results to dataframe
        result <- as.data.frame(do.call(rbind, result))
        return(result)
    })
    
    # Combine all correlation results into a single dataframe
    corr.df <- do.call(rbind, corr.df)
    
    # Ensure all columns are properly formatted
    for (col in 1:ncol(corr.df)){
        corr.df[, col] <- unlist(corr.df[, col])
    }
    
    return(corr.df)
}

##############################
# Function: correlation.plot
# Input:
#    - exp: Gene expression matrix with genes as rows and samples as columns
#    - met: Methylation matrix with probes as rows and samples as columns
#    - anot: Annotation matrix of methylation probes with associated genes
#    - base: Clinical data or sample metadata with samples as rows
#    - genes: List of gene symbols to plot
#    - outdir: Directory path for saving the output plots
#    - group: Variable (column in `base`) used for grouping the samples in the plots
#    - color.group: Colors for each group
#    - flag: Additional identifier for output filenames
#
# Output:
#    - Correlation plots saved as PDF files.
#
# Description:
#    This function generates scatter plots with linear regression lines to visualize the correlation 
#    between gene expression and methylation probe data, faceted by the group variable.
##############################
correlation.plot <- function(exp, met, anot, base, genes, outdir, group, color.group, flag) {
    
    # Replace "-" with "_" in gene names for consistency
    genes <- gsub("-", "_", genes)
    rownames(exp) <- gsub("-", "_", rownames(exp))
    anot$genes <- gsub("-", "_", anot$genes)
    
    # Filter expression and methylation data to include only samples present in `base`
    exp <- exp[, rownames(base)]
    met <- met[, rownames(base)]
    
    # Create directory for saving correlation plots if it doesn't exist
    if (!dir.exists(paste0(outdir, "/corr.plots"))) {
        dir.create(paste0(outdir, "/corr.plots"))
    }
    
    # Iterate over each gene and generate plots for associated probes
    for (gene in genes){
        print(gene)
        
        # Get associated probes for the current gene
        probes <- anot$ID[anot$genes == gene]
        probes <- probes[probes %in% rownames(met)]
        
        for (probe in probes){
            
            # Create dataframe combining expression, methylation, and group data
            df <- data.frame(
                group = base[, group],
                "exp" = as.numeric(exp[gene,]),
                "met" = as.numeric(met[probe,]),
                row.names = rownames(base)
            )
            
            gene_met <- paste(gene, probe, sep = "_")
            colnames(df) <- c(group, gene, gene_met)
            df[, group] <- as.factor(df[, group])
            
            # Generate scatter plot with linear regression line and correlation coefficient
            plt1 <- ggplot(df, aes_string(x = gene, y = gene_met, colour = group, group = 1)) +
                theme_classic() +
                geom_point(size = 3) + 
                scale_color_manual(values = color.group) +
                stat_cor(size = 10) +  # Add correlation coefficient
                stat_smooth(method = "lm", color = "black") +  # Linear regression line
                theme(
                    text = element_text(size = 26),
                    legend.key.size = unit(2, "cm"),
                    legend.text = element_text(size = 26)
                )
            
            # Save the plot as a PDF
            ggsave(
                plot = plt1, 
                filename = paste0(outdir, "/corr.plots/", gene, "_", probe, "_", flag, ".pdf"),
                device = "pdf", 
                dpi = 500, 
                width = 17, 
                height = 12
            )
            
            # Generate faceted plot with linear regression and correlation by group
            plt2 <- ggplot(df, aes_string(x = gene, y = gene_met)) +
                theme_classic() + 
                stat_cor(size = 10) +  # Add correlation coefficient
                stat_smooth(method = "lm") +  # Linear regression line
                geom_point(size = 3) +
                facet_wrap(as.formula(paste0("~ ", group)), drop = TRUE, scales = "free_x", nrow = 1) +
                theme(
                    text = element_text(size = 26),
                    legend.key.size = unit(2, "cm"),
                    legend.text = element_text(size = 26)
                )
            
            # Save the faceted plot as a PDF
            ggsave(
                plot = plt2, 
                filename = paste0(outdir, "/corr.plots/", gene, "_", probe, "_", flag, ".wrapped.pdf"),
                device = "pdf", 
                dpi = 500, 
                width = 25, 
                height = 12
            )
        }
    }
}
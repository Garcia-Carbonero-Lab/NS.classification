# Function: theme_general
# Input:
#      plot: ggplot object
# Output:
#      ggplot object with modified theme
# Description: Applies general theme settings to a ggplot object. 
#              Specifically, it modifies the text size and the legend size for clarity.


theme_general <- function(plot) {
    themed <- plot +
        theme(
            text = element_text(size = 15), # Set global text size to 15
            legend.key.height = unit(10, "mm") # Set the height of the legend key to 10 mm
        )
        
}


# Function: pca_corr
# Input:
#      data: Matrix with features as rows and samples as columns
#      base: Data frame with clinical information where samples are rows and features are columns
#      columns: Vector with the names of continuous columns in clinical data to correlate
#      outpath: Output directory path to save results
#      scale: Logical. If TRUE, data is scaled by removing median and dividing by standard deviation
#      flag: Identifier string added to the result files
#      anotate1: Coordinates (numeric vector) to place the text annotation for PCA1 correlation
#      anotate2: Coordinates (numeric vector) to place the text annotation for PCA2 correlation
# Output:
#      A matrix with features as rows and the first 10 principal components as columns.
#      It includes correlation coefficients and p-values calculated using Pearson correlation.
# Description: Performs a PCA on the input data, correlates  continuous variables with PCA axes,
#              and creates a scatterplot of PC1 vs PC2. Correlation values and p-values for continous variables
#              are shown in the plots.


pca_corr <- function(data, base, columns, outpath, scale = T, flag,
                     anotate1 = c(80, 150),
                     anotate2 = c(80, 140)) {
    if (scale == T) {
        data <- zscore.rows2(data)
    }

    # Perform PCA
    pcal <- prcomp(t(data), scale. = F)
    pca <- as.data.frame(pcal$x[, 1:10])
    prop_varianza <- pcal$sdev^2 / sum(pcal$sdev^2)

    # Prepare clinical data for correlation
    base_num <- base[, columns,drop = F]
    base_num <- merge(base_num, pca, by = "row.names")

    corrl <- list() # Store correlation results
    est <- list() # Store correlation estimates


    
    # Correlate each clinical variable with each principal component
    for (c in columns) {
        for (p in colnames(pca)) {

            # We calculate correlation between PCAs and variables
            bs <- base_num[complete.cases(base_num[c]), ] # Remove missing data
            corr <- cor.test(bs[, p], bs[, c]) # Perform Pearson correlation
            corrl[[c]][[p]] <- paste0(
                "corr=", round(corr$estimate[[1]], 2),
                "/pval=", round(corr$p.value, 6)
            )
            est[[c]][[p]] <- round(corr$estimate, 2) # Save the correlation estimate
        }
        # Create scatterplot of PC1 vs PC2, coloring by clinical variable value
        plt <- ggplot(aes_string("PC1", "PC2", color = c), data = bs) +
            geom_point() +
            theme_classic() +
            xlab(paste("PCA1", round(prop_varianza[[1]], 2))) +
            ylab(paste("PCA2", round(prop_varianza[[2]], 2))) +
            scale_color_gradient2(
                midpoint = median(base_num[, c], na.rm = T),
                low = "blue", mid = "white", high = "red",
                name = c
            ) +
            annotate("text",
                x = anotate1[1], y = anotate1[2], size = 4,
                label = paste0(
                    "correlation PC1 = ",
                    est[[c]][["PC1"]]
                )
            ) +

            annotate("text",
                x = anotate2[1], y = anotate2[2], size = 4,
                label = paste0(
                    "correlation PC2 = ",
                    est[[c]][["PC2"]]
                )
            ) +

    theme(text=element_text(size=20),
          legend.key.height = unit(1,'cm'),
          aspect.ratio=1/1
          )

    
    # Save plot as a PDF file

        ggsave(plot = plt, filename = paste0(
            outpath, "/", c,
            "_", flag, "_PCA.pdf"
        ), height = 8, width = 8, device = "pdf"
           ,dpi =  300)
    }


    df_corr <- do.call(rbind, corrl)

    # Save correlation results to a text file
    write.table(df_corr,
    paste0(outpath,"/Correlation_",flag,".txt"),
    sep = "/t",
    row.names = T,
    col.names = NA)
    return(plt)
}


# Function: generate_pca
#
# Input:
# - data: A matrix with features as rows and samples as columns. Each 
#         row represents a different feature (e.g., gene expression levels), 
#         and each column corresponds to a sample.
#
# Output:
# A data frame containing the principal component scores for the first 
# principal component (PC1) for each sample. The output will have samples 
# as rows and a single column corresponding to PC1.
#
# Description:
# This function performs Principal Component Analysis (PCA) on the input 
# data matrix. It first transposes the matrix to ensure that PCA is applied 
# across samples. The data is scaled before PCA is computed to standardize 
# the variance across features. The output consists of the first principal 
# component scores for each sample

generate_pca <- function(data) {
    pca_data <- prcomp(t(data), scale. = TRUE)
    pca_data <- as.data.frame(t(pca_data$x))
    pca_data <- t(pca_data)
    pca_data <- pca_data[, 1, drop = FALSE]
    return(pca_data)
}



# Function: filter_tissue
#
# Input:
# - data: A matrix with features as rows and samples as columns, representing the 
#         expression levels of various genes across different samples.
# - base: A data frame containing clinical information with samples as rows and 
#         clinical features (e.g., biopsy location, biopsy site) as columns.
# - tissues: A character vector specifying the tissues from which features should 
#            be selected (e.g., c("Lung", "Colorectal")).
# - n: An integer indicating the number of features to select based on their 
#      Mean Absolute Deviation (MAD). Default is 5000.
# - flag: A string to be included in the output filename for identification 
#         purposes.
#
# Output: 
# A matrix with features as rows and selected features based on MAD. The output will 
# also include a file containing the MAD values for each tissue, with a global 
# median across all tissues.
#
# Description: 
# This function filters features based on their MAD values computed from the specified 
# tissues in primary tumor samples. It writes the MAD results to a text file and 
# returns the subset of features with the highest MAD values for further analysis.

filter_tissue <- function(data, base, tissues, n = 5000, flag) {
    
    # Apply selection for each primary tumor separately
    mad_tissue <- lapply(tissues, function(x) {
    # Identify patients corresponding to the specified tissue and primary tumor sites
    patients <- rownames(base)[base$BIOPSY_LOCATION == x & base$BIOPSY_SITES == "PRIMARY"]
        
        # Filter data for identified patients
        df <- data[, colnames(data) %in% patients]
        
        # Calculate Mean Absolute Deviation (MAD) for each feature across samples
        v <- apply(df, 1, mad)
        
        # Order features by MAD in descending order
        df <- df[order(v, decreasing = T), ]


        # Create a data frame of ordered features
            
        res <- data.frame(rownames(df), 1:nrow(df))
        colnames(res) <- c("genes", x)
        res <- res[order(res$genes), ]
        return(res)
    })

    # Combine MAD results from all tissues
    mad_tissue <- bind_cols(mad_tissue)
    rownames(mad_tissue) <- mad_tissue$genes...1
    mad_tissue <- mad_tissue[, grep("genes", invert = T, colnames(mad_tissue))]
    mad_tissue$GLOBAL <- apply(mad_tissue, 1, median)
    
    
    # Write MAD results to a text file
    write.table(mad_tissue, paste0(flag, "_MAD_TISSUE.txt"),
        sep = "\t", row.names = T, col.names = NA
    )


    # Select top n genes based on global median MAD
    genes <- rownames(mad_tissue[order(mad_tissue$GLOBAL), ])[1:n]
    data <- data[genes, ]
    return(data)
}
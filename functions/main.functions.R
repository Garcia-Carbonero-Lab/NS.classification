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

# function: create.met.anot
# input: 
#   - epic: EPIC arrays annotation (data.frame)
#
# output: 
#   - data.frame with probe names as row.names and columns with information 
#     about chromosome, start and end position, Island, fgroup, DMR genes, 
#     and UCSC reference group.
#
# description: 
#   Takes the EPIC annotation file and transforms it into a simplified 
#   annotation matrix that includes details about the probes, genes, 
#   genomic coordinates, CpG island information, regulatory feature groups, 
#   DMR (differentially methylated region) genes, and UCSC reference groups.
#   The result is a data.frame with one row per probe and the relevant metadata.

create.met.anot <- function(epic) {

    # Extract gene names and UCSC reference groups from the annotation, splitting by ';'
    genes <- strsplit(as.character(epic$UCSC_RefGene_Name), ';')
    genes[lengths(genes) == 0] <- NA  # Set empty lists to NA

    group <- strsplit(as.character(epic$UCSC_RefGene_Group), ';')
    group[lengths(group) == 0] <- NA  # Set empty lists to NA

    # Create the main data.frame with relevant columns from the annotation
    met_anot <- data.frame(
        'ID' = rep(epic$IlmnID, sapply(genes, length)),  # Probe IDs
        "genes" = unlist(genes),  # Unlist gene names
        "CpG_chrm" = rep(epic$CHR_hg38, sapply(genes, length)),  # Chromosome
        "CpG_beg" = rep(epic$Start_hg38, sapply(genes, length)),  # Start position
        "CpG_end" = rep(epic$End_hg38, sapply(genes, length)),  # End position
        "UCSC_RefGene_Group" = unlist(group),  # UCSC reference gene group
        "ISLAND" = rep(epic$Relation_to_UCSC_CpG_Island, sapply(genes, length)),  # CpG Island relation
        "FGROUP" = rep(epic$Regulatory_Feature_Group, sapply(genes, length)),  # Regulatory feature group
        "DMR" = rep(epic$DMR, sapply(genes, length))  # DMR (differentially methylated region) information
    )

    # Remove duplicate rows where both probe ID and gene are identical
    met_anot <- met_anot[!duplicated(paste(met_anot$ID, met_anot$genes)),]

    # Summarize the data by probe ID, combining genes and UCSC reference groups into a single entry
    met_anot <- met_anot %>% 
        group_by(ID, CpG_chrm, CpG_beg, CpG_end, ISLAND, FGROUP, DMR) %>%
        dplyr::summarise(
            genes = toString(genes),  # Concatenate gene names
            UCSC_RefGene_Group = toString(UCSC_RefGene_Group)  # Concatenate UCSC gene groups
        )

    # Replace commas with underscores in the gene names to ensure unique formatting
    met_anot$genes <- gsub(', ', '_', met_anot$genes)

    # Set the row names of the data frame as the probe IDs
    rownames(met_anot) <- met_anot$ID

    # Return the final annotation data.frame
    return(met_anot)
}

# Function to apply a general theme to ggplot2 plots

# function: theme_general
# input:
#   - plot: a ggplot2 object (the plot to which the theme will be applied)
#
# output:
#   - themed: a ggplot2 object with the applied theme settings
#
# description:
#   This function customizes the appearance of a given ggplot2 plot by applying 
#   a general theme. It modifies text size, adjusts the angle and alignment of 
#   x-axis text, sets the height of legend keys, and defines the aspect ratio of 
#   the plot to be square (1:1). The function returns the modified plot object 
#   for further use or display.

theme_general <- function(plot) {
    # Apply custom theme settings to the ggplot2 plot
    themed <- plot +
        theme(text = element_text(size = 40),  # Set base text size for the plot
              axis.text.x = element_text(angle = 0, vjust = 1, hjust = 0.5),  # Customize x-axis text orientation and positioning
              legend.key.height = unit(2, 'cm'),  # Set height of legend keys
              aspect.ratio = 1/1  # Set aspect ratio of the plot to be square
              )
    return(themed)  # Return the themed plot object
}

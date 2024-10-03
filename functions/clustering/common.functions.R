
####################################################
# Function: heatmap
####################################################
# Input:
#      data: matrix with features as rows and samples as columns
#      base: data frame with samples as rows and features as columns (used for annotations)
#      outdir: path to save the heatmap PDF output
#      col: list indicating colors of annotations
#      size: vector indicating width and height of the heatmap plot
#      anot_param_legend: list indicating parameters of the annotation legend
#      flag: string to append to the output file name
#      max: upper limit for color scaling; if NULL, it selects the maximum in the matrix
#      min: lower limit for color scaling; if NULL, it selects the minimum in the matrix
#      names.size: size of row names in the heatmap
#      column: column in 'base' to split samples for heatmap annotation
#      anot_h: height for annotation elements
#      annotation_name_gp: size of annotation names in the plot
#      show_row_names: boolean indicating whether to show row names (TRUE or FALSE)
#      color.max: color for the maximum values in the heatmap
#      color.min: color for the minimum values in the heatmap
#      z.score: boolean indicating whether to apply z-score normalization on rows
#      row.split: factor to split rows for grouping (optional)
# Output:
#      - A heatmap plot saved as a PDF
# Description:
#      This function generates a heatmap with annotated sample groups and customizable 
#      color scales. It supports z-score normalization, automatic min/max scaling, and 
#      annotation for specific groupings. The heatmap is saved as a PDF file.
####################################################

heatmap <- function(data,
base,
outdir, col, size,
anot_param_legend,
flag,
max = NULL,
min = NULL,
names.size = 3,
column,
anot_h,
annotation_name_gp = 60,
show_row_names = T,
color.max = "red",
color.min = "blue",
z.score = T,
row.split = NULL) {

# Ensure data is ordered to match rows in base
data <- data[, rownames(base)]

# If z-score normalization is requested, apply it to the rows
if(z.score == T){
data <- zscore.rows2(data)
}

# Prepare the annotation data using the provided base and color list
anot <- base[, names(col)]

# Create the column annotation object for the heatmap
ha_column <- HeatmapAnnotation(
        df = anot,
        show_legend = TRUE,
        col = col,
        annotation_name_gp = gpar(fontsize = annotation_name_gp),  # Size of annotation names
        border = unit(40, "cm"),  # Border size for annotation
        annotation_height = anot_h,  # Annotation height
        height = unit(10, "cm"),  # General height for the annotation
        gp = gpar(fontsize = 80),  # Font size for annotation
        simple_anno_size_adjust = TRUE,  # Adjust annotation size dynamically
        annotation_legend_param = anot_param_legend  # Legend parameters
    )

# Set max and min values for color scaling if they are not provided
    if (is.null(max)) {
        max <- max(data)
    }
    if (is.null(min)) {
        min <- min(data)
    }

    # Define the color gradient for the heatmap using the provided max and min colors
    col <- colorRamp2(
        c(max, 0, min),
        c(color.max, "white", color.min)
    )


    # Generate the heatmap without row splitting
    if (is.null(row.split)) {
        hm1 <- Heatmap(as.matrix(data),
                       top_annotation = ha_column,
                       show_row_names = show_row_names,  # Control visibility of row names
                       show_column_names = FALSE,
                       col = col,  # Apply color scaling
                       gap = unit(5, "mm"),  # Set gap between heatmap sections
                       cluster_rows = TRUE,  # Enable row clustering
                       show_row_dend = FALSE,  # Hide row dendrogram
                       show_column_dend = FALSE,  # Hide column dendrogram
                       column_dend_side = "top",  # Position column dendrogram at the top
                       column_dend_height = unit(6, "cm"),  # Set height of column dendrogram
                       cluster_columns = TRUE,  # Enable column clustering
                       clustering_distance_rows = "pearson",  # Distance metric for rows
                       clustering_distance_columns = "pearson",  # Distance metric for columns
                       heatmap_legend_param = list(
                           title = "Z score",  # Legend title
                           title_position = "topcenter",  # Title position in legend
                           legend_direction = "horizontal",  # Legend layout direction
                           title_gp = gpar(fontsize = 60),  # Title font size
                           labels_gp = gpar(fontsize = 60,  # Label font settings
                                            family = "Times", face = "bold"),
                           legend_width = unit(30, "cm"),  # Width of the legend
                           grid_height = unit(2.5, "cm")  # Height of legend grid
                       ),
                       row_names_side = "left",  # Position row names to the left
                       row_names_gp = gpar(fontsize = 60),  # Row names font size
                       column_split = anot[, column],  # Split columns by the selected annotation
                       row_names_max_width = unit(names.size, 'cm'),  # Control max row names width
                       cluster_column_slices = FALSE  # Disable clustering of column slices
        )
    } else {
        # Generate the heatmap with row splitting
        hm1 <- Heatmap(as.matrix(data),
                       top_annotation = ha_column,
                       show_row_names = show_row_names,
                       show_column_names = FALSE,
                       col = col,
                       gap = unit(5, "mm"),
                       cluster_rows = TRUE,
                       show_row_dend = FALSE,
                       show_column_dend = FALSE,
                       column_dend_side = "top",
                       column_dend_height = unit(6, "cm"),
                       cluster_columns = TRUE,
                       clustering_distance_rows = "pearson",
                       clustering_distance_columns = "pearson",
                       heatmap_legend_param = list(
                           title = "Z score",
                           title_position = "topcenter",
                           legend_direction = "horizontal",
                           title_gp = gpar(fontsize = 60),
                           labels_gp = gpar(fontsize = 60,
                                            family = "Times", face = "bold"),
                           legend_width = unit(30, "cm"),
                           grid_height = unit(2.5, "cm")
                       ),
                       row_names_side = "left",
                       row_names_gp = gpar(fontsize = 60),
                       column_split = anot[, column],
                       row_split = row.split,  # Specify row splitting
                       row_names_max_width = unit(names.size, 'cm'),
                       cluster_column_slices = FALSE
        )
    }

    # Combine the heatmap into a list for plotting
    hlist <- hm1

    # Draw the heatmap with customized legend positions and padding
    p <- draw(hlist, 
              heatmap_legend_side = "bottom",  # Position heatmap legend at the bottom
              annotation_legend_side = "top",  # Position annotation legend at the top
              padding = unit(c(20, 240, 20, 20), "mm")  # Set padding for the plot
    )

    # Save the heatmap as a PDF file
    pdf(paste0(outdir, "/Heatmap_anotated_", flag, ".pdf"),
        width = size[["width"]], size[["height"]])
    print(p)
    dev.off()
}



####################################################
# Function: wilcox_plots
####################################################
# Inputs:
#      data: Matrix with features as rows and samples as columns
#      base: Data frame with samples as rows and features as columns (used for group information)
#      group.x: Name of the column in 'base' to use for the x-axis in the boxplots
#      group.fill: Name of the column in 'base' to use for coloring the boxplots (grouping factor)
#      color.fill: Vector of colors to use for each level in 'group.fill'
#      outdir: Path to save the output plots and results
#      flag: String to append to file names (used to distinguish different plots)
# Outputs:
#      - Boxplots for each feature, showing comparisons between groups
#      - Wilcoxon rank-sum test results saved as a TSV file
# Description:
#      This function generates boxplots to compare groups for each feature in the data
#      using the Wilcoxon rank-sum test. Adjusts p-values using FDR and annotates the 
#      plots with significance levels. Saves the results and plots as PDFs and a TSV.
####################################################

wilcox_plots <- function(data,
base,
group.x,
group.fill,
color.fill,
outdir,
flag) {
 # Ensure that the rows of the data matrix correspond to the samples in 'base'
data <- data[rownames(base),,drop = F]

# Combine the data matrix with the group information from 'base'
df<- cbind(data,base[,c(group.x,group.fill)])

# Convert group columns to factors
df[,group.x] <- as.factor(as.character(df[,group.x]))
df[,group.fill] <- as.factor(as.character(df[,group.fill]))

# List to store the results from Wilcoxon tests
stats <- list()

# Loop over each feature in the data (columns)
for (row in colnames(data)) {

# Perform the Wilcoxon test based on whether group.x and group.fill are different
if (group.x != group.fill){
stat.test <- df %>%
  group_by(.dots=group.x) %>%
  wilcox_test(as.formula(paste0(row, "~", group.fill))) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj")

   median <- df[,c(row,group.x, group.fill)] %>% 
    group_by_at(vars(one_of(c(group.x,group.fill)))) %>%
    get_summary_stats(type = "median")
}else{
stat.test <- df %>%
  wilcox_test(as.formula(paste0(row, "~", group.fill))) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj")

  median <- df[,c(row, group.fill)] %>% 
    group_by_at(vars(group.fill)) %>%
    get_summary_stats(type = "median") 

}

 # Prepare median results for both groups
median_l <- lapply(1:nrow(stat.test), function(x){
dfmedian <- data.frame(median_group_1 = 
median$median[median[,group.fill] == stat.test[[x,"group1"]]],
median_group_2 = 
median$median[median[,group.fill] == stat.test[[x,"group2"]]]
)
})

median <- do.call(rbind,median_l)
# Store test statistics (Wilcoxon results + median values)
stats[[row]] <- as.data.frame(cbind(stat.test,median))
# Add p-value positions for the plot
stat.test <- stat.test %>%
  add_xy_position(x = group.x, dodge = 0.8)

# Generate a boxplot for the feature
bxp <- ggboxplot(
  df, x = group.x, y = row,
  fill = group.fill,
  width = 0.6 ,
  palette = color.fill
  ) +


  stat_pvalue_manual(
  stat.test,  label = "p.adj.signif", tip.length = 0, hide.ns = T
) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +

theme(
        text = element_text(size = 26),
        legend.key.size = unit(2, "cm"),
        legend.text = element_text(size = 26)
        ) +

        labs(x = "", y = row)

    bxp$layers[[2]]$aes_params$size <- 1
    bxp$layers[[2]]$aes_params$label.size <- 20
# Save the boxplot to a PDF file
ggsave(
        plot = bxp, filename = paste0(
            outdir, "/Boxplot.primary.subtype_", row, ".pdf"
        ),
        device = "pdf",
        dpi = 300, width = 16,
        height = 12
    )


}

stats <- do.call(rbind,stats)

# Combine all statistics into a single data frame and save as a TSV
write.table(stats,
paste0(outdir, "/Wilcox_", flag, ".txt"),
sep = "\t",
row.names = F,
col.names = T)

}

####################################################
# Function: heatmap.molecular
####################################################

# Input:
#      exp: matrix with genes as rows and samples as columns (gene expression data)
#      met: matrix with methylation probes as rows and samples as columns (methylation data)
#      base: data frame with samples as rows and annotations as columns
#      outdir: path to save the heatmap PDF output
#      col: list indicating colors for annotations
#      size: vector with width (first term) and height (second term) for the heatmap plot
#      anot_param_legend: list indicating parameters for annotation legend
#      flag: string to append to the output file name
#      max.exp: max value for color scaling in expression data; NULL auto-selects the maximum
#      min.exp: min value for color scaling in expression data; NULL auto-selects the minimum
#      max.met: max value for color scaling in methylation data; NULL auto-selects the maximum
#      min.met: min value for color scaling in methylation data; NULL auto-selects the minimum
#      names.size: size of row names in the heatmap
#      column: column in `base` to split samples for heatmap annotation
#      anot_h: height for annotation elements
#      annotation_name_gp: size of annotation names
#      show_row_names: boolean, TRUE to show row names, FALSE to hide
#      scale_met: boolean, TRUE to apply z-score normalization on methylation data
#      scale_exp: boolean, TRUE to apply z-score normalization on expression data
#      height_anotation: height for annotation block
#      cluster: choose either "exp" or "met" to select which dataset to use for clustering
#      order: order of the dendrogram for sample grouping
# Output:
#      - A heatmap plot saved as a PDF file with both expression and methylation data
# Purpose:
#   This function generates a heatmap that visualizes gene expression and 
#   methylation data, incorporating sample annotations. It applies z-score 
#   normalization to the data if specified, clusters samples based on either 
#   expression or methylation data, and outputs the results as a PDF file.
#

####################################################

heatmap.molecular <- function(exp,
met,
base,
outdir, col, size,
anot_param_legend,
flag,
max.exp = NULL,
min.exp = NULL,
max.met = NULL,
min.met = NULL,
names.size = 3,
column,
anot_h,
annotation_name_gp = 60,
show_row_names = T,
scale_met = T,
scale_exp = T,
height_anotation = 25,
cluster = "exp",
order) {


# Ensure expression data matches the sample order in 'base'
exp <- exp[, rownames(base)]

# Apply z-score normalization to expression data if requested
base[,column] <- as.factor(as.character(base[,column]))


if (scale_exp == T){
exp <- zscore.rows2(exp)
}


# Ensure methylation data matches the sample order in 'base'
met <- met[, rownames(base)]
# Apply z-score normalization to methylation data if requested
if (scale_met == T){
met <- zscore.rows2(met)
}

# Cluster samples based on either expression or methylation data
if (cluster == "exp"){
cluster_cols <- cluster_within_group(exp,base[,column])

}

if (cluster == "met"){
cluster_cols <- cluster_within_group(met,base[,column])
}

# Reorder base data and expression/methylation matrices based on dendrogram order
base <- base[labels(cluster_cols),]
base <- base[c(grep(order[[1]],base[,column]),
grep(order[[2]],base[,column]),
grep(order[[3]],base[,column])),]
exp <- exp[, rownames(base)]
met <- met[, rownames(base)]

# Create annotation data for the heatmap
anot <- base[, names(col)]

# Define the heatmap annotation block for samples
ha_column <- HeatmapAnnotation(
            df = anot,
            show_legend = T,
            col = col,
            annotation_name_gp = gpar(fontsize = annotation_name_gp),
            border = unit(40, "cm"),
            annotation_height = anot_h,
            height = unit(height_anotation, "cm"),
            gp = gpar(fontsize = 80),
            simple_anno_size_adjust = TRUE,
            annotation_legend_param = anot_param_legend,
            simple_anno_size = unit(4, "cm"))

# Set maximum and minimum values for color scaling in expression and methylation data
if (is.null(max.exp)){

    max.exp <- max(exp)
}
if (is.null(max.met)){

    max.met <- max(met)
}

if (is.null(min.exp)){

    min.exp <- min(exp)
}

if (is.null(min.met)){

    min.met <- min(met)
}

# Define color gradients for expression and methylation data
col.exp <- colorRamp2(
                c(max.exp, 0, min.exp),
                c("red", "white", "blue")
)



col.met <- colorRamp2(
                c(max.met, 0, min.met),
                c("green", "white", "blue")
)

# Create the heatmap for gene expression data
hm1 <- Heatmap(as.matrix(exp),
            top_annotation = ha_column,
            show_row_names = show_row_names,
            show_column_names = T,
            col = col.exp,
            gap = unit(5, "mm"),
            cluster_rows = T,
            show_row_dend = F,
            show_column_dend = F,
            cluster_columns = F,
            clustering_distance_rows = "pearson",
            clustering_distance_columns = "pearson",
            heatmap_legend_param = list(
                title = "expression Z score",
                title_position = "topcenter",
                legend_direction = "horizontal",
                title_gp = gpar(fontsize = 60),
                labels_gp = gpar(fontsize = 60,
                family = "Times", face = "bold"),
                legend_width = unit(30, "cm"),
                grid_height = unit(2.5, "cm")
            ),
            row_names_side = "left",
            row_names_gp = gpar(fontsize = 100),
            row_names_max_width = unit(names.size,'cm'),
            cluster_column_slices = FALSE,
            column_split = anot[, column],
            height = unit(40, "cm")

        )
# Create the heatmap for methylation data
hm2 <- Heatmap(as.matrix(met),
            show_row_names = show_row_names,
            show_column_names = F,
            col = col.met,
            gap = unit(5, "mm"),
            cluster_rows = T,
            show_row_dend = F,
            show_column_dend = F,
            cluster_columns = F,
            clustering_distance_rows = "pearson",
            clustering_distance_columns = "pearson",
            heatmap_legend_param = list(
                title = "methylation Z score",
                title_position = "topcenter",
                legend_direction = "horizontal",
                title_gp = gpar(fontsize = 60),
                labels_gp = gpar(fontsize = 60,
                family = "Times", face = "bold"),
                legend_width = unit(30, "cm"),
                grid_height = unit(2.5, "cm")
            ),
            row_names_side = "left",
            row_names_gp = gpar(fontsize = 100),
            row_names_max_width = unit(names.size,'cm'),
            cluster_column_slices = FALSE,
            height = unit(40, "cm")
        )

# Combine the expression and methylation heatmaps vertically
hlist <- hm1 %v% hm2

# Draw and output the heatmap with custom legend positioning and padding
p <- draw(hlist, heatmap_legend_side = "top",
        annotation_legend_side = "top",
        padding = unit(c(20, 240, 20, 20), "mm"),
        ht_gap = unit(3, "cm"))

# Save the heatmap to a PDF file
pdf(paste0(
            outdir, "/Heatmap_anotated_",flag, ".pdf"
        ),
        width = size[["width"]], size[["height"]]
        )
        print(p)
        dev.off()

}


####################################################
# Function: anotation
####################################################
# Input:
#      dep: Data frame containing differential expression results of methylation data (must include `adj.P.Val` column)
#      anotation: Data frame containing methylation annotations to merge with `dep` (must have an "ID" column)
#      outdir: Directory path to save the output file
#      flag: String to append to the output file name
# Output:
#      dep: Merged and filtered data frame (significant annotated results)
#################################################
# Description:
#   This function merges annotation data with differential expression results 
#   for methylation data, filtering to retain only significant results (based 
#   on adjusted p-values). The results are then written to a specified output file.
#
####################################################
anotation <- function(dep, anotation, outdir, flag) {

    # Merge the annotation data with the differential expression data
    dep <- merge(anotation, dep, by.y = "row.names", by.x = "ID")
    
    # Rename the first column to "ID" for consistency
    colnames(dep)[1] <- "ID"
    
    # Sort the results by adjusted p-value (adj.P.Val)
    dep <- dep[order(dep$adj.P.Val),]
    
    # Filter the results to include only significant genes (adj.P.Val < 0.05)
    dep <- dep[dep$adj.P.Val < 0.05,]
    
    # Write the filtered results to a tab-delimited text file
    write.table(dep,
                sep = "\t",
                paste0(outdir, "/Significant_anotated_dep_", flag, ".txt"),
                row.names = FALSE,
                col.names = TRUE)
    
    # Return the filtered data frame
    return(dep)
}

####################################################
# Function: kruskal
####################################################
# Input:
#      base: data.frame with samples as rows and features as columns
#      features: column names to analyze (continuous variables)
#      group: column name of the feature used to stratify samples into groups (categorical variable)
#      outdir: directory path where the results should be written
#      flag: string identifier to append to output file names
# Output:
#      - A matrix with the p-value and adjusted p-value (using the Benjamini-Hochberg method) 
#        for each feature tested using the Kruskal-Wallis test.
#      - Results are saved as a tab-delimited file in the specified output directory.
# Description:
#      This function evaluates differences in continuous variables (features) among groups 
#      using the Kruskal-Wallis test. For each feature, a Kruskal-Wallis test is performed 
#      to assess if the median values differ significantly between the groups. The results 
#      are saved in a matrix and output to a text file.
####################################################

kruskal <- function(base,
group, 
features,
outdir,
flag){


 # Apply Kruskal-Wallis test for each feature
    # lapply iterates over each feature and applies the Kruskal-Wallis test
vf <- lapply(features, function(x) {
            
            
            # Filter the data to include only rows with complete cases for the current feature

            df <- base[complete.cases(base[,x]),]
            
            # Perform the Kruskal-Wallis test
            # This tests whether the distribution of the feature differs across groups
            k <- kruskal.test(as.formula(paste(x, " ~ ", group)),
                data = base)
            
            # Extract the p-value from the test result
            p <- k[[3]]

            return(c("p.value" = p))
        })

# Convert the list of p-values into a data frame
df <- as.data.frame(do.call(rbind,vf))

# Set the row names of the data frame to the names of the features
rownames(df) <- features

# Adjust p-values using the fdr method to control the false discovery rate
df$p.adjust <- p.adjust(df$p.value, method = "fdr")

# Save the results as a tab-delimited file
write.table(df,
paste0(outdir, "/Kruskal_table_", flag, ".txt"),
row.names = T,
col.names = NA,
sep = "\t")

# Return the data frame with p-values and adjusted p-values
return(df)
}

####################################################
# Function: anotation
#
# Description:
#      This function merges differential expression results with their corresponding 
#      annotations, filters the results based on adjusted p-values, and saves the 
#      annotated results to a specified output file.
#
# Input:
#      dep: data.frame containing differential expression results with rows corresponding
#           to features and columns including adjusted p-values.
#      anotation: data.frame containing annotations for the features, with a column named 'ID'.
#      outdir: Path to the directory where results will be saved.
#      flag: Suffix to append to the output file name for identification.
#
# Output:
#      A filtered data.frame containing significant annotated differential expression results.
####################################################

anotation <- function(dep, anotation, outdir, flag) {
  
  # Merge the differential expression results with their annotations using 'ID' as the key
  dep <- merge(anotation, dep, by.y = "row.names", by.x = "ID")
  
  # Rename the first column to 'ID' for clarity
  colnames(dep)[1] <- "ID"
  
  # Order the results by adjusted p-values
  dep <- dep[order(dep$adj.P.Val),]
  
  # Filter the results to retain only those with adjusted p-value < 0.05
  dep <- dep[dep$adj.P.Val < 0.05,]
  
  # Save the significant annotated results to a text file
  write.table(dep,
              sep = "\t",
              paste0(outdir, "/Significant_anotated_dep_", flag, ".txt"),
              row.names = FALSE,
              col.names = TRUE)
  
  # Return the filtered annotated results
  return(dep)
}


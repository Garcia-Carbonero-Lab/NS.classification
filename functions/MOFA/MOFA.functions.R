# Function: Mofa.model
# Input:
#   - input: Named list with the input matrices where samples are columns and features are rows.
#   - file: Path where the generated MOFA model will be saved.
#   - num_factors: Number of latent factors to infer in the MOFA model.
#   - convergence_mode: Convergence mode for model training, can be "slow" or "fast".
#   - groups: (Optional) Factor indicating sample groups for group mode. Must be in the same order as in the input matrix.
#
# Output:
#   - MOFA model saved in the specified file path.
#
# Description:
#   This function creates and trains a MOFA model using input data. It can work with or without sample grouping.
#   The model is configured with the specified number of factors and convergence mode, and it is trained using
#   MOFA's default options, which are modified for the specific requirements of the user.
#
# Example Usage:
#   Mofa.model(input = data.list, file = "output/Mofa.model.hdf5", num_factors = 8, convergence_mode = "slow", groups = sample_groups)
#


Mofa.model <- function(input,
                       file,
                       num_factors,
                       convergence_mode,
                       groups = NULL) {
    
    
    # Create MOFA model (with or without groups)
    
        if (is.null(groups) == T) {
        mofa1 <- create_mofa(input)
    } else {
        mofa1 <- create_mofa(input, groups = groups)
    }
    
    
        # Retrieve default model options
    ModelOptions <- get_default_model_options(mofa1)
    ModelOptions$num_factors <- num_factors
    
        # Retrieve default data options
    data_Options <- get_default_data_options(mofa1)
    
    # Retrieve and modify training options
    TrainOptions <- get_default_training_options(mofa1)
    TrainOptions$seed <- 57
    TrainOptions$convergence_mode <- convergence_mode

    # Prepare the MOFA model with the specified options
    mofa1 <- prepare_mofa(
        mofa1,
        data_options = data_Options,
        model_options = ModelOptions,
        training_options = TrainOptions
    )

       # Train the MOFA model and save it to the specified file
    mofa1 <- run_mofa(mofa1, outfile = file, use_basilisk = TRUE)
}

# Functions to create and analyze MOFA models

# function: anotate.metilation
# input:
#   - model: MOFA model (object of class MOFA)
#   - anotation: annotation matrix with information of EPIC CpG probes 
#
# output:
#   - MOFA model with genes annotated in CpG probe names (updated feature names)
#
# description:
#   This function updates the MOFA model by adding gene annotations to the 
#   feature names of the methylation probes. It retrieves the current feature 
#   names, matches them with the corresponding genes in the annotation matrix, 
#   and modifies the feature names to include gene symbols. The updated model 
#   is then returned for further analysis.

anotate.metilation <- function(model, anotation) {
    
    # Retrieve the current feature names from the MOFA model
    updated_features_names <- features_names(model)

    # Extract the annotation information for methylation features
    mt <- anotation[updated_features_names$Methylation, ]
    mt <- as.data.frame(mt)  # Convert to data.frame for easy manipulation

    # Change feature names by concatenating gene symbol with probe ID
    updated_features_names[["Methylation"]] <- as.vector(paste(mt$ID, mt$genes, sep = "_"))

    # Update the model with the new feature names
    features_names(model) <- updated_features_names

    # Return the updated MOFA model
    return(model)
}
# Function to calculate the mean variance explained by each factor in a MOFA model

# function: calculate.mean.variance
# input:
#   - model: a MOFA model object containing the data and results
#   - outdir: a string representing the output directory where results will be saved
#
# output:
#   - fve: a matrix with factors as rows and views as columns, indicating the mean 
#          variance explained by each factor in each view across all tissues.
#
# description:
#   This function calculates the mean variance explained by each factor in a 
#   given MOFA model for each view. It aggregates the variance explained data 
#   across all tissues and saves the results to a text file in the specified 
#   output directory.

calculate.mean.variance <- function(model, outdir) {
    # Calculate variance explained for each factor in the MOFA model
    ve <- calculate_variance_explained(model)

    # Combine the variance explained data into a single matrix
    ve <- Reduce(rbind, ve$r2_per_factor)

    # Get the names of the factors from the model
    factors <- colnames(model@expectations$Z[[1]])

    # Calculate the mean variance explained for each factor across all views
    fve <- lapply(factors, function(x) {
        apply(ve[rownames(ve) == x, ], 2, mean)  # Calculate mean for each view
    })

    # Combine the mean variances into a matrix
    fve <- do.call(rbind, fve)
    rownames(fve) <- factors  # Assign factor names as row names

    # Write the mean variance explained results to a text file
    write.table(fve,
        paste0(outdir, "/Variance_Mean_views.txt"),  # Construct output file path
        sep = "\t",  # Use tab as the separator
        row.names = TRUE,  # Include row names in the output
        col.names = NA  # Do not include column names in the output
    )
}

# Function to analyze associations between discrete variables and MOFA model factors

# function: discrete.features
# input:
#   - mod: a MOFA model object 
#   - p.th: a numeric value representing the p-value threshold for selecting variables in the Kruskal-Wallis test (default is 0.05)
#   - columns: a vector of column names representing the discrete variables to analyze
#   - outdir: a string representing the output directory where results will be saved
#   - flag: a string that will be appended to the names of output files to differentiate them
#   - colors: a list of named vectors indicating the colors for each variable
#   - select: a logical value (TRUE or FALSE); if TRUE, selects variables with Kruskal-Wallis p-values below the threshold
#   - group: a string indicating the column name used to stratify samples (default is "PRIMARY_TUMOR")
#   - levels: a vector of labels that define the order of the groups in the specified column
##   - we: logical, TRUE to calculate effect size (Wilcoxon Effect Size).

# output: 
#   - Generates various files and plots, including scatter plots and box plots, that display the results 
#     of the Kruskal-Wallis and Wilcoxon tests to analyze associations between discrete variables and model factors.
#
# description:
#   This function assesses the associations between discrete variables and the factors of a MOFA model. 
#   It performs Kruskal-Wallis tests to identify significant associations, generates plots to visualize these associations, 
#   and saves the statistical results to text files.

discrete.features <- function(mod,
    p.th = 0.05,
    columns,
    outdir,
    flag,
    colors,
    select = TRUE,
    group = "PRIMARY_TUMOR",
    levels = c("SMALL INTESTINE", "LUNG", "COLORECTAL", "GASTRIC", "PANCREAS"),
    we = TRUE) {
    # Create a directory for discrete variable outputs if it doesn't exist
    if (!dir.exists(paste0(outdir, "/discrete.variables/", flag))) {
        dir.create(paste0(outdir, "/discrete.variables/", flag))
    }

    # Obtain the factor matrix from the model
    Z <- Reduce(rbind, mod@expectations$Z)
    Z <- Z[mod@samples_metadata$sample, ]
    
    # Access sample metadata
    df <- mod@samples_metadata

    # Perform Kruskal-Wallis test and obtain results
    p.df <- kruskal(
        data = df, Z = Z, columns = columns,
        flag = flag, outdir = outdir, p.th = p.th
    )

    heatmap.kruskal (
        mod = mod, columns = columns,
        flag = flag, outdir = outdir, p.th = p.th
    )

    # Identify columns with significant p-values
    num <- apply(p.df, 1, sum)
    if (select) {
        columns <- rownames(p.df)[num != ncol(p.df)]
    } else {
        columns <- rownames(p.df)
    }

    # Remove grouping column from analysis
    columns <- columns[!columns %in% group]

    # Iterate over each column to create plots and save results
    for (column in columns) {
        print(column)

        # Create a directory for each variable if it doesn't exist
        if (!dir.exists(paste0(outdir, "/discrete.variables/", flag, "/", column))) {
            dir.create(paste0(outdir, "/discrete.variables/", flag, "/", column))
        }

        # Convert the column to a factor with specified levels
        mod@samples_metadata[, column] <- factor(mod@samples_metadata[, column],
            levels = names(colors[[column]])
        )

        # Select factors based on p-value thresholds
        if (select) {
            factors <- colnames(p.df)[p.df[column, ] < 1]
        } else {
            factors <- colnames(p.df)
        }

        # Plot all factors together
        plt <- plot_factor(mod,
            factors = as.numeric(gsub("Factor", "", factors)),
            color_by = column
        ) +
            scale_fill_manual(values = colors[[column]]) +
            theme(
                text = element_text(size = 18),
                axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
            ) +
            guides(color = guide_legend(override.aes = list(size = 15)))

        # Save the plot to a PDF file
        ggsave(
            plot = plt, filename = paste0(
                outdir, "/discrete.variables/",
                flag, "/", column, "/Factors_Scatterplot_",
                flag, "_", column, ".pdf"
            ),
            device = "pdf",
            dpi = 500, width = 15,
            height = 5
        )

        # Initialize lists to store statistical results
        stat_f <- list()
        stat_c <- list()
        
        # Analyze each factor separately
        for (factor in factors) {
            plt <- plot_factor(mod,
                factors = factor,
                color_by = column
            ) +
                scale_fill_manual(values = colors[[column]]) +
                theme(
                    text = element_text(size = 18),
                    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
                ) +
                guides(color = guide_legend(override.aes = list(size = 15)))

            # Save the factor-specific plot to a PDF file
            ggsave(
                plot = plt, filename = paste0(
                    outdir, "/discrete.variables/",
                    flag, "/", column, "/Factors_Scatterplot_",
                    flag, "_", column, "_", factor, ".pdf"
                ),
                device = "pdf",
                dpi = 500, width = 10,
                height = 5
            )

            # Generate boxplots stratified by group and obtain Wilcoxon results
            stat_df <- wrap.boxplot(
                mod = mod, column = column, factor = factor,
                flag = flag, colors = colors, outdir = outdir,
                group = group, levels = levels, we = we
            )

            stat_f[[factor]] <- as.data.frame(stat_df)

            # Generate boxplots using all samples together and obtain Wilcoxon results
            stat_complete <- complete.boxplot(
                mod = mod, column = column, factor = factor,
                flag = flag, colors = colors, outdir = outdir
            )

            stat_c[[factor]] <- as.data.frame(stat_complete)
        }

        # Combine results from all selected variables into data frames
        stat_f <- do.call(rbind, stat_f)
        stat_c <- do.call(rbind, stat_c)

        # Adjust p-values using the false discovery rate method
        stat_f$p.adj.all <- p.adjust(stat_f$p, method = "fdr")
        stat_c$p.adj.all <- p.adjust(stat_c$p, method = "fdr")

        # Save the results of the Wilcoxon tests for the factors
        write.table(stat_f,
            paste0(
                outdir, "/discrete.variables/",
                flag, "/", column, "/Factor_Wilcox_", flag, "_", column, ".txt"
            ),
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE
        )

        # Save the results of the Wilcoxon tests for complete analysis
        write.table(stat_c,
            paste0(
                outdir, "/discrete.variables/",
                flag, "/", column, "/Factor_Wilcox_Complete_", flag,
                "_", column, ".txt"
            ),
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE
        )
    }
}
# Function to apply the Kruskal-Wallis test between factors and discrete variables

# function: kruskal
# input:
#   - data: a data.frame containing discrete variables and their corresponding samples
#   - Z: a matrix with samples in rows and factors in columns
#   - columns: a vector of column names representing the discrete variables to analyze
#   - flag: a string to label the output files
#   - outdir: a string indicating the output directory for saving results
#   - p.th: a numeric value representing the p-value threshold for filtering significant results
#
# output:
#   - Returns a data frame of adjusted p-values for each discrete variable against factors.
#   - Generates a text file with results from the Kruskal-Wallis tests.

kruskal <- function(data, Z, columns, flag, outdir, p.th) {
    # Calculate the Kruskal-Wallis test for each discrete variable
    ano_df <- lapply(columns, function(x) {
        vf <- lapply(colnames(Z), function(y) {
            # Perform the Kruskal-Wallis test
            k <- kruskal.test(as.formula(paste(y, " ~ ", x)), data = data)
            f <- k[[1]]  # Extract Chi-squared statistic
            p <- k[[3]]  # Extract p-value

            return(c("Chi-squared" = f, "p.value" = p))
        })

        names(vf) <- colnames(Z)  # Name the results by factor names
        vf <- as.data.frame(do.call(rbind, vf))  # Combine results into a data frame
        colnames(vf) <- paste(x, colnames(vf))  # Rename columns to include the variable name
        return(as.data.frame(vf))  # Return the result as a data frame
    })

    # Combine results from all discrete variables into a single data frame
    ano_df <- as.data.frame(t(do.call(cbind, ano_df)))

    # Extract p-values and adjust for multiple comparisons using the FDR method
    p.df <- ano_df[grep("p.value", rownames(ano_df)), ]
    for (col in colnames(p.df)) {
        p.df[, col] <- p.adjust(p.df[, col], method = "fdr")
    }

    # Prepare adjusted p-value data frame
    p.adj.df <- p.df
    rownames(p.adj.df) <- gsub("p.value", "p.adj", rownames(p.adj.df))  # Rename rows for clarity
    ano_df <- rbind(ano_df, p.adj.df)  # Combine adjusted p-values with the results
    ano_df <- t(ano_df[order(rownames(ano_df)), ])  # Transpose and order results

    # Create output directories if they do not exist
    if (!dir.exists(paste0(outdir, "/discrete.variables/"))) {
        dir.create(paste0(outdir, "/discrete.variables/"))
    }

    if (!dir.exists(paste0(outdir, "/discrete.variables/", flag))) {
        dir.create(paste0(outdir, "/discrete.variables/", flag))
    }

    # Write the results to a text file
    write.table(ano_df, paste0(
        outdir,
        "/", flag, "_Kruskal.txt"
    ),
    sep = "\t",
    row.names = TRUE,
    col.names = NA
    )

    # Prepare p-values for visualization
    p.df[p.df > p.th] <- 1  # Set non-significant p-values to 1
    p.df <- p.df[apply(p.df, 1, sum) != nrow(p.df), ]  # Remove rows with all non-significant results
    rownames(p.df) <- gsub(" p.value", "", rownames(p.df))  # Clean up row names

    return(p.df)  # Return the adjusted p-value data frame
}

# Function to generate a heatmap from Kruskal-Wallis test results

# function: heatmap.kruskal
# input:
#   - mod: the MOFA model
#   - columns: names of the discrete variables to plot
#   - outdir: output directory for saving results
#   - flag: a string to label the output files
#   - p.th: numeric value for p-value threshold (default is 0.05)
#
# output: 
#   - Generates a heatmap indicating the p-values of features associated with factors.

# description: Function to generate the heatmap with the results of the Kruskal-Wallis test

heatmap.kruskal <- function(mod, columns, flag, outdir, p.th = 0.05) {
    # Combine expectations from the MOFA model to form a factor matrix
    Z <- Reduce(rbind, mod@expectations$Z)
    Z <- Z[mod@samples_metadata$sample, ]  # Align samples with metadata

    df <- mod@samples_metadata  # Extract the sample metadata

    # Perform the Kruskal-Wallis test
    p.df <- kruskal(
        data = df, Z = Z, columns = columns,
        flag = flag, outdir = outdir, p.th = p.th
    )
    
    # Define a color gradient for the heatmap
    col <- colorRamp2(
        c(max(-log10(t(p.df))), 0),
        c("red", "grey")
    )

    # Generate the heatmap
    hm1 <- Heatmap(as.matrix(-log10(t(p.df))),  # Use -log10 of p-values for visualization
        show_row_names = TRUE,
        show_column_names = TRUE,
        col = col,
        cluster_rows = FALSE,  # Do not cluster rows
        cluster_columns = TRUE,  # Cluster columns
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        heatmap_legend_param = list(
            title = "-log10(q-value)",  # Legend title
            title_position = "topcenter",
            legend_direction = "horizontal",
            title_gp = gpar(fontsize = 60),
            labels_gp = gpar(fontsize = 60, family = "Times", face = "bold"),
            legend_width = unit(20, "cm"),
            grid_height = unit(2, "cm")
        ),
        row_names_side = "right",
        row_names_gp = gpar(fontsize = 60),
        column_names_gp = gpar(fontsize = 60),
        column_names_rot = 45,
        rect_gp = gpar(col = "black", lwd = 2)  # Border properties
    )

    # Draw the heatmap
    p <- draw(hm1, heatmap_legend_side = "bottom",
        padding = unit(c(2, 10, 2, 20), "mm"))

    # Save the heatmap as a PDF
    pdf(paste0(
        outdir, "/Heatmap_anotated_", flag, ".pdf"
    ),
    width = 15, height = 20
    )
    print(p)  # Print the heatmap
    dev.off()  # Close the PDF device

    return(p)  # Return the heatmap object
}

# Function to create boxplots comparing the levels of a specified factor across different groups.

# function: wrap.boxplot
# input:
#   - mod: MOFA model or clinical data containing factor information.
#   - column: column name for the feature to analyze.
#   - factor: name of the factor to analyze.
#   - flag: a string to label the output files.
#   - colors: a named vector indicating colors for the labels.
#   - outdir: path to save output files.
#   - is.model: logical, TRUE if mod is a MOFA model, FALSE if it is a clinical data frame.
#   - group: column name to stratify samples.
#   - levels: desired order for the group levels.
#   - we: logical, TRUE to calculate effect size (Wilcoxon Effect Size).
#
# output:
#   - Returns a data frame with results from the Wilcoxon tests, including adjusted p-values.

# description: Function to create boxplots comparing factor levels between groups defined by a specified column feature.

wrap.boxplot <- function(mod, column, factor, flag, colors, outdir, is.model = TRUE, group, levels, we = TRUE) {

    # Obtain clinical data frame from MOFA model or use provided clinical data
    if (is.model) {
        cli.df <- mod@samples_metadata
    } else {
        cli.df <- mod
    }

    # Clean and prepare the data
    cli.df[, column] <- gsub(" ", "_", cli.df[, column])  # Replace spaces with underscores
    cli.df[, group] <- gsub(" ", "_", cli.df[, group])  # Same for group column
    levels <- gsub(" ", "_", levels)  # Clean levels
    cli.df[, group] <- factor(cli.df[, group], levels = levels)  # Set factor levels

    # Eliminate samples with missing values
    df <- table(cli.df[, group], cli.df[, column])  # Create contingency table
    l <- rownames(df)[df[, 1] == 0 | df[, 2] == 0]  # Identify groups with no samples
    if (length(l) > 0) {
        cli.df <- cli.df[cli.df[, group] != l, ]  # Remove these groups from the data frame
    }

    # Convert specified column to factor and drop NAs
    cli.df[, column] <- as.factor(cli.df[, column])
    cli.df <- cli.df[complete.cases(cli.df[, column]), ]  # Remove rows with NAs
    cli.df[, column] <- droplevels(cli.df[, column])  # Drop unused factor levels

    # Obtain Wilcoxon test results
    stat_df <- compare_means(as.formula(paste0(factor, " ~ ", column)), group.by = group, data = cli.df, p.adjust.method = "fdr")

    # Calculate median values for groups
    median <- cli.df[, c(column, group, factor)] %>% 
        group_by_at(vars(one_of(c(group, column)))) %>%
        get_summary_stats(type = "median")

    # Save median values to a TSV file
    write.table(median,
        paste0(outdir, "/discrete.variables/", flag, "/", column, "/Medians", flag, "_", column, "_", factor, ".tsv"),
        sep = "\t",
        col.names = TRUE)

    if (we) {  # Calculate effect size if requested
        # Calculate Wilcoxon Effect Size
        we <- cli.df %>% 
            group_by_at(vars(group)) %>%
            wilcox_effsize(as.formula(paste0(factor, " ~ ", column)))
        
        we$group <- gsub(" ", "_", we$group)  # Clean up group names
        we$names <- paste(we$.y., we$group1, we$group2, we$group, sep = "_")  # Combine group names for merging
        
        # Prepare Wilcoxon statistics for merging
        stat_df$names <- paste(stat_df$.y., stat_df$group1, stat_df$group2, stat_df[, group][[1]], sep = "_")
        
        # Merge effect size data with Wilcoxon results
        stat_df <- merge(stat_df, we[, c(4:9)], by = "names")
    } else {
        stat_df$names <- paste(stat_df$.y., stat_df$group1, stat_df$group2, stat_df[, group][[1]], sep = "_")
    }

    # Generate boxplot
    plt <- ggboxplot(cli.df,
        y = factor,
        x = group, fill = column,
        width = 0.6
    ) +
        scale_fill_manual(values = colors[[column]]) +  # Apply specified colors
        theme_classic() +  # Use classic theme

        # Add p-value annotations to the plot
        geom_pwc(aes(group = .data[[column]]),
            method = "wilcox.test", label = "p.adj.signif",
            hide.ns = TRUE,
            method.args = list(p.adjust.method = "BH", exact = FALSE)
        ) +

        # Separate boxplots by group column
        facet_wrap(as.formula(paste0("~ ", group)), drop = TRUE, scales = "free_x", nrow = 1) +

        # Customize theme sizes and appearances
        theme(
            text = element_text(size = 40),
            axis.text.x = element_blank(),  # Hide x-axis text
            legend.key.size = unit(2, "cm"),
            legend.text = element_text(size = 40),
            strip.background = element_rect(color = "black", size = 3, linetype = "solid")  # Customize strip appearance
        )

    # Customize the size of the p-value annotations
    plt$layers[[2]]$aes_params$size <- 1
    plt$layers[[2]]$aes_params$label.size <- 26

    # Save the boxplot as a PDF
    ggsave(
        plot = plt, 
        filename = paste0(outdir, "/discrete.variables/", flag, "/", column, "/Factors_Boxplot_", flag, "_", column, "_", factor, ".pdf"),
        device = "pdf",
        dpi = 500, 
        width = 20,
        height = 15
    )

    return(stat_df)  # Return the statistics data frame
}

# Function: complete.boxplot
# Purpose: Create boxplots that compare the levels of a factor between 
#          the labels of a column, with the option to use a MOFA model or clinical data.

# Input parameters:
#   mod: MOFA model (or clinical data frame with factor information).
#   column: The column name to analyze (independent variable).
#   factor: The name of the factor to analyze (dependent variable).
#   flag: A name to add to the output files for identification.
#   colors: List of named vectors indicating the colors for the labels.
#   outdir: Path to the directory for saving output files.
#   is.model: Logical, if TRUE, `mod` is a MOFA model, if FALSE, `mod` is a data frame.

# Output: A boxplot that compares the levels of the factor between labels of the column, 
#         including results from Wilcoxon test and summary statistics.

complete.boxplot <- function(mod, column, factor, flag, colors, outdir, is.model = TRUE) {
    
    # If `mod` is a MOFA model, extract the metadata, else treat `mod` as a clinical data frame
    if (is.model == TRUE) {
        cli.df <- mod@samples_metadata
    } else {
        cli.df <- mod
    }

    # Replace spaces with dots in column values and color names for compatibility
    cli.df[, column] <- gsub(" ", ".", cli.df[, column])
    names(colors) <- gsub(" ", ".", names(colors))

    # Set the factor levels in `column` according to the provided colors vector
    cli.df[, column] <- factor(cli.df[, column], levels = names(colors[[column]]))

    # Convert the column to a factor, remove rows with missing values, and drop unused levels
    cli.df[, column] <- as.factor(cli.df[, column])
    cli.df <- cli.df[complete.cases(cli.df[, column]), ]  # Remove missing values
    cli.df[, column] <- droplevels(cli.df[, column])      # Drop unused factor levels

    # Perform pairwise comparisons between the groups in `column` using Wilcoxon test
    stat_df <- compare_means(as.formula(paste0(factor, " ~ ", column)),
                             data = cli.df, p.adjust.method = "BH")  # Adjust p-values

    # Calculate the median summary statistics for each group in `column`
    median <- cli.df[, c(column, factor)] %>%
        group_by_at(vars(column)) %>%
        get_summary_stats(type = "median")  # Get summary statistics (median)

    # Save the median statistics to a file
    write.table(median,
                paste0(outdir, "/discrete.variables/", flag, "/", column, 
                       "/Medians_Complete", flag, "_", column, "_", factor, ".tsv"),
                sep = "\t", col.names = TRUE)

    # Calculate effect sizes for the Wilcoxon test using `wilcox_effsize`
    we <- wilcox_effsize(data = cli.df, as.formula(paste0(factor, " ~ ", column)))
    
    # Create a combined name for each comparison group to merge results later
    we$names <- paste(we$.y., we$group1, we$group2, sep = "_")
    stat_df$names <- paste(stat_df$.y., stat_df$group1, stat_df$group2, sep = "_")
    
    # Merge the Wilcoxon results (p-values) with the effect sizes (Cohen's d)
    stat_df <- merge(stat_df, we[, c(4:8)], by = "names")

    # Create the boxplot using `ggboxplot`
    plt <- ggboxplot(cli.df,
                     y = factor,      # Dependent variable (y-axis)
                     x = column,      # Independent variable (x-axis)
                     fill = column,   # Fill colors based on the column levels
                     width = 0.6      # Width of the boxes
    ) +
        scale_fill_manual(values = colors[[column]]) +  # Set fill colors
        theme_classic() +  # Use the classic theme for the plot

        # Add pairwise comparisons (p-values) from the Wilcoxon test
        geom_pwc(
            aes(group = .data[[column]]),  # Group the comparisons by column
            method = "wilcox.test",        # Use Wilcoxon test
            label = "p.adj.signif",        # Add significance labels (adjusted p-values)
            hide.ns = TRUE,                # Hide non-significant results
            method.args = list(p.adjust.method = "BH", exact = FALSE)  # Adjust p-values
        ) +

        # Customize theme settings for text, axis, and legend appearance
        theme(
            text = element_text(size = 40),               # Text size
            axis.text.x = element_blank(),                # Remove x-axis text
            legend.key.size = unit(2, "cm"),              # Set legend key size
            legend.text = element_text(size = 40)         # Set legend text size
        )

    # Set aesthetics for the significance labels
    plt$layers[[2]]$aes_params$size <- 1              # Line thickness for p-value annotations
    plt$layers[[2]]$aes_params$label.size <- 26       # Label text size for significance annotations

    # Save the plot to a PDF file
    ggsave(
        plot = plt, filename = paste0(
            outdir, "/discrete.variables/", flag, "/", column, 
            "/Factors_Boxplot_Complete_", flag, "_", column, "_", factor, ".pdf"
        ),
        device = "pdf", dpi = 500, width = 20, height = 15
    )

    # Return the data frame with statistical results
    return(stat_df)
}
####################################################
# Function: continous.features
####################################################
# Input:
#      mod: MOFA model object containing the expectations matrix and metadata
#      p.th: p-value threshold for selecting variables in Correlation test
#      columns: column names of the variables to analyze (continuous variables)
#      outdir: output directory path to write results
#      flag: string identifier to append to output file names
#      group: column name used to stratify samples (default is "PRIMARY_TUMOR")
#      levels: vector specifying label order for the group column
# Output:
#      Files and plots that show results of Correlation tests 
#      analyzing the associations between continuous variables and model factors
# Description:
#      This function studies associations between continuous variables 
#      (given by columns) and the latent factors from the MOFA model. 
#      It performs statistical tests and creates plots to visualize the 
#      relationships between the variables and the MOFA factors.
####################################################

continous.features <- function(mod,
p.th = 0.05,
columns,
outdir,
flag,
group = "PRIMARY_TUMOR",
levels = c("COLORECTAL", "STOMACH", "LUNG", "PANCREAS", "SMALL INTESTINE")) {


    # Check if output directory exists, if not, create it
if (dir.exists(paste0(outdir, "/continous.variables/", flag)) == F) {

    dir.create(paste0(outdir, "/continous.variables/", flag))

}

# Extract latent factors (Z matrix) from MOFA model

    Z <- Reduce(rbind, mod@expectations$Z)
    Z <- Z[mod@samples_metadata$sample, ]

# Create a data frame from the model's sample metadata
    df <- mod@samples_metadata

# Perform correlation tests between metadata columns and MOFA factors
# Returns a data frame with correlation coefficients and p-values
    c.df <- correlation.test(
        data = df, Z = Z, columns = columns,
        flag = flag, outdir = outdir, p.th = p.th
    )

    # Remove '.rho' suffix from row names for clarity
    rownames(c.df) <- gsub("\\.rho", "", rownames(c.df))

    # Get the factor names from the Z matrix columns
    factors <- colnames(Z)

    
    
    
    # Identify columns that have significant associations with any factor


    
    num <- apply(c.df,1,sum)
    columns <- rownames(c.df)[num != 0]

    # Loop through each column that has significant associations
    for (column in columns) {

        # Check if output sub-directory for the current column exists, create if not
        if (dir.exists(paste0(
            outdir, "/continous.variables/",
            flag, "/", column
        )) == F) {
            dir.create(paste0(
                outdir, "/continous.variables/",
                flag, "/", column
            ))
        }


        ####################################################
        # Plot all factors together (scatterplots colored by column value)
        ####################################################

        # Create a scatterplot for all factors, colored by the current column
        plt <- plot_factor(mod,
            factors = as.numeric(gsub("Factor", "", factors)),
            color_by = column
        ) +
            theme(
                text = element_text(size = 18),
                axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
            )

        # Save the scatterplot for all factors in a single PDF file
        ggsave(
            plot = plt, filename = paste0(
                outdir, "/continous.variables/",
                flag, "/", column, "/Factors_Scatterplot_",
                flag, "_", column, ".pdf"
            ),
            device = "pdf",
            dpi = 500, width = 15,
            height = 5
        )



        ####################################################
        # Plot each factor separately (only significant ones)
        ####################################################


        # Loop through each factor and plot it separately if it has significant correlation
        for (factor in factors) {

        if (abs(c.df[column, factor]) > 0){
        # Create a scatterplot for the current factor, colored by the current column
            plt <- plot_factor(mod,
                factors = factor,
                color_by = column
            ) +
                theme(
                    text = element_text(size = 18),
                    axis.text.x = element_text(angle = 90,
                    vjust = 0.5, hjust = 1)
                ) 
        # Save the scatterplot for the current factor
            ggsave(
                plot = plt, filename = paste0(
                    outdir, "/continous.variables/",
                    flag, "/", column, "/Factors_Scatterplot_",
                    flag, "_", column, "_", factor, ".pdf"
                ),
                device = "pdf",
                dpi = 500, width = 10,
                height = 5
            )


          ####################################################
            # Plot significant variables (correlation plot)
          ####################################################

            # Plot the correlation between the current column and factor

            correlation.plot(
                mod = mod, column = column, factor = factor,
                flag = flag, outdir = outdir, group = group, levels = levels
            )

        }

        }

    }
}

####################################################
# Function: correlation.plot
####################################################
# Input:
#      mod: MOFA model object or clinical data with factor information
#      factor: the name of the factor to analyze (MOFA latent factor)
#      column: the name of the column/variable to analyze (continuous variable)
#      outdir: path to write the output files (plots)
#      flag: string identifier to append to the output file names
#      group: column name used to stratify samples (factor levels)
#      levels: vector specifying label order for the group column
# Output:
#      Correlation plots between factors and variables for all samples,
#      and stratified by the specified group column. Two plots are generated:
#      one for all samples and one faceted by groups.
# Description:
#      This function generates correlation plots between a selected latent
#      factor and a variable (from the column). It produces two types of plots:
#      (1) a scatter plot showing correlation across all samples, and
#      (2) a scatter plot with samples stratified by the group column using facets.
####################################################


correlation.plot <- function(mod, column, factor,
flag, outdir, group, levels) {

    # Extract clinical data or sample metadata from the MOFA model
 cli.df <- mod@samples_metadata

    # Remove rows with missing values in the specified column
    cli.df <- cli.df[complete.cases(cli.df[, column]), ]

    # Convert the group column to a factor 
    cli.df[,group] <- as.factor(cli.df[,group])

####################################################
    # Plot 1: Correlation plot for all samples
    ####################################################
    # Create scatter plot of the column vs. factor with linear regression
    # and correlation statistics

plt1 <- ggplot(cli.df, aes_string(x = column, y = factor, group = 1)) +
        theme_classic() +  # Clean theme
        geom_point(size = 3) +  # Add points to scatter plot
        stat_cor(size = 10) +  # Display correlation coefficient and p-value
        stat_smooth(method = "lm", color = "black") +  # Add linear regression line
        theme(
            text = element_text(size = 26),
            legend.key.size = unit(2, "cm"),
            legend.text = element_text(size = 26)
        )

    ggsave(
        plot = plt1, filename = paste0(
            outdir, "/continous.variables/",
            flag, "/", column, "/Factors_Corr_plot_", flag,
            "_", column, "_", factor, ".pdf"
        ),
        device = "pdf",
        dpi = 500, width = 17,
        height = 12
    )



    # Save the first plot (correlation plot for all samples)
    ggsave(
        plot = plt1,
        filename = paste0(outdir, "/continous.variables/", flag, "/", column, "/Factors_Corr_plot_", flag, "_", column, "_", factor, ".pdf"),
        device = "pdf",
        dpi = 500,
        width = 17,
        height = 12
    )

    ####################################################
    # Plot 2: Correlation plot faceted by group
    ####################################################
    # Create scatter plot of the column vs. factor, faceted by group column
    plt2 <- ggplot(cli.df, aes_string(x = column, y = factor)) +
        theme_classic() +  # Clean theme
        stat_cor(size = 10) +  # Display correlation coefficient and p-value
        stat_smooth(method = "lm") +  # Add linear regression line
        geom_point(size = 3) +  # Add points to scatter plot
        facet_wrap(as.formula(paste0("~ ", group)),  # Create facet for each group level
                   drop = TRUE, scales = "free_x", nrow = 1) +  # Free x-axis scale
        theme(
            text = element_text(size = 26),
            legend.key.size = unit(2, "cm"),
            legend.text = element_text(size = 26)
        )

    # Save the second plot (faceted correlation plot by group)
    ggsave(
        plot = plt2,
        filename = paste0(outdir, "/continous.variables/", flag, "/", column, "/Factors_Corr_plot_wrapper_", flag, "_", column, "_", factor, ".pdf"),
        device = "pdf",
        dpi = 500,
        width = 22,
        height = 12
    )
}

####################################################
# Function: correlation.test
####################################################
# Input:
#      data: clinical data frame, containing both metadata and continuous variables
#      Z: matrix with samples as rows and latent factors as columns (from MOFA)
#      columns: vector of column names representing continuous variables to analyze
#      flag: string identifier to append to output file names
#      outdir: path to write output files (results and plots)
#      p.th: p-value threshold to select variables in Spearman correlation
# Output:
#      - Data frame with Spearman correlation estimates and adjusted p-values
#      - Heatmap of significant correlations between factors and continuous features
#      - Text file containing full correlation results
# Description:
#      This function performs Spearman correlation tests between MOFA factors
#      and continuous variables. It adjusts p-values using the Benjamini-Hochberg
#      method and generates a heatmap of significant correlations, saving the
#      results to output files.
####################################################

correlation.test <- function(data, Z, columns, flag, outdir, p.th) {

    # Loop over each column (continuous variable) and calculate correlation with factors
    ano_df <- lapply(columns, function(x) {

        # Subset data to remove rows with missing values for the current column
        cdata <- data[complete.cases(data[, x]),]

        # For each factor (columns of Z), calculate Spearman correlation with the current variable
        vf <- lapply(colnames(Z), function(y) {

            # Perform Spearman correlation test between variable 'x' and factor 'y'
            corr <- cor.test(cdata[, x], cdata[, y], method = "spearman")
            f <- corr$estimate[[1]]  # Correlation coefficient (Spearman's rho)
            p <- corr$p.value  # P-value of the correlation test

            # Return correlation coefficient and p-value as a named vector
            return(c("estimate" = f, "p.value" = p))
        })

        # Name the results after the factors
        names(vf) <- colnames(Z)

        # Combine the results into a data frame for the current variable
        vf <- as.data.frame(do.call(rbind, vf))
        colnames(vf) <- paste(x, colnames(vf))  # Add variable name to column headers

        return(as.data.frame(vf))  # Return correlation data frame for current variable
    })

    # Combine the correlation results for all variables
    ano_df <- as.data.frame(t(do.call(cbind, ano_df)))

    ####################################################
    # Adjust p-values using Benjamini-Hochberg correction
    ####################################################
    # Extract only the p-values from the results
    p.df <- ano_df[grep("p.value", rownames(ano_df)), ]

    # Adjust p-values for multiple testing using the Benjamini-Hochberg method
    for (col in colnames(p.df)) {
        p.df[, col] <- p.adjust(p.df[, col], method = "BH")
    }

    # Replace original p-values with adjusted p-values in the main result data frame
    ano_df[grep("p.value", rownames(ano_df)), ] <- p.df

    ####################################################
    # Save the full correlation results to a file
    ####################################################
    write.table(
        ano_df, paste0(outdir, "/continous.variables/", flag, "/Correlation_Features_Factors_", flag, ".txt"),
        sep = "\t", row.names = TRUE, col.names = NA
    )

    ####################################################
    # Filter results based on the p-value threshold and create a heatmap
    ####################################################
    # Extract the correlation coefficients (estimates)
    c.df <- ano_df[grep("estimate", rownames(ano_df)), ]
    rownames(c.df) <- gsub("\\.cor", "", rownames(c.df))  # Clean row names

    # Set correlation values to 0 where the adjusted p-value exceeds the threshold
    c.df[p.df > p.th] <- 0

    # Define color breaks for the heatmap (from -1 to 1 with small intervals)
    breaksList <- seq(from = -1, to = 1, by = 0.01)
    col <- colorRampPalette(c("blue", "lightgrey", "red"))(length(breaksList))

    # Clean row names by removing "estimate" from them
    rownames(c.df) <- gsub(" estimate", "", rownames(c.df))

    # Generate the heatmap for the correlation results
    pheatmap::pheatmap(
        t(c.df),  # Transpose the matrix for proper orientation
        cluster_rows = FALSE,  # Do not cluster rows
        color = col,  # Use the color palette defined above
        breaks = breaksList,  # Apply the color breaks
        filename = paste0(outdir, "/continous.variables/", flag, "/Heatmap_Correlation_", flag, ".pdf")
    )

    # Return the final data frame with correlation estimates
    return(c.df)
}


####################################################
# Function: immune_lm
####################################################
# Input:
#      data: matrix with samples as rows and immune cells as columns (immune cell infiltration data)
#      mod: MOFA model (contains factors and metadata)
#      flag: string identifier to append to output file names
#      outdir: path to write output files (results and plots)
# Output:
#      - Linear model results associating factors with immune infiltration
#      - p-values and adjusted p-values (FDR) for each immune cell type and factor
#      - R-squared values and estimated coefficients from the linear models
#      - Heatmap of significant associations (log-transformed p-values)
# Description:
#      This function fits a linear model to study associations between MOFA factors
#      and immune cell infiltration data. It adjusts p-values using the FDR method,
#      calculates R-squared values, and generates a heatmap showing the strength of
#      these associations.
####################################################

immune_lm <- function(data, mod, flag, outdir) {

    ####################################################
    # Prepare data: Merge MOFA factor matrix (Z) with immune data
    ####################################################
    # Extract latent factor matrix (Z) from the MOFA model
    Z <- obtain_Z(mod)
    
    # Merge MOFA factors (Z) with immune cell infiltration data by row (samples)
    df <- merge(Z, data, by = "row.names")
    
    # Initialize lists to store p-values, R-squared values, and regression estimates
    p.val <- list()
    r2 <- list()
    estimate <- list()

    ####################################################
    # Perform linear models for each MOFA factor
    ####################################################
    for (f in colnames(Z)) {
        # Create a linear model for each factor f, using all immune cells as independent variables
        lmx <- lm(as.formula(paste0(f, " ~ ", paste0(colnames(data), collapse = " + "))), data = df)
        
        # Extract summary statistics of the linear model
        sum <- summary(lmx)
        
        # Store p-values, R-squared values, and regression estimates for each factor
        p.val[[f]] <- sum$coefficients[,"Pr(>|t|)"]
        r2[[f]] <- sum$r.squared
        estimate[[f]] <- sum$coefficients[,"Estimate"]
    }

    ####################################################
    # Save linear model results to output files
    ####################################################
    # Combine p-values into a matrix
    p.val <- do.call(cbind, p.val)

    # Adjust p-values using FDR (False Discovery Rate) method
    p.adjust <- apply(p.val, 2, p.adjust, method = "fdr")

    # Write adjusted p-values to a file
    write.table(p.adjust, paste0(outdir, "/continous.variables/", flag, "/Linear_model_p.value.adjust.txt"),
                sep = "\t", row.names = TRUE, col.names = NA)
    
    # Combine R-squared values and regression estimates
    r2 <- do.call(rbind, r2)
    estimate <- do.call(rbind, estimate)

    # Assign column name "R2" to the R-squared matrix
    colnames(r2) <- "R2"

    # Write unadjusted p-values, R-squared values, and regression estimates to files
    write.table(p.val, paste0(outdir, "/continous.variables/", flag, "/Linear_model_p.value.txt"),
                sep = "\t", row.names = TRUE, col.names = NA)
    
    write.table(r2, paste0(outdir, "/continous.variables/", flag, "/Linear_model_r2.txt"),
                sep = "\t", row.names = TRUE, col.names = NA)

    write.table(estimate, paste0(outdir, "/continous.variables/", flag, "/Linear_model_estimate.txt"),
                sep = "\t", row.names = TRUE, col.names = NA)

    ####################################################
    # Create heatmap of -log(p-values) for significant associations
    ####################################################
    # Replace non-significant p-values (p > 0.05) with 1
    p.val[p.val > 0.05] <- 1
    # Convert p-values to -log(p-value) for heatmap visualization
    p.val <- -log(p.val)

    # Remove the intercept row from the p-value matrix for the heatmap
    p.val <- p.val[-1, ]

    ####################################################
    # Prepare heatmap annotations and parameters
    ####################################################
    # Define color scale for the R-squared values
    col <- list("R2" = colorRamp2(c(1, 0), c("red", "white")))

    # Define legend parameters for the R-squared annotation in the heatmap
    anot_legend_param <- list(
        "R2" = list(
            title_gp = gpar(fontsize = 40),
            labels_gp = gpar(fontsize = 40),
            legend_height = unit(8, "cm"),
            grid_width = unit(2, "cm"),
            title_position = "leftcenter-rot"
        )
    )

    # Create heatmap annotation object for R-squared values
    ha_column <- HeatmapAnnotation(
        df = as.data.frame(r2),
        show_legend = TRUE,
        col = col,
        annotation_name_gp = gpar(fontsize = 40),
        border = unit(40, "cm"),
        simple_anno_size = unit(2, "cm"),
        gp = gpar(fontsize = 40),
        annotation_legend_param = anot_legend_param
    )

    ####################################################
    # Generate the heatmap with -log(p-values)
    ####################################################
    hm <- Heatmap(as.matrix(p.val),
                  top_annotation = ha_column,
                  show_row_names = TRUE,
                  show_column_names = TRUE,
                  col = colorRamp2(c(max(p.val), 0), c("red", "white")),
                  name = "Heatmap",
                  gap = unit(5, "mm"),
                  cluster_rows = FALSE,
                  show_row_dend = FALSE,
                  show_column_dend = FALSE,
                  cluster_columns = FALSE,
                  clustering_distance_rows = "spearman",
                  heatmap_legend_param = list(
                      title = "-log(p.value)",
                      title_position = "topcenter",
                      legend_direction = "horizontal",
                      title_gp = gpar(fontsize = 40),
                      labels_gp = gpar(fontsize = 40, family = "Times", face = "bold"),
                      legend_width = unit(20, "cm"),
                      grid_height = unit(2, "cm")
                  ),
                  row_names_side = "left",
                  row_names_gp = gpar(fontsize = 36),
                  column_names_gp = gpar(fontsize = 36),
                  rect_gp = gpar(col = "black", lwd = 2)
    )

    # Draw the heatmap and arrange the legend
    p <- draw(hm, heatmap_legend_side = "bottom", padding = unit(c(2, 120, 2, 2), "mm"))

    ####################################################
    # Save the heatmap to a PDF file
    ####################################################
    pdf(paste0(outdir, "/continous.variables/", flag, "/immune.linear.model.pdf"), width = 20, height = 20)
    print(p)
    dev.off()
}
####################################################
# Function: surv.model
####################################################
# Input:
#      mod: MOFA model (contains factors and metadata for samples)
#      time: string specifying the column name containing time information (e.g., survival time)
#      event: string specifying the column name containing event information (e.g., censoring status)
#      flag: string identifier to append to output file names (used to differentiate runs or datasets)
#      outdir: path to write output files (location where results and plots will be saved)
# Output:
#      - A table containing the results of the Cox proportional hazards model:
#          * factor names
#          * p-values for each factor
#          * hazard ratios (exponentiated coefficients)
#          * lower and upper 95% confidence intervals for hazard ratios
#      - A forest plot visualizing the hazard ratios with confidence intervals
# Description:
#      This function fits a Cox proportional hazards model to analyze the association 
#      between MOFA factors and survival outcomes. The function constructs a survival 
#      object from the time and event columns, extracts factors from the MOFA model,
#      fits the Cox model, and summarizes the results in a table. It also generates 
#      a forest plot to visually represent the hazard ratios and their confidence intervals.
####################################################
surv.model <- function(mod, time, event, flag, outdir) {
  
    # Check if the output directory exists; if not, create it
    if (dir.exists(paste0(outdir, "/survival/", flag)) == F) {
        dir.create(paste0(outdir, "/survival/", flag))
    }

    # Create a survival object using time and event data
    SurvObject <- Surv(
        time = mod@samples_metadata[, time],
        event = mod@samples_metadata[, event]
    )

    # Obtain factors data from the MOFA model
    Z <- Reduce(rbind, mod@expectations$Z)
    Z <- Z[mod@samples_metadata$sample, ]

    # Fit the Cox proportional hazards model using the factors
    fit <- coxph(SurvObject ~ Z)

    # Summarize the fit and extract coefficients
    s <- summary(fit)
    coef <- s[["coefficients"]]

    # Create a results table with hazard ratios and confidence intervals
    sur.df <- data.frame(
        factor = factor(rownames(coef), levels = rev(rownames(coef))),
        p      = coef[, "Pr(>|z|)"],
        coef   = coef[, "exp(coef)"],
        lower  = s[["conf.int"]][, "lower .95"],
        higher = s[["conf.int"]][, "upper .95"]
    )

    # Write the results table to a text file
    write.table(sur.df,
        sep = "\t",
        paste0(outdir, "/survival/", flag, "/Hazard_ratio_factors.txt"),
        row.names = F,
        col.names = T
    )

    # Generate a forest plot for visualizing hazard ratios
    ggplot(sur.df, aes(x = factor, y = coef, ymin = lower, ymax = higher)) +
        geom_pointrange(col = "#619CFF", size = 1.3) +
        coord_flip() +
        scale_x_discrete() +
        labs(y = "Hazard Ratio", x = "") +
        geom_hline(aes(yintercept = 1), linetype = "dotted") +
        theme_bw() + 
        theme(text = element_text(size = 26))

    # Save the forest plot to a PDF file
    ggsave(
        filename = paste0(outdir, "/survival/", flag,
        "/Hazard_ratio_factors.pdf"),
        device = "pdf",
        dpi = 900, width = 10, height = 10
    )
}

####################################################
# Function: important.features
####################################################
# Input:
#      mod: MOFA model 
#      gmt: string specifying the path to a GMT file containing pathway information
#      nfeatures: integer indicating the number of top features to show in plots (default is 20)
#      th: numeric threshold for feature importance (default is 0.6)
#      met_anot: annotation matrix of the CpGs (methylation sites) information
# Output:
#      - Several text files containing important features extracted from the MOFA model for each omic type:
#          * Importance features for Expression and Methylation
#          * Top genes with positive and negative weights based on the specified threshold
#      - A summary table of the important features across different factors and pathways
# Description:
#      This function identifies and visualizes the most important features (genes, CpGs)
#      from the MOFA model by examining the weights associated with each factor.
#      It scales the weights, applies pathway annotation, and generates plots for the top
#      features. The results are saved to specified output directories.
####################################################
important.features <- function(mod, outdir, gmt,
                               nfeatures = 20, th = 0.6, met_anot) {

    genes_l <- list()  # Initialize a list to store important genes for each factor

    # Loop over each omic type in the MOFA model
    for (v in names(mod@dimensions$D)) {
        # Obtain weights and scale them within each omic and factor
        w <- get_weights(mod)[[v]]
        w <- as.data.frame(apply(w, 2, function(x) {
            x / max(abs(x))  # Scale weights to be between -1 and 1
        }))
        wi <- as.data.frame(w)  # Create a data frame from scaled weights

        # If the omic type is "Expression", annotate pathways
        if (v == "Expression") {
            wi <- Anot_path(wi, gmt)
        }

        # If the omic type is "Methylation", merge with annotation and process
        if (v == "Methylation") {
            wi$id <- gsub("_.*", "", rownames(wi))  # Extract the ID from row names
            wi <- merge(wi, met_anot, by.x = "id", by.y = "ID")  # Merge with methylation annotations
            
            wi <- wi %>% separate_rows(genes, sep = "_")  # Split genes into separate rows
            wi <- Anot_path(wi, gmt)  # Annotate pathways
        }

        # Write the importance features to a text file
        write.table(wi,
            paste0(outdir, "/features/", v, "_Importance_features.txt"),
            sep = "\t",
            col.names = T,
            row.names = F
        )

        # Loop through each factor and generate plots for the top features
        for (f in colnames(w)) {
            nfactor <- as.numeric(gsub("Factor", "", f))  # Extract factor number
            
            # Plot top weights for the current factor
            plot_top_weights(mod,
                view = v,
                factors = nfactor,
                nfeatures = nfeatures,  # Top number of features to highlight
                scale = T  # Scale weights from -1 to 1
            ) +
            theme(text = element_text(size = 18))

            # Create directory for the current factor if it doesn't exist
            if (dir.exists(paste0(outdir, "/features/", f)) == F) {
                dir.create(paste0(outdir, "/features/", f))
            }

            # Save the plot to a PDF file
            ggsave(
                filename = paste0(outdir, "/features/", f, "/", v,
                "_top_features.pdf"),
                device = "pdf", dpi = 600, width = 10, height = 10
            )

            # Select features by threshold for Expression omics
            if (v == "Expression") {
                genes_l[[f]][["Expression Positive"]] <- 
                    as.character(unique(wi$genes[wi[, f] >= th]))  # Features with positive weights
                genes_l[[f]][["Expression Negative"]] <- 
                    as.character(unique(wi$genes[wi[, f] <= -th]))  # Features with negative weights
                
                write.table(wi[order(abs(wi[, f]), decreasing = T),
                    c("genes", f, "Pathways")],
                    paste0(outdir, "/features/", f,"/", v,
                    "_top_genes.txt"),
                    sep = "\t",
                    row.names = F,
                    col.names = T
                )
            }

            # Select features by threshold for Methylation omics
            if (v == "Methylation") {
                genes_l[[f]][["Methylation Positive"]] <-
                    as.character(unique(wi$genes[wi[, f] >= th]))  # Features with positive weights
                genes_l[[f]][["Methylation Negative"]] <-
                    as.character(unique(wi$genes[wi[, f] <= -th]))  # Features with negative weights

                write.table(wi[order(abs(wi[, f]), decreasing = T),
                    c("id", f, "Pathways",
                    colnames(met_anot)[2:ncol(met_anot)])],
                    paste0(outdir, "/features/", f,"/", v,
                    "_top_genes.txt"),
                    sep = "\t",
                    row.names = F,
                    col.names = T
                )
            }
        }
    }

    # Create a summary table of important features across factors and pathways
    tabla_final <- MOFA_list2DF(genes_l, gmt)
    write.table(tabla_final,
        sep = "\t",
        paste0(outdir, "/features/Resume_table_features"),
        row.names = T, col.names = NA
    )
}
####################################################
# Function: gsea_factors
####################################################
# Input:
#      omic: type of omic data ("expression" or "methylation")
#      mod: MOFA model 
#      gmts: list of GMT files with pathway information
#      outdir: path to write output files
#      type: string indicating the type of analysis ("expression" or "methylation")
#      met_anot: annotation matrix of the CpGs (if type is "methylation", this is required)
#      flag: string identifier to append to output file names
# Output:
#      - Data frames containing GSEA results for each factor and pathway
#      - Text files documenting correlation results
#      - PDF dot plots visualizing GSEA results
# Description:
#      This function performs Gene Set Enrichment Analysis (GSEA) for factors in the 
#      MOFA model based on omic data. It calculates correlations between omic features 
#      and factors, applies multiple testing correction, and generates GSEA results 
#      for specified pathways.
####################################################
gsea_factors <- function(omic,
                         mod,
                         gmts,
                         outdir,
                         type,
                         met_anot = NULL,
                         flag) {
    
    # Validate the type input
    if (!type %in% c("expression", "methylation")) {
        stop("type should be expression or methylation")
    }

    # If the type is methylation, check if the annotation matrix is provided
    if (type == "methylation") {
        if (is.null(met_anot)) {
            stop("provide methylation annotation")
        }
    }

    # Extract the factors matrix from the MOFA model
    Z <- Reduce(rbind, mod@expectations$Z)
    Z <- Z[mod@samples_metadata$sample, ]
    Z <- as.data.frame(Z)
    Z <- Z[colnames(omic), ]  # Ensure Z has the same columns as the omic data

    # Loop through each factor in the model
    for (f in colnames(Z)) {
        
        # Calculate the correlation between the omic data and the current factor
        df_corr <- Corr_factor(omic, Z, f)
        
        # If the analysis type is methylation, merge with the annotation data
        if (type == "methylation") {
            df_corr <- merge(df_corr, met_anot, by.x = "row.names", by.y = "ID")
            df_corr <- df_corr %>% separate_rows(genes, sep = "_")  # Split multiple genes into separate rows
            df_corr <- df_corr[!duplicated(df_corr$genes), ]  # Remove duplicate genes
            df_corr <- as.data.frame(df_corr)
            rownames(df_corr) <- df_corr$genes  # Set row names to genes
        }
        
        # Adjust p-values using the FDR method
        df_corr$p.adj <- p.adjust(df_corr$p.value, method = "fdr")
        df_corr <- df_corr[order(df_corr[,'p.value'], decreasing = FALSE), ]  # Order by p-value
        
        # Write the correlation results to a text file
        write.table(df_corr,
                    sep = '\t',
                    paste0(outdir, "/features/", f, "/",
                           "Correlation_", flag, ".txt"),
                    row.names = TRUE,
                    col.names = NA)

        # Prepare the ranking for GSEA
        rnk <- df_corr[order(df_corr[,'cor'], decreasing = TRUE), 'cor']
        names(rnk) <- rownames(df_corr[order(df_corr[,'cor'], decreasing = TRUE), ])
        rm(df_corr)  # Free up memory by removing df_corr
        
        # Loop through each GMT file and perform GSEA
        for (gmt in gmts) {
            name <- gsub(".*/|.gmt", "", gmt)  # Extract the name of the GMT file

            # Perform GSEA
            gseahk <- gsea(rnk = rnk,
                           gmt = gmt,
                           title = paste0(outdir, "/features/", f, "/", name,
                                          "_dotplot_", flag, ".pdf"))

            gseahk$Factor <- f  # Add the factor to the GSEA results

            # Write the GSEA results to a text file
            write.table(gseahk,
                        paste0(outdir, "/features/", f, "/", name,
                               "_GSEA_", flag, ".txt"),
                        sep = '\t',
                        row.names = FALSE,
                        col.names = TRUE)
        }
    }
}


####################################################
# Function: obtain_Z
####################################################
# Input:
#      mod: MOFA model (contains factor expectations)
# Output:
#      - A matrix with samples as rows and factors as columns
# Description:
#      This function extracts the expected values of factors from a MOFA model
#      and organizes them into a matrix format, where each row corresponds 
#      to a sample and each column corresponds to a factor.
####################################################
obtain_Z <- function(mod) {
    # Combine the factor expectations into a single matrix
    Z <- Reduce(rbind, mod@expectations$Z)
    
    # Reorder the matrix to match the sample order in the metadata
    Z <- Z[mod@samples_metadata$sample, ]
    
    # Return the resulting matrix
    return(Z)
}


#function:Anot_path
#input:
#      df: matrix with genes as rownamaes or a column with genes to anotate in pathways
#      gmt: gmt file 
#     

#output: # df matrix with pathway anotations

#description:Function to anotate the genes of df matrix
####################################################
# Function: Anot_path
####################################################
# Input:
#      df: A data frame or matrix with genes as row names or a column containing genes to annotate with pathways
#      gmt: A GMT file (gene set file) where each gene is mapped to its corresponding pathways
# Output:
#      - A data frame that includes the original data along with pathway annotations for each gene
# Description:
#      This function annotates the genes in the input data frame or matrix with their corresponding
#      pathways as defined in the provided GMT file. It merges the pathway information with the
#      input data based on gene matches.
####################################################
Anot_path <- function(df, gmt) {
    
    # Extract unique genes from the GMT file
    genes <- unlist(gmt)
    genes <- unique(genes)

    # Initialize a list to hold pathway annotations
    genes_l <- list()

    # Loop through each gene and associate it with its pathways
    for (x in genes) {
        for (y in names(gmt)) {
            if (x %in% gmt[[y]]) {
                genes_l[[x]] <- append(genes_l[[x]], y)  # Append pathway names to the gene
            }
        }
    }

    # Convert the list of pathways to a string format
    genes_l <- lapply(genes_l, toString)

    # Combine the gene-pathway associations into a data frame
    genes <- do.call(rbind, genes_l)
    colnames(genes) <- "Pathways"
    genes <- as.data.frame(genes)

    # Ensure that the 'genes' column exists in the gene-pathway data frame
    if ("genes" %in% colnames(genes) == FALSE) {
        genes$genes <- rownames(genes)
    }

    # Ensure that the 'genes' column exists in the input data frame
    if (!"genes" %in% colnames(df)) {
        df$genes <- rownames(df)
    }
    
    # Merge the original data frame with the pathway annotations
    return(merge(df, genes, by.x = "genes", by.y = "genes", all.x = TRUE))
}
####################################################
# Function: MOFA_list2DF
####################################################
# Input:
#      genes_l: A list of genes to be annotated with pathways
#      gmt: A GMT file (gene set file) mapping genes to their respective pathways
# Output:
#      - A data frame containing gene annotations, including pathways and factors associated with each gene
# Description:
#      This function takes a list of genes and maps them to their corresponding pathways 
#      using the provided GMT file. It also organizes the data into a structured format 
#      that includes information about the factors and omic types associated with each gene.
####################################################
MOFA_list2DF <- function(genes_l, gmt) {

    # Convert the list of genes into a data frame
    df_all <- as.data.frame(unlist(genes_l))
    
    # Split the row names to extract factors and omic types
    features <- lapply(rownames(df_all), function(x) {
        strsplit(x, "\\.")
    })
    
    # Extract factor names and omic types
    factor <- unlist(lapply(features, function(x) {
        return(x[[1]][1])
    }))
    omic <- unlist(lapply(features, function(x) {
        return(x[[1]][2])
    }))
    omic <- gsub("[0-9]*", "", omic)  # Remove numerical characters from omic types

    # Add factor and omic type information to the data frame
    df_all$Factor <- factor
    df_all$Omic <- omic
    colnames(df_all) <- c("genes", colnames(df_all)[2:ncol(df_all)])
    rownames(df_all) <- NULL
    
    # Combine omic types into a single string for each gene and factor
    df_all <- transform(df_all, lab = ave(Omic, genes, Factor, FUN = toString))

    # Keep only relevant columns and remove duplicates
    df_all <- df_all[, c(1, 2, 4)]
    df_all <- df_all[!duplicated(df_all), ]

    # Reshape the data frame to have genes as rows and factors as columns
    df_all <- dcast(data = df_all, formula = genes ~ Factor, value.var = "lab")

    # Count the total number of non-missing values for each gene
    df_all$Total <- apply(df_all, 1, function(x) {
        return(length(x[complete.cases(x)]))
    })

    # Annotate the data frame with pathways using the Anot_path function
    df_all <- Anot_path(df_all, gmt)

    # Order the data frame by the Total column in descending order
    df_all <- df_all[order(df_all$Total, decreasing = TRUE), ]
    df_all$Total <- df_all$Total - 1  # Adjust Total count
    rownames(df_all) <- df_all$genes

    # Remove the genes column from the output
    df_all <- df_all[, -1]

    return(df_all)  # Return the annotated data frame
}

####################################################
# Function: Corr_factor
####################################################
# Input:
#      omic: A matrix where rows represent features (e.g., genes) and columns represent samples
#      Z: A matrix with rows representing samples and columns representing factors
#      f: A specific factor to correlate with features
# Output:
#      - A data frame containing the correlation coefficients and p-values 
#        between each feature in the omic matrix and the specified factor
# Description:
#      This function computes the correlation between features in the omic matrix 
#      and a specified factor in the Z matrix, returning the correlation coefficients 
#      and associated p-values for each feature.
####################################################
Corr_factor <- function(omic, Z, f) {
    # Apply correlation test for each feature against the specified factor
    crexp <- lapply(as.data.frame(t(omic)), function(x) {
        crt <- cor.test(x, Z[, f])  # Perform correlation test
        return(list("cor" = crt$estimate[[1]], "p.value" = crt$p.value[[1]]))
    })

    # Convert the list of results into a data frame
    dfcorr <- as.data.frame(do.call(rbind, crexp))
    dfcorr$cor <- unlist(dfcorr$cor)         # Extract correlation coefficients
    dfcorr$p.value <- unlist(dfcorr$p.value) # Extract p-values

    return(dfcorr)  # Return the resulting data frame
}


####################################################
# Function: gsea
####################################################
# Input:
#      rnk: A named vector containing the ranking of features (correlation values, for example)
#      gmt: A GMT file containing pathway information (gene set)
#      title: The name of the output file to save the GSEA results plot
# Output:
#      - A data frame containing GSEA results for the given gene rankings
# Description:
#      This function performs Gene Set Enrichment Analysis (GSEA) on a ranked list 
#      of features based on their correlation with specified factors, using pathways defined in a GMT file.
####################################################
gsea <- function(rnk, gmt, title) {
    
    # Load the pathways from the GMT file
    gmt <- gmtPathways(gmt)
    
    # Perform GSEA using the fgsea package
    gsea <- fgsea(
        pathways = gmt,   # The gene sets from the GMT file
        stats = rnk,      # The ranked gene statistics
        minSize = 15,     # Minimum size of the pathways to consider
        maxSize = 500     # Maximum size of the pathways to consider
    )

    # Create a leading edge column by collapsing the leading edge genes into a single string
    gsea$leadingEdge <- unlist(lapply(gsea$leadingEdge, function(x) {
        paste(x, collapse = "/")
    }))
    
    # Check if there are significant pathways (adjusted p-value < 0.05)
    if (nrow(gsea[gsea$padj < 0.05, ]) > 0) {
        n <- nrow(gsea[gsea$padj < 0.05, ])
        
        # Limit the number of pathways displayed to a maximum of 10
        if (n > 10) {
            n <- 10
        }

        # Filter significant pathways and order them by p-value
        gsea <- gsea[which(gsea$padj < 0.05), ]
        gsea <- gsea[order(abs(gsea$padj)), ]
        gsea$pathway <- factor(gsea$pathway, levels = gsea$pathway)

        # Generate a scatter plot for the top pathways
        ggplot(data = gsea[1:n, ], aes(x = NES, y = pathway, size = size,
                                        color = -1 * log10(padj))) +
            geom_point() +
            scale_size_area(max_size = 10) +
            scale_colour_gradient(low = "green", high = "red") +
            theme(text = element_text(size = 26))

        # Save the plot to a PDF file
        ggsave(
            filename = title,  # Name of the output file
            device = "pdf", 
            dpi = 500, 
            width = 22, 
            height = 12
        )
    }

    return(gsea)  # Return the GSEA results
}
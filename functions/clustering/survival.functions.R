####################################################
# Function: Hazard_Cox
#
# Description:
#      This function creates Cox proportional hazards models to evaluate the impact
#      of clinical features on overall survival (OS) in a given dataset.
#
# Input:
#      base: data.frame with samples as rows and clinical information in columns
#      class: column name of feature of interest (e.g., neuroendocrine subtypes)
#      OS: overall survival time from diagnosis to last interview
#      OS_event: Overall survival status (0 = alive, 1 = dead)
#      covariates: vector of column names for additional covariates to adjust the model (default: NULL)
#      outdir: Path to write results (output directory)
#      flag: Name to add in the name of output files
#
# Output:
#      A table containing hazard ratios, confidence intervals, p-values for each level in each variable,
#      and a forest plot saved as a PDF file.
####################################################


Hazard_Cox <- function(base,
class,
OS,
OS_event,
covariates = NULL,
outdir,
flag
) {

  # Rename OS and OS_event columns in the dataframe for clarity
  colnames(base) <- gsub(OS, 'OS', colnames(base))
  colnames(base) <- gsub(OS_event, 'OS_event', colnames(base))

  # Ensure OS and OS_event are numeric
  
     base$OS <- as.numeric(base$OS)
    base$OS_event <- as.numeric(base$OS_event)
    
    
    # Generate the survival object
    base$Surv_OS <- with(base, Surv(base$OS, base$OS_event == 1))
    
  # Fit Cox proportional hazards model
  
    if (is.null(covariates)){
    fit <- coxph( as.formula(paste('Surv_OS ~ ', class)), data = base)

    }else{
      fit <- coxph( as.formula(paste('Surv_OS ~ ', class,' + ',paste(covariates, collapse = ' + '))) , data = base)
    }

      # Summary of the fitted model
    sum <- summary(fit)
    coef <- sum[["coefficients"]]

      # Get model coefficients and variance-covariance matrix
    coef_model <- coef(fit)
    cov <- vcov(fit)

      # Global p-value calculation using Wald test
   global_p <- lapply(c(class,covariates), function(x) {
        
        indices <- grep(paste0("^",x),names(coef_model))

        wald <- wald.test(b = coef_model[indices],
        Sigma = as.matrix(cov[indices,indices]),
        L = diag(length(indices)))
        
        return(rep(wald$result$chi2[["P"]],length(indices)))
    })

  # Generate results table
   
      sur.df <- data.frame(
        cluster = factor(rownames(coef), levels = rev(rownames(coef))),
        p      = coef[, "Pr(>|z|)"],
        coef   = coef[, "exp(coef)"],
        lower  = sum[["conf.int"]][, "lower .95"],
        higher = sum[["conf.int"]][, "upper .95"],
        global_p = unlist(global_p)
    )

  # Save results table to a text file
    write.table(sur.df,
        sep = "\t",
        paste0(outdir, "/", flag, "_Hazard_ratio_factors.txt"),
        row.names = F,
        col.names = T
    )
    
      # Create a forest plot for hazard ratios
    ggplot(sur.df, aes(x = cluster, y = coef, ymin = lower, ymax = higher)) +
        geom_pointrange(col = "#619CFF", size = 1.3) +
        coord_flip() +
        scale_x_discrete() +
    labs(y = "Hazard Ratio", x = "") +
        geom_hline(aes(yintercept = 1), linetype = "dotted") +
        theme_bw() + 
        theme(text = element_text(size = 26))
  # Save forest plot as a PDF
    ggsave(
        filename = paste0(outdir, "/", flag, "Hazard_ratio_factors.pdf"),
        device = "pdf",
        dpi = 900, width = 10, height = 10
    )

    # Return results data frame
    return <- sur.df
    }

####################################################
# Function: Kapplan_Meyer
#
# Description:
#      This function creates a Kaplan-Meier survival plot, calculates the survival probabilities
#      at specified time points, and performs a log-rank test to compare survival distributions.
#
# Input:
#      base: data.frame with samples as rows and clinical information in columns.
#      class: column name of the feature of interest (used for stratifying the survival analysis).
#      OS: overall survival time from diagnosis to last interview.
#      OS_event: Overall survival status, where 0 indicates alive and 1 indicates dead.
#      outdir: Path to write results.
#      flag: Name to add to the names of output files.
#      custcolor: vector of colors to use for levels in the class variable.
#      width: width of the output PDF plot.
#      height: height of the output PDF plot.
#
# Output:
#      A Kaplan-Meier plot with p-value from the log-rank test and survival probabilities
#      at specified time points (60 and 120 months).
####################################################

Kapplan_Meyer <- function(base,
class,
OS,
OS_event,
outdir,
flag,
custcolor,
width= 8,
height=8
){
  
    # Replace column names for consistency
  colnames(base) <- gsub(OS, 'OS', colnames(base))
  colnames(base) <- gsub(OS_event, 'OS_event', colnames(base))
  
  # Remove rows with missing values in OS or OS_event
  base <- base[complete.cases(base$OS),]
  base <- base[complete.cases(base$OS_event),]

    # Ensure survival time and event are numeric  
    base$OS <- as.numeric(base$OS)
    base$OS_event <- as.numeric(base$OS_event)
    
    # Create survival object
    base$Surv_OS <- with(base, Surv(base$OS, base$OS_event))

    # Generate survival formula for Kaplan-Meier analysis
    form <- as.formula (paste( "Surv_OS ~ ", class))
  
    # Fit Kaplan-Meier survival curves
    km.by.SG <- surv_fit (form,data = base , conf.type = "log-log")

   # Perform log-rank test to compare survival curves
       
        diff_SG <- survdiff(as.formula (paste( 'base$Surv_OS ~',
        class)),
        data = base, rho = 0)
  # Calculate overall survival for the entire population
        svc <- survfit(as.formula ( "Surv_OS ~ 1"),data = base , conf.type = "log-log")
        dfc <- summary(svc, c(60,120))
        dfc <-  data.frame("time" = dfc[[2]],
        "supervivencia" = dfc[[6]],
        "clase" = "Population")

  # Calculate survival probabilities at 60 and 120 months for the stratified data
        sv <- survfit(form, data= base)
        df <- summary(sv, c(60,120))
        df <- data.frame("time" = df[[2]],
        "supervivencia" = df[[6]],
        "clase" = df[[10]])
  # Combine overall population and stratified survival data
        df <- rbind(df,dfc)

  # Save survival probabilities to a text file
        write.table(df,
        sep = "\t",
        paste0(outdir, "/", flag, "_Supervivencia_tabla.txt"),
        row.names = F,
        col.names = T
        )
  # Calculate p-value from log-rank test
        pval <- signif(pchisq(diff_SG$chisq,
        length(diff_SG$n)-1,
        lower.tail = FALSE),
        digits = 3)
  # Create PDF for Kaplan-Meier plot        
        pdf(paste0(outdir, "/", flag, "_",
        "Kapplan_Meyer.pdf"), width = width, height = height)
  # Generate and customize Kaplan-Meier plot
        p <- ggsurvplot(km.by.SG,
        data = base,
        pval = pval,
        risk.table = T,
        xlab = "Time",
        ylab = "Overall Survival",
        palette = custcolor,
        font.main = c(18, "bold"),
        font.x = c(18, "bold"),
        font.y = c(18, "bold"),
        font.caption = c(18, "bold"), 
        font.legend = c(18, "bold"), 
        font.tickslab = c(18, "bold"),
        pval.size = 10,
        fontsize = 7) 


# Print the plot
        print(p)

# Close the PDF device     
        dev.off()
}

####################################################
# Function: cut_survival
#
# Description:
#      This function modifies the overall survival time and event status 
#      for samples in the provided dataset, cutting the survival time 
#      to a specified maximum (in months). If the overall survival time 
#      exceeds this maximum, the event status is set to 0 (indicating alive).
#
# Input:
#      base: data.frame with samples as rows and clinical information in columns.
#      cut: numeric value representing the maximum months to retain in survival time.
#      time: character string indicating the name of the column containing overall survival time.
#      event: character string indicating the name of the column containing overall survival status.
#
# Output:
#      A modified version of the input data frame (base) where:
#      - The survival time is capped at the specified cut value.
#      - The event status is adjusted accordingly based on the cut value.
####################################################

cut_survival <- function(base, cut,
                         time = "OS.time",
                         event = "EXITUS") {
  
  # Adjust the event status for samples where the survival time exceeds the cut-off
  base[, event][base[, time] > cut] <- 0  # Set event to 0 (alive) if time exceeds cut
  
  # Cap the survival time at the specified cut value
  base[, time][base[, time] > cut] <- cut
  
  # Return the modified data frame
  return(base)
}


####################################################
# Function: generate_medians 
#
# Description:
#      This function determines whether each sample's feature value is above or below 
#      the median value of specified features. It outputs a data frame with the samples 
#      as rows and the features as columns, indicating whether each sample is "Alto" (High) 
#      or "Bajo" (Low) relative to the median.
#
# Input:
#      data: matrix with features as rows and samples as columns.
#      features: character vector containing the names of the features for which to calculate medians.
#
# Output:
#      A data.frame where:
#      - Rows represent samples.
#      - Columns represent features.
#      - Each cell contains "Alto" if the sample's value is above the median, 
#        or "Bajo" if it is below the median for that feature.
####################################################

generate_medians <- function(data, features) {
  
    # Calculate medians and determine high/low status for specified features
    medians <- lapply(features, function(x) {
        
        # Compute the median for the current feature
        med <- median(as.numeric(data[x, ]))
        
        # Compare each sample's value to the median and label accordingly
        result <- ifelse(as.numeric(data[x, ]) > med, "Alto", "Bajo")
        
        return(result)  # Return the high/low status for the feature
    })

    # Name the list elements with feature names
    names(medians) <- features
    
    # Combine the results into a data frame
    medians <- do.call(cbind, medians)
    
    # Set the row names of the resulting data frame to the sample names (column names of the original data)
    rownames(medians) <- colnames(data)
    
    return(medians)  # Return the final data frame
}


####################################################
# Function: cox.multiple 
#
# Description:
#      This function conducts multiple Cox proportional hazards regression analyses 
#      for a set of genes using a provided survival dataset. It adjusts for covariates 
#      and outputs the hazard ratios and p-values for each gene.
#
# Input:
#      base: data.frame containing survival data with samples as rows and clinical information in columns.
#      data: matrix of gene expression values with genes as rows and samples as columns.
#      OS: character string indicating the column name for overall survival time.
#      OS_event: character string indicating the column name for overall survival event status (0 = alive, 1 = dead).
#      covariates: vector of column names to adjust the model (optional).
#      outdir: directory path for writing results.
#      flag: name prefix for output files.
#      th: threshold for p-values to determine significance (default = 0.05).
#      adj: logical indicating whether to adjust p-values for multiple testing (default = TRUE).
#
# Output:
#      A data.frame with hazard ratios, confidence intervals, and p-values for each gene.
####################################################



cox.multiple <- function(base,
data,
OS,
OS_event,
covariates = NULL,
outdir,
flag,
th = 0.05,
adj = T) {


  # Filter out complete cases for OS and OS_event
base <- base[complete.cases(base[,OS]),]
base <- base[complete.cases(base[,OS_event]),]

  
# Z-score the data
data <- zscore.rows2(data)

  # Ensure base and data are in the same order
data <- data[,rownames(base)]
data <- as.data.frame(t(data))

  # Create survival object
base$Surv_OS <- with(base, Surv(base[,OS], base[,OS_event] == 1))

  # Initialize a list to store results for each gene
sur.df <- lapply(colnames(data), function(gene){

# Merge survival data with gene expression data for the current gene
df <- merge(base,data[,gene,drop = F], by = "row.names")
#We generate the rurvival objectbase

    # Fit the Cox model
    if (is.null(covariates)){
    fit <- coxph( as.formula(paste('Surv_OS ~ ', gene)), data = df)

    }else{
      fit <- coxph( as.formula(paste('Surv_OS ~ ', gene,' + ', 
      paste(covariates, collapse = ' + '))) , data = df)
    }
    
    
# Get the summary of the Cox model
   sum <- summary(fit)
   coef <- sum[["coefficients"]]

    # Create a data frame to store results for the current gene
   result <- data.frame(
      gene = gene,
      p = coef[gene, "Pr(>|z|)"],                # p-value for the gene
      coef = coef[gene, "exp(coef)"],            # Hazard ratio
      lower = sum[["conf.int"]][gene, "lower .95"],  # Lower CI
      higher = sum[["conf.int"]][gene, "upper .95"]  # Upper CI
    )
    
    return(result)  # Return results for the current gene

})

  # Combine results from all genes into a single data frame
sur.df <- as.data.frame(do.call(rbind,sur.df))

  # Adjust p-values for multiple testing
sur.df$p.adj <- p.adjust(sur.df$p, method = "fdr")

  # Filter results based on adjusted p-values if adj is TRUE
if (adj == T){
sur.df <- sur.df[sur.df$p.adj < th,]
}else{
sur.df <- sur.df[sur.df$p < th,]
}

  # Order results by adjusted p-values
sur.df <- sur.df[order(sur.df$p.adj),]


  # Write the results to a text file
write.table(
    sur.df,
    sep = "\t",
    paste0(outdir,"/",flag, ".txt"),
    row.names = F,
    col.names = T)

return(sur.df) # Return the final results data frame
}

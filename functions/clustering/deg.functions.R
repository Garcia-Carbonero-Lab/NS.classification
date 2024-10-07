####################################################
# Function: deg.function
#
# Input:
#      exp: A matrix where rows represent features (e.g., genes) and columns represent samples.
#      base: A data frame where rows are samples and columns are various features (e.g., clinical data).
#      cov: Covariates to adjust for in the differential expression analysis (optional).
#      group: Column name in 'base' used to define the groups for differential expression analysis.
#      control: Label of the control group in 'group' (optional; default is the first factor level).
#      outdir: Directory path to save the results and plots.
#      flag: String identifier to append to file names (useful for distinguishing different analyses).
#
# Output:
#      A result table of the differential expression analysis performed using the limma package.
#      The function writes the top results to text files in the specified output directory.
#
# Description:
#      This function performs differential expression analysis using the limma package. It generates 
#      a design matrix based on the input data and performs linear modeling, followed by contrast
#      fitting and hypothesis testing to identify differentially expressed features between groups.
####################################################


deg.function <- function(exp,
base, group,
cov = NULL,
outdir,
control = NULL,
flag) {


# Convert the group column in 'base' to a factor, ensuring it's treated as categorical
base[,group] <- as.factor(as.character(base[,group]))

# Subset 'base' to only include samples present in the expression matrix
base <- base[colnames(exp),]

# Set control group as the first level if not specified
if (is.null(control)){

    control <- levels(base[,group])[1]
}

############### Design Matrix Creation ###############
# If no covariates are provided, create a simple design matrix with only the group variable
if (is.null(cov)){
  design <- data.frame(Group = base[,group])
  # Relevel the group factor so that the control group is the reference level
  design$Group <- relevel(design$Group, control)
  
  # Create the design matrix with no intercept (to compare group levels directly)
  design <- model.matrix( ~ 0 + Group, data = as.data.frame(design))
  }else{
    # If covariates are provided, add them to the design matrix
    design <- data.frame(Group = base[,group])
    # Add each covariate to the design dataframe
    for (f in cov){design[,f] <- base[,f]}
    
    # Relevel the group factor so that the control group is the reference level
    design$Group <- relevel(design$Group, control)
    # Create the design matrix with covariates included
    design <- model.matrix(as.formula(paste0("~ 0 + Group + ",
    paste0(cov, collapse = " + "))), data = as.data.frame(design))
  }
  
  ############### Contrast Matrix Creation ###############
    # Generate all possible pairwise comparisons between group levels

  CN <- combn(colnames(design)[1:length(levels(base[,group]))],2,
  simplify = FALSE)
  
  # Combine the group levels into contrast strings (e.g., GroupA - GroupB)
  CN <- unlist(lapply(CN,function(x){paste0(x, collapse = "-")}))
  
   # Set the row names of the design matrix to match the sample names
  rownames(design) <- rownames(base)

  # Fit the linear model to the expression data using the design matrix
  fit<-lmFit(exp, design)
  
  
############### Contrast Fitting and Hypothesis Testing ###############
    # Create a contrast matrix based on the pairwise comparisons (CN)

  contrast.matrix <- makeContrasts(contrasts= CN,
  levels= design)
  
# Apply the contrasts to the fitted model  
  fit2 <- contrasts.fit(fit, contrast.matrix)
 
# Perform empirical Bayes moderation to compute statistics for each contrast 
  fit2 <- eBayes(fit2)

 ############### Result Extraction and Output ###############
res <- list() # Initialize a list to store the results for each comparison

# Loop through each contrast (pairwise comparison)
for (comb in CN){

# Extract the top table of differentially expressed features for the current comparison
result <- topTable(fit2, coef=comb, adjust="BH",number = nrow(exp))


# Write the result to a file in the specified output directory
write.table(result,
paste0(outdir,"/TopTable_",flag, "_",comb ,".txt"),
sep = "\t",
row.names = T,
col.names = NA)

# Store the result in the list
res[[comb]] <- result
}

# Combine the individual results into a single data frame
res <- do.call(rbind, res)

# Return the complete result set
return(res)
}

####################################################
# Function: deg.plot
#
# Input:
#      data: A matrix where rows represent features (e.g., genes) and columns represent samples.
#      base: A data frame where rows represent samples and columns represent features (e.g., clinical data).
#      group.x: The column name in 'base' that will be used for the x-axis in the boxplot.
#      group.fill: The column name in 'base' that will be used for the fill color in the boxplot.
#      color.fill: A vector of colors to represent the levels of the 'group.fill' variable.
#      cov: Optional covariates to adjust for in the limma model (NULL by default).
#      outdir: The path where the output plots and results will be saved.
#      flag: A string identifier used to create unique filenames for the output.
#      genes: A vector of gene names to include in the boxplots.
#      filter.table: A boolean value to determine if the result table should be filtered to include only the provided genes.
#
# Output:
#      - Boxplots showing the distribution of the selected genes across different subtypes, annotated with p-values.
#      - A table summarizing the limma results, including p-values, significance levels, and gene names.
#
# Description:
#      This function generates boxplots for selected genes, using the groupings defined by 'group.x' and
#      'group.fill', and performs differential expression analysis using the limma package. The results
#      are displayed in the boxplots with adjusted p-values, and the significance levels are marked.
#      The function saves the plots as PDF files and writes the differential expression results to a text file.
####################################################

deg.plot <- function(data,
base,
group.x,
group.fill,
color.fill,
cov = NULL,
flag,
outdir,
genes,
filter.table = T) {

############### Data Preparation ###############
    # Replace problematic characters in gene and sample names


genes <- gsub("-","_",genes)
colnames(data) <- gsub("-", "_", colnames(data))
genes <- gsub("\\.","_",genes)
colnames(data) <- gsub("\\.", "_", colnames(data))


# Select only the genes present in the dataset
genes <- genes[genes %in% colnames(data)]



 # Create an output directory if it doesn't exis
if (dir.exists(paste0(outdir, "/", flag)) == F) {

    dir.create(paste0(outdir, "/", flag))

}


# Ensure the rows of 'data' match the rownames of 'base'
data <- data[rownames(base),,drop = F]

# Merge the data with the relevant columns from 'base'
df<- cbind(data,base[,c(group.x,group.fill)])
df[,group.x] <- as.factor(as.character(df[,group.x]))
df[,group.fill] <- as.factor(as.character(df[,group.fill]))

############### Differential Expression and Plot Generation ###############
    # Check if 'group.x' and 'group.fill' are different, and handle accordingly

if (group.x != group.fill){

# Convert group variables to factors in 'base'
base[,group.x] <- as.factor(as.character(base[,group.x]))
base[,group.fill] <- as.factor(as.character(base[,group.fill]))

# Perform differential expression analysis for each level of 'group.x'
degs <- lapply(levels(base[,group.x]), function(x){

dfp <- base[base[,group.x] == x,]

# Call deg.function to perform limma analysis
degp <- deg.function(as.data.frame(t(data))[,rownames(dfp)],
dfp,
group = group.fill,
cov = cov,
outdir = paste0(outdir, "/", flag),
flag = x)

## Prepare the limma results for use in ggpubr plotting
degp$Row.names <- rownames(degp)
degp <- as_tibble(degp)
fdegp <- degp %>%
separate(Row.names, c("Group", "gene"), sep = "\\.") %>%
separate(Group, c("group1", "group2"), sep = "-") %>%
mutate(
  group1 = str_replace(group1, "Group", ""),
  group2 = str_replace(group2, "Group", "")

) %>%
rename(
  p = P.Value,
  p.adj = adj.P.Val 
) 

# Filter the table to only include the selected genes, if applicable
if (filter.table == T) {
fdegp <- fdegp%>%
filter(
  gene %in% genes
)
}

# Annotate the adjusted p-values with significance stars
fdegp[,"p.adj.signif"] <- rep("ns", nrow(fdegp))

if (nrow(fdegp[fdegp$p.adj < 0.001,]) > 0){
fdegp[fdegp[, "p.adj"] < 0.001,"p.adj.signif"] <- "***"
}

if (nrow(fdegp[fdegp$p.adj < 0.01,]) > 0){
fdegp[fdegp[, "p.adj"] < 0.01 & fdegp[, "p.adj"] > 0.001,"p.adj.signif"] <- "**"
}

if (nrow(fdegp[fdegp$p.adj < 0.05,]) > 0){
fdegp[fdegp[, "p.adj"] < 0.05 & fdegp[, "p.adj"] > 0.01,"p.adj.signif"] <- "*"
}

# Assign the current group level to the 'group.x' column
fdegp[,group.x] <- as.factor(rep(x,nrow(fdegp)))
return(fdegp)

})

# Combine the differential expression results from each group level
degs <- degs %>%
bind_rows(degs) %>%
distinct

# Write the combined limma results to a file
write.table(degs,
paste0(outdir, "/",flag, "/Limma_genes_primary_tumor_Subtype.txt"),
sep = "\t",
row.names = F,
col.names = T)


}else{


# Handle case when 'group.x' equals 'group.fill'
base[,group.fill] <- as.factor(as.character(base[,group.fill]))

# Perform limma analysis on the entire dataset

degp <- deg.function(as.data.frame(t(data)),
base,
group = group.fill,
cov = cov,
outdir = paste0(outdir, "/", flag),
flag = "Complete")
degp$Row.names <- rownames(degp)

degp <- as_tibble(degp)

# Prepare limma results for plotting
fdegp <- degp %>%
separate(Row.names, c("Group", "gene"), sep = "\\.") %>%
separate(Group, c("group1", "group2"), sep = "-") %>%
mutate(
  group1 = str_replace(group1, "Group", ""),
  group2 = str_replace(group2, "Group", "")

) %>%
rename(
  p = P.Value,
  p.adj = adj.P.Val 
)

if (filter.table == T) {
fdegp <- fdegp%>%
filter(
  gene %in% genes
)
}

# Annotate p-values with significance stars

fdegp[,"p.adj.signif"] <- rep("ns", nrow(fdegp))

if (nrow(fdegp[fdegp$p.adj < 0.001,]) > 0){
fdegp[fdegp[, "p.adj"] < 0.001,"p.adj.signif"] <- "***"
}

if (nrow(fdegp[fdegp$p.adj < 0.01,]) > 0){
fdegp[fdegp[, "p.adj"] < 0.01 & fdegp[, "p.adj"] > 0.001,"p.adj.signif"] <- "**"
}

if (nrow(fdegp[fdegp$p.adj < 0.05,]) > 0){
fdegp[fdegp[, "p.adj"] < 0.05 & fdegp[, "p.adj"] > 0.01,"p.adj.signif"] <- "*"
}

degs <- fdegp

# Write the results to a file
write.table(degs,
paste0(outdir, "/", flag, "/Limma_genes_primary_tumor_Subtype.txt"),
sep = "\t",
row.names = F,
col.names = T)

}

############### Boxplot Creation ###############
    # Generate boxplots for each gene

for (row in genes) {

#Prepare data for ggpurb
if (group.x != group.fill){
# Prepare data for plotting
dfa <- df[,c(row,group.x,group.fill)]

# Perform Wilcoxon test using ggpubr
stat.test <- dfa %>%
  group_by_at(vars(group.x)) %>%
  wilcox_test(as.formula(paste0(row, "~", group.fill))) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj")

stat.deg <- degs %>%
filter( gene == row)

# Add p-values from limma analysis to stat.test for use in plotting
stat.test <- stat.test %>%
  add_xy_position(x = group.x, dodge = 0.8)



stat.test <- stat.test %>% 
unite(col= "order", c( as.name(group.x), group1:group2), sep = "_",
remove = FALSE) 

stat.deg <- stat.deg %>% 
unite(col= "order", c(
group1:group2),
sep = "_",
remove = FALSE)%>%
distinct()


}else{

dfa <- df[,c(row,group.fill)]
 
stat.test <- dfa %>%
  wilcox_test(as.formula(paste0(row, "~", group.fill))) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj")


stat.deg <- degs %>%
filter( gene == row)


stat.test <- stat.test %>%
  add_xy_position(x = group.x, dodge = 0.8)



stat.test <- stat.test %>% 
unite(col= "order", c(group1:group2), sep = "_",
remove = FALSE) 



stat.deg <- stat.deg %>% 
unite(col= "order", c(
group1:group2),
sep = "_",
remove = FALSE)%>%
distinct()



}

# Add p-values from limma analysis to stat.test for use in plottingstat.test <- stat.test %>%
mutate(
  p.adj.signif = stat.deg$p.adj.signif
)

stat.test[[group.fill]] <- NA

# Generate the boxplot
bxp <- ggplot(
  df, aes(x = .data[[group.x]],
  y = .data[[row]], fill= .data[[group.fill]])) +
  geom_boxplot(width = 0.6, outlier.shape = NA,
  aes(fill= .data[[group.fill]])) +

  scale_fill_manual(values = as.list(c(color.fill,"blue", "green"))) +

  geom_point(position = position_jitterdodge(jitter.width = 0.2,
                                             dodge.width = 0.7)) +
  
  stat_pvalue_manual(
  stat.test,  label = "p.adj.signif", tip.length = 0, hide.ns = T
) +


  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +

theme_classic() +

theme(
        text = element_text(size = 26),
        legend.key.size = unit(2, "cm"),
        legend.text = element_text(size = 26)
        ) +

        labs(x = "", y = row)


bxp$layers[[3]]$aes_params$size <- 1.2
bxp$layers[[3]]$aes_params$label.size <- 12 
bxp

# Save the plot as a PDF
ggsave(
        plot = bxp, filename = paste0(
            outdir, "/", flag, "/Boxplot.subtype_", row, ".pdf"
        ),
        device = "pdf",
        dpi = 300, width = 16,
        height = 12
    )

}


}
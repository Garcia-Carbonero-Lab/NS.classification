


####################################################
# Function: anotate.gmt
# Input:
#      genesets: List. A list of gene sets where each element contains a list of genes.
#      expression: Matrix. Gene expression matrix with gene symbols as row names and samples as columns.
#
# Output:
#      Annotated gene sets with gene symbols replaced by those matching the expression matrix.
#
# Description:
#      This function checks for synonyms between the gene symbols in the gene sets and the gene expression matrix.
#      If a gene symbol in the gene set is not present in the expression matrix, it attempts to find a matching symbol
#      using a gene symbol conversion table (e.g., `alias2SymbolTable`). The gene set is then updated with the matching symbols.
####################################################


anotate.gmt <- function(genesets, expression){

  # Convert the gene set object into a list (if not already in list form)
  genesets <- geneIds(genesets)

  # Obtain a unique list of all genes found in the gene sets
genes_geneset <- unique(unlist(lapply(genesets,function(x){return(x)})))

  # Identify genes from the expression matrix that are not present in the gene sets
  genes_expression <- rownames(expression)[!rownames(expression) %in% genes_geneset]
  
# Identify genes in the gene set that are not present in the expression matrix
genes <- genes_geneset[!genes_geneset %in% rownames(expression)]

  # Convert gene symbols using a synonym table (alias2SymbolTable)
symbol_genes <- alias2SymbolTable(genes)

  # Create dataframes to map genes from the gene set to the expression matrix
df_geneset <- data.frame("Geneset" = genes, "Symbol" = symbol_genes)

df_expression <- data.frame("Expression" = genes_expression,
"Symbol" = alias2SymbolTable(genes_expression))

# Merge gene set data with expression data to align symbols
df <- merge(df_geneset, df_expression, by="Symbol")
geneset_complete <- data.frame("Geneset" = genes_geneset[!genes_geneset %in% df$Geneset])

df <- merge(df,geneset_complete, by = "Geneset", all = T)
df <- df[!duplicated(df$Geneset),]
rownames(df) <- df[,1]

  # Apply changes in gene symbols to the gene sets
  genesets_new <- lapply(genesets, function(x){
dft <- df[x,]
x[complete.cases(dft$Expression)] <- dft$Expression[complete.cases(dft$Expression)]
return(x)
})

names(genesets_new) <- names(genesets)
return(genesets)
}
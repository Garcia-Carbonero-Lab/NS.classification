# Function: clinical.table
# Input:
#      data: data.frame with samples as rows and features as columns.
#      features: vector of column names (features) to analyze.
#      group: column name representing the group variable.
#      na: boolean indicating whether to include NA samples (TRUE) or exclude them (FALSE).
#
# Output: data.frame indicating the percentage of each level of features in each group and the p-value from Fisher's exact test.
#
# Description: This function evaluates the distribution of categorical features across different group levels and calculates Fisher's exact test p-values.
clinical.table <- function(data, features, group, na = F) {

    # Apply the function to each feature
    dfl <- lapply(features, function(feature) {

        # Create contingency table for the feature and group
        if (na == F) {
            count <- table(data[, feature], data[, group])
        } else {
            count <- table(data[, feature], data[, group], useNA = "a")
        }

        # Perform Fisher's exact test
        f <- try(fisher.test(count, simulate.p.value = T))

        # Handle errors in Fisher's test
        if (class(f) == 'try-error') {
            p <- 'Not calculated'
        } else {
            p <- round(f$p.value, 5)
        }

        # Generate a total count table for the feature
        if (na == F) {
            total <- table(data[, feature])
        } else {
            total <- table(data[, feature], useNA = "a")
        }

        # Combine count and total, then calculate percentages
        count <- cbind(count, total)
        pc <- round(t(t(count) / apply(count, 2, sum)) * 100, 2)

        # Prepare output data frame with counts and percentages
        df <- do.call(rbind, lapply(seq_len(nrow(count)), function(x) {
            paste0(count[x, ], " (", pc[x, ], "%)")
        }))

        # Add p-value to the result
        df <- cbind(df, c(p, rep('', nrow(df) - 1)))

        # Assign column names with group information and total counts
        colnames(df) <- c(paste0(colnames(count), " (", c(table(data[, group]), sum(nrow(data))), ")"), 'P.value')
        rownames(df) <- paste(feature, rownames(count), sep = ".")

        return(df)
    })

    # Combine results from all features into a single data frame
    table <- as.data.frame(do.call(rbind, dfl))
    
    return(table)
}



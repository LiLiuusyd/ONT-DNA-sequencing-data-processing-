# Define input and output file paths
input_file <- "/srv/scratch/z5395183/project_Ga_pUc19/analysis/real_dinucleotides_pos_neg_counts.csv"
output_file <- "/srv/scratch/z5395183/project_Ga_pUc19/analysis/real_dinucleotides_pos_neg_counts_pairwise_comparisons.csv"

# Load data (CSV)
df <- read.csv(input_file, header = TRUE)

# Identify all condition columns (exclude the dinucleotide column)
condition_cols <- setdiff(colnames(df), "dinucleotide")

# Convert all condition columns to numeric
df[condition_cols] <- lapply(df[condition_cols], as.numeric)

# Prepare a list to store results for each pair
results_list <- list()

# Loop over all unique pairs of conditions
pair_idx <- 1
for (i in 1:(length(condition_cols) - 1)) {
  for (j in (i + 1):length(condition_cols)) {
    cond_A <- condition_cols[i]
    cond_B <- condition_cols[j]
    n_A <- sum(df[[cond_A]])
    n_B <- sum(df[[cond_B]])
    
    # Chi-squared test for each dinucleotide
    pvals <- sapply(1:nrow(df), function(k) {
      present_A <- df[[cond_A]][k]
      absent_A  <- n_A - present_A
      present_B <- df[[cond_B]][k]
      absent_B  <- n_B - present_B
      mtx <- matrix(c(present_A, absent_A, present_B, absent_B), nrow = 2, byrow = TRUE)
      chisq.test(mtx)$p.value
    })
    
    # Adjust for multiple testing
    adj_pvals <- p.adjust(pvals, method = "BH")
    
    # Frequencies
    freq_A <- df[[cond_A]] / n_A
    freq_B <- df[[cond_B]] / n_B
    
    # Significance
    significant <- adj_pvals < 0.05
    
    # Store results in a data frame
    res <- data.frame(
      dinucleotide = df$dinucleotide,
      cond_A = cond_A,
      cond_B = cond_B,
      count_A = df[[cond_A]],
      count_B = df[[cond_B]],
      freq_A = freq_A,
      freq_B = freq_B,
      p_value = pvals,
      adj_p_value = adj_pvals,
      significant = significant
    )
    res <- res[order(res$adj_p_value), ]
    results_list[[pair_idx]] <- res
    pair_idx <- pair_idx + 1
  }
}

# Combine all results into one data frame
all_results <- do.call(rbind, results_list)

# View results
print(all_results)

# Write to output file as CSV
write.csv(all_results, output_file, row.names = FALSE, quote = TRUE)
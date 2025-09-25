# Install Biostrings if not already installed
if (!requireNamespace("Biostrings", quietly = TRUE)) {
    install.packages("BiocManager")
    BiocManager::install("Biostrings")
}

library(Biostrings)

# Read the pUC19 sequence from FASTA
seq <- readDNAStringSet("/srv/scratch/z5395183/ProjectGa/analysis/pUC19.fasta")
seq <- as.character(seq[[1]])

# Get all 3-mers
trimers <- substring(seq, 1:(nchar(seq)-2), 3:nchar(seq))

# Define the 3-mers of interest
target_trimers <- c("TTT", "AAA", "CCC", "GGG")

# Count each 3-mer
trimer_table <- table(trimers)

# Make sure all 4 target trimers are present
trimer_counts <- setNames(rep(0, length(target_trimers)), target_trimers)
trimer_counts[names(trimer_table)[names(trimer_table) %in% target_trimers]] <- as.integer(trimer_table[names(trimer_table) %in% target_trimers])

# Calculate frequencies
total_trimers <- sum(trimer_counts)
trimer_freqs <- trimer_counts / sum(trimer_counts)

# Print results
result <- data.frame(
    Trimer = names(trimer_counts),
    Count = as.integer(trimer_counts),
    Frequency = as.numeric(trimer_freqs)
)
print(result)

# 保存为CSV文件
write.csv(result, file = "/srv/scratch/z5395183/ProjectGa/analysis/pUC19_trimer_freq.csv", row.names = FALSE)
# Install Biostrings if not already installed
if (!requireNamespace("Biostrings", quietly = TRUE)) {
    install.packages("BiocManager")
    BiocManager::install("Biostrings")
}

library(Biostrings)

# Read the pUC19 sequence from FASTA
seq <- readDNAStringSet("/srv/scratch/z5395183/ProjectGa/analysis/pUC19.fasta")
seq <- as.character(seq[[1]])

# Get all 5-mers
pentamers <- substring(seq, 1:(nchar(seq)-4), 5:nchar(seq))

# Define the 5-mers of interest
target_pentamers <- c("TTTTT", "AAAAA", "CCCCC", "GGGGG")

# Count each 5-mer
pentamer_table <- table(pentamers)

# Make sure all 4 target pentamers are present
pentamer_counts <- setNames(rep(0, length(target_pentamers)), target_pentamers)
pentamer_counts[names(pentamer_table)[names(pentamer_table) %in% target_pentamers]] <- as.integer(pentamer_table[names(pentamer_table) %in% target_pentamers])

# Calculate frequencies
total_pentamers <- sum(pentamer_counts)
pentamer_freqs <- pentamer_counts / sum(pentamer_counts)

# Print results
result <- data.frame(
    Pentamer = names(pentamer_counts),
    Count = as.integer(pentamer_counts),
    Frequency = as.numeric(pentamer_freqs)
)
print(result)

# 保存为CSV文件
write.csv(result, file = "/srv/scratch/z5395183/ProjectGa/analysis/pUC19_pentamer_freq.csv", row.names = FALSE)
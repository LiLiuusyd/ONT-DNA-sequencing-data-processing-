# Install Biostrings if not already installed
if (!requireNamespace("Biostrings", quietly = TRUE)) {
    install.packages("BiocManager")
    BiocManager::install("Biostrings")
}

library(Biostrings)

# Read the pUC19 sequence from FASTA
seq <- readDNAStringSet("/srv/scratch/z5395183/ProjectGa/analysis/pUC19.fasta")
seq <- as.character(seq[[1]])

# Get all 4-mers
quadmers <- substring(seq, 1:(nchar(seq)-3), 4:nchar(seq))

# Define the 4-mers of interest
target_quadmers <- c("TTTT", "AAAA", "CCCC", "GGGG")

# Count each 4-mer
quadmer_table <- table(quadmers)

# Make sure all 4 target quadmers are present
quadmer_counts <- setNames(rep(0, length(target_quadmers)), target_quadmers)
quadmer_counts[names(quadmer_table)[names(quadmer_table) %in% target_quadmers]] <- as.integer(quadmer_table[names(quadmer_table) %in% target_quadmers])

# Calculate frequencies
total_quadmers <- sum(quadmer_counts)
quadmer_freqs <- quadmer_counts / sum(quadmer_counts)

# Print results
result <- data.frame(
    Quadmer = names(quadmer_counts),
    Count = as.integer(quadmer_counts),
    Frequency = as.numeric(quadmer_freqs)
)
print(result)

# 保存为CSV文件
write.csv(result, file = "/srv/scratch/z5395183/ProjectGa/analysis/pUC19_quadmer_freq.csv", row.names = FALSE)
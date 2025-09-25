# Install Biostrings if not already installed
if (!requireNamespace("Biostrings", quietly = TRUE)) {
    install.packages("BiocManager")
    BiocManager::install("Biostrings")
}

library(Biostrings)

# Read the pUC19 sequence from FASTA
seq <- readDNAStringSet("/srv/scratch/z5395183/ProjectGa/analysis/pUC19.fasta")
seq <- as.character(seq[[1]])

# Get all dinucleotides
dinucs <- substring(seq, 1:(nchar(seq)-1), 2:nchar(seq))

# Count each dinucleotide
dinuc_table <- table(dinucs)

# Make sure all 16 dinucleotides are present
all_dinucs <- c("AA","AC","AG","AT","CA","CC","CG","CT",
                "GA","GC","GG","GT","TA","TC","TG","TT")
dinuc_counts <- setNames(rep(0, length(all_dinucs)), all_dinucs)
dinuc_counts[names(dinuc_table)] <- as.integer(dinuc_table)

# Calculate frequencies
total_dinucs <- sum(dinuc_counts)
dinuc_freqs <- dinuc_counts / total_dinucs

# Print results
result <- data.frame(
    Dinucleotide = names(dinuc_counts),
    Count = as.integer(dinuc_counts),
    Frequency = as.numeric(dinuc_freqs)
)
print(result)


write.csv(result, file = "/srv/scratch/z5395183/ProjectGa/analysis/pUC19_dinucleotide_freq.csv", row.names = FALSE)
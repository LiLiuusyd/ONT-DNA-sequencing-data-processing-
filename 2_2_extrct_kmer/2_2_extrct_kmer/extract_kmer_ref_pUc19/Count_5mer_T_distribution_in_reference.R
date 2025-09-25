# Load the Biostrings package
library(Biostrings)
library(ggplot2)

print(sessionInfo())

# Set the file path to pUC19.fasta (adjust if necessary)
fasta_file <- "/srv/scratch/z5395183/ProjectGa/analysis/pUC19.fasta"

# Read the FASTA file
seq <- readDNAStringSet(fasta_file)
puc19_seq <- seq[[1]]  # Keep as DNAString

# Convert to character string for length calculation
puc19_seq_char <- as.character(puc19_seq)

# Check sequence length (should be 2686 bp for pUC19)
seq_length <- nchar(puc19_seq_char)
cat("Sequence length:", seq_length, "\n")

# Define window size and step
window_size <- 500
step <- 2

# Prepare sliding windows
start_positions <- seq(1, seq_length - window_size + 1, by=step)

# Initialize result storage
data_patterns <- data.frame(start = start_positions,
                           T = integer(length(start_positions)),
                           TT = integer(length(start_positions)),
                           TTT = integer(length(start_positions)),
                           TTTT = integer(length(start_positions)),
                           TTTTT = integer(length(start_positions)))

# Count patterns in each window
for (i in seq_along(start_positions)) {
  window_seq <- subseq(puc19_seq, start=start_positions[i], width=window_size)
  data_patterns$T[i]      <- countPattern("T", window_seq, fixed=TRUE)
  data_patterns$TT[i]     <- countPattern("TT", window_seq, fixed=TRUE)
  data_patterns$TTT[i]    <- countPattern("TTT", window_seq, fixed=TRUE)
  data_patterns$TTTT[i]   <- countPattern("TTTT", window_seq, fixed=TRUE)
  data_patterns$TTTTT[i]  <- countPattern("TTTTT", window_seq, fixed=TRUE)
}

# Save results
write.csv(data_patterns, "/srv/scratch/z5395183/ProjectGa/analysis/t_distribution_patterns.csv", row.names = FALSE)



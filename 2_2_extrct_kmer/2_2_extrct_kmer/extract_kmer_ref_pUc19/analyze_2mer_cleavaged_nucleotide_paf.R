# Required packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("Biostrings", quietly = TRUE)) BiocManager::install("Biostrings")
library(Biostrings)

# ---- User Inputs ----
paf_file <- "/srv/scratch/z5395183/ProjectGa/20250506_demux_fastq/demu_SQK-NBD114-96_barcode25.paf"
ref_file <- "/srv/scratch/z5395183/ProjectGa/analysis/pUC19.fasta"
output_file <- "/srv/scratch/z5395183/ProjectGa/20250506_demux_fastq/demu_barcode25_05092000.txt"

# ---- Read Reference ----
ref_seq <- readDNAStringSet(ref_file)
if (length(ref_seq) != 1) stop("Reference FASTA must contain one sequence")
ref <- as.character(ref_seq[[1]])
ref_len <- nchar(ref)

# Read PAF file without specifying col.names
paf <- read.table(paf_file, sep="\t", header=FALSE, fill=TRUE, comment.char = "", quote = "")

# Assign names to the first 12 columns only
if (ncol(paf) < 12) stop("PAF file has fewer than 12 columns.")
col_names <- c("query_name", "query_length", "query_start", "query_end",
               "strand", "target_name", "target_length", "target_start", "target_end",
               "residue_matches", "alignment_block_length", "mapping_quality")
names(paf)[1:12] <- col_names

# Now use the first 12 columns as before
start_pos <- paf$target_start + 1  # 0-based to 1-based
end_pos <- paf$target_end          # end is exclusive, so this is the first base after the alignment

# ---- Filter by boundaries (with exclusions) ----
valid_start_idx <- (start_pos >= 20 & start_pos <= 2666 & start_pos != 234)
valid_end_idx <- (end_pos >= 20 & end_pos <= 2666 & end_pos != 237)

# Get strand info
strands <- paf$strand

# Separate indices for positive and negative strands
pos_idx <- strands == "+"
neg_idx <- strands == "-"

# Valid starts and ends for each strand
valid_starts_pos <- start_pos[valid_start_idx & pos_idx]
valid_ends_pos <- end_pos[valid_end_idx & pos_idx]
valid_starts_neg <- start_pos[valid_start_idx & neg_idx]
valid_ends_neg <- end_pos[valid_end_idx & neg_idx]

# ---- Dinucleotide Extraction ----
get_dinucs <- function(start, end, ref) {
  ref_len <- nchar(ref)
  before_start <- ifelse(start > 1, 
                         substr(ref, start - 1, start), 
                         paste0("N", substr(ref, start, start)))
  after_end <- ifelse(end <= ref_len, 
                      substr(ref, end, min(end + 1, ref_len)), 
                      paste0(substr(ref, ref_len, ref_len), "N"))
  c(before_start, after_end)
}

dinucs_start_pos <- sapply(valid_starts_pos, function(start) get_dinucs(start, start, ref)[1])
dinucs_end_pos <- sapply(valid_ends_pos, function(end) get_dinucs(end, end, ref)[2])
dinucs_pos <- c(dinucs_start_pos, dinucs_end_pos)
dinucs_pos <- dinucs_pos[!grepl("N", dinucs_pos)]

dinucs_start_neg <- sapply(valid_starts_neg, function(start) get_dinucs(start, start, ref)[1])
dinucs_end_neg <- sapply(valid_ends_neg, function(end) get_dinucs(end, end, ref)[2])
dinucs_neg <- c(dinucs_start_neg, dinucs_end_neg)
dinucs_neg <- dinucs_neg[!grepl("N", dinucs_neg)]

# ---- Count Frequencies ----
dinuc_all <- c("AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT",
               "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT")
freq_pos <- table(factor(dinucs_pos, levels = dinuc_all))
total_pos <- sum(freq_pos)
prop_pos <- freq_pos / total_pos

freq_neg <- table(factor(dinucs_neg, levels = dinuc_all))
total_neg <- sum(freq_neg)
prop_neg <- freq_neg / total_neg

# ---- Output ----
cat("Total cleavage dinucleotides extracted (positive strand):", total_pos, "\n")
cat("\nFrequency table for all dinucleotides (positive strand):\n")
print(cbind(Frequency = freq_pos, Proportion = round(prop_pos, 4)))

cat("\nTotal cleavage dinucleotides extracted (negative strand):", total_neg, "\n")
cat("\nFrequency table for all dinucleotides (negative strand):\n")
print(cbind(Frequency = freq_neg, Proportion = round(prop_neg, 4)))

# Save to file
sink(output_file)
cat("Total cleavage dinucleotides extracted (positive strand):", total_pos, "\n")
cat("\nFrequency table for all dinucleotides (positive strand):\n")
print(cbind(Frequency = freq_pos, Proportion = round(prop_pos, 4)))

cat("\nTotal cleavage dinucleotides extracted (negative strand):", total_neg, "\n")
cat("\nFrequency table for all dinucleotides (negative strand):\n")
print(cbind(Frequency = freq_neg, Proportion = round(prop_neg, 4)))
sink()
cat("Results saved to:", output_file, "\n")
# Required packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("Biostrings", quietly = TRUE)) BiocManager::install("Biostrings")
library(Biostrings)

# ---- User Inputs ----
paf_file <- "/srv/scratch/z5395183/ProjectGa/20250506_demux_fastq/demu_SQK-NBD114-96_barcode26.paf"
ref_file <- "/srv/scratch/z5395183/ProjectGa/analysis/pUC19.fasta"
output_file <- "/srv/scratch/z5395183/ProjectGa/20250506_demux_fastq/demu_barcode26_3mer.csv"

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
valid_start_idx <- (start_pos >= 20 & start_pos <= 2666)
valid_end_idx <- (end_pos >= 20 & end_pos <= 2666)

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

# ---- 3-mer Extraction (all four windows) ----
get_trimer_window <- function(pos, ref, offset_start, offset_end) {
  ref_len <- nchar(ref)
  start <- pos + offset_start
  end <- pos + offset_end
  if (start < 1 || end > ref_len) {
    return(NA)
  }
  substr(ref, start, end)
}

# Define the four windows as offsets relative to the cleavage site
# (offset_start, offset_end)
trimer_windows <- list(
  window1 = c(0, 2),    # start to start+2
  window2 = c(-1, 1),   # start-1 to start+1
  window3 = c(-2, 0),   # start-2 to start
  window4 = c(-3, -1)   # start-3 to start-1
)

# Helper to extract all windows for a set of positions
extract_trimers_all_windows <- function(positions, ref) {
  lapply(trimer_windows, function(offsets) {
    sapply(positions, function(pos) get_trimer_window(pos, ref, offsets[1], offsets[2]))
  })
}

# For positive strand
start_trimers_pos_all <- extract_trimers_all_windows(valid_starts_pos, ref)
end_trimers_pos_all   <- extract_trimers_all_windows(valid_ends_pos, ref)
# For negative strand
start_trimers_neg_all <- extract_trimers_all_windows(valid_starts_neg, ref)
end_trimers_neg_all   <- extract_trimers_all_windows(valid_ends_neg, ref)

# Only count AAA, TTT, CCC, GGG
bases <- c("A", "T", "C", "G")
trimer_targets <- as.vector(outer(outer(bases, bases, paste0), bases, paste0))
count_trimers <- function(trimers, targets) {
  trimers <- trimers[!is.na(trimers)]
  tab <- table(factor(trimers, levels = targets))
  total <- sum(tab)
  prop <- if (total > 0) tab / total else rep(0, length(targets))
  list(counts = tab, proportions = prop, total = total)
}

# Count for each window, each position/strand
results <- list(
  pos_start = lapply(start_trimers_pos_all, count_trimers, targets = trimer_targets),
  pos_end   = lapply(end_trimers_pos_all, count_trimers, targets = trimer_targets),
  neg_start = lapply(start_trimers_neg_all, count_trimers, targets = trimer_targets),
  neg_end   = lapply(end_trimers_neg_all, count_trimers, targets = trimer_targets)
)

# ---- Output ----
window_names <- c(
  "Window 1: pos to pos+2",
  "Window 2: pos-1 to pos+1",
  "Window 3: pos-2 to pos",
  "Window 4: pos-3 to pos-1"
)
strand_pos_names <- c("Positive strand, start", "Positive strand, end", "Negative strand, start", "Negative strand, end")

cat("\n================ Cleavage 3-mer Extraction Results ================\n")
for (strand_pos in names(results)) {
  cat("\n---", gsub("_", ", ", strand_pos), "---\n")
  for (i in seq_along(window_names)) {
    res <- results[[strand_pos]][[i]]
    cat(window_names[i], "\n")
    cat("Total 3-mers:", res$total, "\n")
    print(data.frame(Trimer = trimer_targets, Frequency = as.integer(res$counts), Proportion = round(res$proportions, 4)))
    cat("\n")
  }
}

# ---- Summarize by strand (sum all 4 windows for each strand) ----
sum_trimer_counts <- function(res_list) {
  # res_list: list of 4 window results for a strand (each is a list with counts)
  Reduce(`+`, lapply(res_list, function(x) as.integer(x$counts)))
}

# Sum for positive and negative strands
pos_total_counts <- sum_trimer_counts(c(results$pos_start, results$pos_end))
neg_total_counts <- sum_trimer_counts(c(results$neg_start, results$neg_end))

# Prepare summary data frame
summary_df <- data.frame(
  Trimer = trimer_targets,
  Positive_Strand = pos_total_counts,
  Negative_Strand = neg_total_counts
)

# ---- Output to CSV ----
# Prepare detailed results for all windows/strands
library(dplyr)
all_results <- list()
for (strand_pos in names(results)) {
  for (i in seq_along(window_names)) {
    res <- results[[strand_pos]][[i]]
    df <- data.frame(
      Strand_Pos = gsub("_", ", ", strand_pos),
      Window = window_names[i],
      Trimer = trimer_targets,
      Frequency = as.integer(res$counts),
      Proportion = round(res$proportions, 4)
    )
    all_results[[paste(strand_pos, i, sep = "_")]] <- df
  }
}
detailed_df <- bind_rows(all_results)

# Write both summary and detailed results to CSV
summary_file <- sub(".csv$", "_summary.csv", output_file)
write.csv(summary_df, summary_file, row.names = FALSE)
write.csv(detailed_df, output_file, row.names = FALSE)

cat("Summary results saved to:", summary_file, "\n")
cat("Detailed results saved to:", output_file, "\n")
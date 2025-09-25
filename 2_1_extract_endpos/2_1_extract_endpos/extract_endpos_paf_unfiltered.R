############################################################
# 1) Install/load required packages
############################################################

if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages("data.table")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  install.packages("ggrepel")
}
### REVISED ### 添加 openxlsx 包以写入 Excel 文件
if (!requireNamespace("openxlsx", quietly = TRUE)) {
  install.packages("openxlsx")
}

library(data.table)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(openxlsx)  # 用于保存 Excel 文件
###

print(sessionInfo())

############################################################
# 2) Define file paths and parameters (WSL-style)
############################################################

paf_file <- "/srv/scratch/z5395183/ProjectGa/20250506_demux_fastq/demu_SQK-NBD114-96_barcode27.paf"
output_plot <- "/srv/scratch/z5395183/ProjectGa/20250506_demux_fastq/demux_27.png"
### REVISED ### 添加 Excel 输出路径
output_excel <- "/srv/scratch/z5395183/ProjectGa/20250506_demux_fastq/demux_27.xlsx"
###
puc19_len <- 2686
MAX_Y <- 171000

if (!file.exists(paf_file)) {
  stop("Input PAF file not found at: ", paf_file)
}

############################################################
# 3) Read the PAF file and extract endpoints
############################################################

paf <- fread(paf_file, fill = TRUE, header = FALSE)
cat("PAF file loaded, rows:", nrow(paf), "\n")
if (ncol(paf) < 9) {
  stop("PAF file has fewer than 9 columns, check file format.")
}

start_positions <- paf[[8]] + 1
stop_positions  <- paf[[9]]

############################################################
# 4) Count start and end endpoints at each base along pUC19
############################################################

start_counts <- integer(puc19_len)
end_counts   <- integer(puc19_len)

for (pos in start_positions) {
  if (pos >= 1 && pos <= puc19_len) {
    start_counts[pos] <- start_counts[pos] + 1
  }
}
for (pos in stop_positions) {
  if (pos >= 1 && pos <= puc19_len) {
    end_counts[pos] <- end_counts[pos] + 1
  }
}

df_endpoints <- data.frame(
  Position = 1:puc19_len,
  StartCount = start_counts,
  EndCount = end_counts
)

############################################################
# 5) Filter out positions 1 and 2686 (boundary)
############################################################

boundary_left <- 1
boundary_right <- 2686
df_filtered <- df_endpoints %>%
  dplyr::filter(Position > boundary_left & Position < boundary_right)
cat("Filtered data, rows:", nrow(df_filtered), "\n")

# Calculate mean, sd, and threshold for start and end
mean_start <- mean(df_filtered$StartCount)
sd_start   <- sd(df_filtered$StartCount)
threshold_start <- mean_start + 1 * sd_start

mean_end <- mean(df_filtered$EndCount)
sd_end   <- sd(df_filtered$EndCount)
threshold_end <- mean_end + 1 * sd_end

cat("Start: Mean =", mean_start, ", SD =", sd_start, ", Threshold =", threshold_start, "\n")
cat("End:   Mean =", mean_end,   ", SD =", sd_end,   ", Threshold =", threshold_end, "\n")

# Save all start and end sites within the boundary to separate sheets in Excel

all_starts <- df_filtered %>%
  select(Position, StartCount) %>%
  rename(X = Position, Y = StartCount)

all_ends <- df_filtered %>%
  select(Position, EndCount) %>%
  rename(X = Position, Y = EndCount)

write.xlsx(
  list(AllStarts = all_starts, AllEnds = all_ends),
  file = output_excel,
  rowNames = FALSE
)
cat("All start and end points saved to: ", output_excel, "\n")

############################################################


cat("Script completed successfully!\n")
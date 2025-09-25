# Install Dependencies for ONT DNA Sequencing Data Processing

# Check if R is installed
if (!requireNamespace("utils", quietly = TRUE)) {
  stop("R is not properly installed. Please install R first.")
}

# Function to install packages if not already installed
install_if_missing <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  if(length(new_packages)) {
    install.packages(new_packages, dependencies = TRUE)
  }
}

# Function to install Bioconductor packages
install_bioc_if_missing <- function(packages) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  if(length(new_packages)) {
    BiocManager::install(new_packages, dependencies = TRUE)
  }
}

cat("Installing R dependencies for ONT DNA sequencing data processing...\n")

# Core R packages
core_packages <- c(
  "ggplot2", "dplyr", "tidyr", "svglite", "readr", "stringr",
  "data.table", "reshape2", "gridExtra", "cowplot", "viridis"
)

cat("Installing core R packages...\n")
install_if_missing(core_packages)

# Bioconductor packages
bioc_packages <- c(
  "Rsamtools", "GenomicAlignments", "GenomicRanges", "IRanges",
  "Biostrings", "ShortRead", "rtracklayer"
)

cat("Installing Bioconductor packages...\n")
install_bioc_if_missing(bioc_packages)

# Verify installations
cat("\nVerifying package installations...\n")
all_packages <- c(core_packages, bioc_packages)
missing_packages <- all_packages[!(all_packages %in% installed.packages()[,"Package"])]

if(length(missing_packages) == 0) {
  cat(" All packages installed successfully!\n")
} else {
  cat(" Some packages failed to install:\n")
  cat(paste(missing_packages, collapse = ", "), "\n")
  cat("Please install these packages manually.\n")
}

cat("\nInstallation complete!\n")
cat("You can now run the analysis pipeline.\n")

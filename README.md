# ONT-DNA-sequencing-data-processing-

This repository contains all code and data necessary to reproduce the ONT DNA sequencing results presented in the paper titled **"Gallium in Liquid State Shows Nuclease-Mimicking Activity"**.

## Research Overview

This work investigates the novel finding that gallium (Ga) in liquid state exhibits nuclease-mimicking activity, capable of cleaving DNA. The analysis pipeline processes Oxford Nanopore Technology (ONT) sequencing data to identify and characterize DNA cleavage patterns when DNA is exposed to Ga droplets, using pUC19 plasmid as a model system.

## Key Findings

- Gallium droplets show sequence-specific DNA cleavage activity
- Cleavage patterns can be analyzed through k-mer analysis
  - Dinucleotide (2-mer) analysis revealed significant cleavage preferences - these results are presented in the main text and supplementary information
  - 3-mer, 4-mer, and 5-mer analyses were performed but did not show significant differences, so these results are not included in the paper but the script is provided here for reference
- Statistical analysis reveals significant cleavage preferences
- The method provides insights into analyzing the cleavage patterns of metal-nanoparticle-based nuclease activity

## Repository Structure
├── README.md
├── LICENSE
├── .gitignore
├── CITATION.cff
├── CONTRIBUTING.md
├── install_dependencies.R
├── environment.yml
├── setup.sh
├── 1_demu_align_bam_baf/
│ └── 1_demu_align_bam_baf.R
├── 2_1_extract_endpos/
│ ├── extract_endpos_paf_filtered_sd.R
│ └── extract_endpos_paf_unfiltered.R
├── 2_2_extrct_kmer/
│ ├── extract_kmer_ref_pUc19/
│ │ ├── Extract_2nucleotide_from_pUc19.R
│ │ ├── Extract_3nucleotide_from_pUc19.R
│ │ ├── Extract_4nucleotide_from_pUc19.R
│ │ ├── Extract_5nucleotide_from_pUc19.R
│ │ ├── Count_5mer_T_distribution_in_reference.R
│ │ └── analyze_2mer_cleavaged_nucleotide_paf.R
│ └── extract_kmer_tested_samples/
│ ├── README.txt
│ ├── analyze_2mer_cleavaged_nucleotide_paf.R
│ ├── analyze_3mer_cleavaged_nucleotide_paf.R
│ └── analyze_4mer_cleavaged_nucleotide_paf.R
└── 3_chi_squared_test/
├── Chi_squaredtest.R
└── 2mer/



## Requirements

### Software Dependencies
- **R** (version 4.0 or higher)
- **Linux environment** (recommended for bioinformatics tools)
- **Bioinformatics tools:**
  - `minimap2` (for sequence alignment)
  - `samtools` (for BAM file processing)
  - `dorado` (for nanopore basecalling and demultiplexing)
  - `htsbox` (for BAM to PAF conversion)

### R Packages
```r
# Core packages
install.packages(c("ggplot2", "dplyr", "tidyr", "svglite", "readr", "stringr"))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("Rsamtools", "GenomicAlignments", "GenomicRanges", "IRanges"))
```

## Installation

### 1. Clone the repository
```bash
git clone https://github.com/LiLiuusyd/ONT-DNA-sequencing-data-processing-.git
cd ONT-DNA-sequencing-data-processing-
```

### 2. Install R dependencies
```bash
Rscript install_dependencies.R
```

### 3. Install bioinformatics tools
```bash
# Using conda (recommended)
conda install -c bioconda minimap2 samtools

# Install dorado (Nanopore basecaller)
# Follow instructions at: https://github.com/nanoporetech/dorado
```
## Data Availability
 All raw sequencing data are available at the NCBI Sequence Read Archive:
    - BioProject: [PRJNA1307762](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1307762)
    - Runs: SRR35033888, SRR35033889, SRR35033890, SRR35033891, SRR35033892

## Usage


### Pipeline Overview
The analysis follows this workflow:
1. **Demultiplexing and Alignment** → Process raw FASTQ files
2. **K-mer Extraction** → Extract sequence patterns from reference and test samples
3. **Statistical Analysis** → Perform chi-squared tests on cleavage patterns
4. **Visualization** → Generate coverage plots and cleavage pattern figures, it is higly recommended to processing the data with the code in this repository but plot it on your own.

### Running the Analysis

#### Step 1: Demultiplexing and Alignment
```bash
cd 1_demu_align_bam_baf/
Rscript 1_demu_align_bam_baf.R
```

#### Step 2: Extract K-mers from Reference
```bash
cd 2_extrct_kmer/extract_kmer_ref_pUc19/
Rscript Extract_2nucleotide_from_pUc19.R
Rscript Extract_3nucleotide_from_pUc19.R
Rscript Extract_4nucleotide_from_pUc19.R
Rscript Extract_5nucleotide_from_pUc19.R
```

#### Step 3: Analyze K-mers in Test Samples
```bash
cd 2_extrct_kmer/extract_kmer_tested_samples/
Rscript analyze_2mer_cleavaged_nucleotide_paf.R
Rscript analyze_3mer_cleavaged_nucleotide_paf.R
Rscript analyze_4mer_cleavaged_nucleotide_paf.R
```

#### Step 4: Statistical Analysis
```bash
cd 3_chi_squared_test/
Rscript Chi_squaredtest.R
```

## Results and Paper Integration

### Published Results
- **Dinucleotide (2-mer) analysis**: Results showing significant cleavage preferences are presented in the main text and supplementary information
- **Statistical significance**: Chi-squared tests revealed significant differences in dinucleotide cleavage patterns

### Additional Analyses (Not Published)
- **3-mer, 4-mer, and 5-mer analyses**: These were performed but did not show significant differences, so these results are not included in the paper
- **Complete analysis scripts**: All k-mer analysis scripts (2-mer through 5-mer) are provided for reproducibility and potential future analysis

## Data Requirements

### Input Data
- **FASTQ files**: ONT sequencing data of DNA exposed to Ga droplets
- **Reference genome**: pUC19 plasmid sequence (FASTA format)
- **Barcode information**: For demultiplexing 

### Expected Outputs
- **BAM files**: Aligned sequencing data
- **Coverage files**: Read coverage across the plasmid
- **K-mer analysis**: Cleavage site patterns (2-mer, 3-mer, 4-mer, 5-mer)
- **Statistical results**: Chi-squared test results for cleavage preferences

## Reproducibility

### Reproducing Results
1. Install all dependencies (see Installation section)
2. Prepare your sequencing data in the required format
3. Run the pipeline in the specified order
4. Check outputs in the `data/results/` directory

### Path Modification Guide
To adapt the scripts for your environment:

1. **For R scripts**: Replace paths like `/srv/scratch/z5395183/...` with your local paths
2. **For bash commands**: Update directory paths in the main pipeline script
3. **Create necessary directories**: Ensure all required input/output directories exist


## Contact
- **Author**: [Li Liu]
- **Email**: [lliu7637@uni.sydney.edu.au]
- **Institution**: [University of Sydney]
- **ORCID**: [0000-0003-0093-5191]

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments
K.K.-Z. acknowledges support from the Australian Research Council (ARC) Discovery Project (DP240101086). C.Z. is grateful for funding from the Royal Society, UK (Grant Nos. RG/R1/241228 and IEC/NSFC/233339). M.J.S.S acknowledges support from the ARC Discovery Project (DP240101215). This research was supported by the Australian Government’s National Collaborative Research Infrastructure Strategy (NCRIS), with computational resources provided by the National Computational Infrastructure (NCI) Facility and the Pawsey Supercomputing Research Centre through the National Computational Merit Allocation Scheme.

We acknowledge Professor Ewa Magdalena Goldys for providing experiential resources. Technical assistance from Lewis Adler at the Solid State & Elemental Analysis Unit (Mark Wainwright Analytical Centre, UNSW Sydney) is gratefully acknowledged. We also acknowledge the facilities and the scientific and technical assistance of Sydney Analytical, a core research facility at the University of Sydney (USYD). The assistance of supercomputing resources from UNSW Katana is also acknowledged. We thank Matthew Wong, Lydia Murphy, and Aravind Manda at the Ramaciotti Centre for Genomics for their technical assistance with Oxford Nanopore Technology. We also thank Yujian Shi from the Sarina Sarina group at USYD for providing a halogen lamp.


## Data Availability
 All raw sequencing data are available at the NCBI Sequence Read Archive:
    - BioProject: [PRJNA1307762](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1307762)
    - Runs: SRR35033888, SRR35033889, SRR35033890, SRR35033891, SRR35033892


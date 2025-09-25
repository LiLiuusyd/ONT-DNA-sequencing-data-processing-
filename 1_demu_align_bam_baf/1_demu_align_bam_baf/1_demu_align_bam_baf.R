# Set working directory
cd /srv/scratch/z5395183/data_raw_fastq_all

# Create reference directory and copy pUC19.fasta
mkdir -p reference
cp /srv/scratch/z5395183/data_raw_fastq_all/pUC19.fasta reference/
ls reference/

# Demultiplex all fastq files (if needed)
# If dorado demux works on a directory, this is fine:
dorado demux \
  --kit-name SQK-NBD114-96 \
  --output-dir /srv/scratch/z5395183/ProjectGa/1 \
  --barcode-both-ends \
  --emit-fastq \
  /srv/scratch/z5395183/ProjectGa/Allfastq


# Alignment, BAM, Index, and Coverage for all fastq(.gz) files
REF="/srv/scratch/z5395183/data_raw_fastq_all/reference/pUC19.fasta"
for fq in *.fastq *.fastq.gz; do
  [ -e "$fq" ] || continue  # skip if no files match
  base=$(basename "$fq" .fastq.gz)
  if [[ "$base" == "$fq" ]]; then
    base=$(basename "$fq" .fastq)
  fi
  BAM="${base}.bam"
  # Align and sort
  minimap2 -a -t 4 "$REF" "$fq" | \
    samtools view -b -@ 2 - | \
    samtools sort -@ 2 -o "$BAM"
  # Index BAM
  samtools index "$BAM"
  # Coverage
  samtools depth -@ 2 -a "$BAM" > "${base}_coverage.txt"

  
done

# Convert BAMs to PAF format
for file in *.bam; do
  ~/apps/htsbox/htsbox samview -p "$file" > "${file%.bam}.paf"
done

# Get read lengths for all fastq(.gz) files
for file in *.fastq *.fastq.gz; do
  [ -e "$file" ] || continue
  base=$(basename "$file" .fastq.gz)
  if [[ "$base" == "$file" ]]; then
    base=$(basename "$file" .fastq)
  fi
  if [[ "$file" == *.gz ]]; then
    zcat "$file" | awk -v f="$base" '{ if (NR % 4 == 2 ) print f"\t"length($1)}'
  else
    awk -v f="$base" '{ if (NR % 4 == 2 ) print f"\t"length($1)}' "$file"
  fi
done > merged_readLens.txt

# Plot read lengths in R
Rscript -e "
library(ggplot2)
d <- read.delim('merged_readLens.txt', header=FALSE)
ggplot(d, aes(V2, colour=factor(V1))) + geom_density(adjust=0.5) + scale_x_continuous(limits=c(500,4000))
"

# Plot mapped reads extremities in R (assuming mapped_coords_merged.tsv exists)
Rscript -e "
library(ggplot2)
d <- read.delim('mapped_coords_merged.tsv', header=FALSE)
ggplot(d, aes(V3)) + geom_histogram(binwidth=25) + facet_wrap(~factor(V1), ncol=1) + scale_y_continuous(limits=c(0,100)) + xlab('Plasmid coordinates')
"
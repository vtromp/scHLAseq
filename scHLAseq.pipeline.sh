#!/usr/bin/env bash

# -------------------------------
# Initialize variables
# -------------------------------

# Empty variables to store user-provided paths
FASTQS=""
READ_LENGTH=150
REF_TRANSCRIPTOME=""
OUTPUT_DIR=""
TEMP_DIR=""


# -------------------------------
# Define a help/usage function
# -------------------------------

# Explain each expected command-line option
usage() {
  echo "Usage: $0 -f <FASTQ_DIR> -l <READ_LENGTH> -r <REF_TRANSCRIPTOME_DIR> -o <OUTPUT_DIR> -t <TEMP_DIR>"
  echo
  echo "  -f FASTQ_DIR              Directory containing sequencing FASTQ files"
  echo "  -l READ_LENGTH            Read length (default: 150)"
  echo "  -r REF_TRANSCRIPTOME_DIR  Path to CellRanger-compatible reference transcriptome"
  echo "  -o OUTPUT_DIR             Output directory"
  echo "  -t TEMP_DIR               Temporary directory"
  echo
  exit 1
}


# -------------------------------
# Parse command-line options using getopts
# -------------------------------

# Parse command-line options
while getopts ":f:r:o:t:h" opt; do
  case $opt in
    f) FASTQS="$OPTARG" ;;
    l) READ_LENGTH="$OPTARG" ;;
    r) REF_TRANSCRIPTOME="$OPTARG" ;;
    o) OUTPUT_DIR="$OPTARG" ;;
    t) TEMP_DIR="$OPTARG" ;;
    h) usage ;;
    \?) echo "Invalid option -$OPTARG"; usage ;;
    :) echo "Option -$OPTARG requires an argument"; usage ;;
  esac
done

# Remove trailing slash (if present)
FASTQS="${FASTQS%/}"
REF_TRANSCRIPTOME="${REF_TRANSCRIPTOME%/}"
OUTPUT_DIR="${OUTPUT_DIR%/}"
TEMP_DIR="${TEMP_DIR%/}"

#
OUTPUT_DIR="$OUTPUT_DIR/$(basename "$FASTQS")_scHLAseq"


# -------------------------------
# Run cellranger count
# ------------------------------- 

# Create output directory for cellranger count
mkdir -p ""$OUTPUT_DIR"/cellranger_count/"

# Run the Cell Ranger count pipeline and write BAM files and all output results to 'cellranger_count' directory in output folder
cellranger count --id cellranger_count \
  --fastqs "$FASTQS" \
  --transcriptome "$REF_TRANSCRIPTOME" \
  --create-bam true \
  --output-dir ""$OUTPUT_DIR"/cellranger_count" \
  --localcores "$(nproc)" \
  --localmem "$(awk '/MemAvailable/ {printf "%.0f\n", $2/1024/1024}' /proc/meminfo)"


# -------------------------------
# Run arcasHLA
# ------------------------------- 

# Create output directory for arcasHLA genotype
mkdir -p ""$OUTPUT_DIR"/arcashla_genotype/"

# Extract HLA-related and unmapped reads from the Cell Ranger BAM file
arcasHLA extract ""$OUTPUT_DIR"/cellranger_count/outs/possorted_genome_bam.bam" \
  --unmapped \
  --single \
  --outdir ""$OUTPUT_DIR"/arcashla_genotype/" \
  --temp "$TEMP_DIR" \
  --threads "$(nproc)"

# Perform HLA genotyping using the extracted reads
arcasHLA genotype ""$OUTPUT_DIR"/arcashla_genotype/extracted_reads.fq.gz" \
  --single \
  --avg "$READ_LENGTH" \
  --outdir ""$OUTPUT_DIR"/arcashla_genotype/" \
  --temp "$TEMP_DIR" \
  --threads "$(nproc)"

# Rename extracted_reads output files to remove the common prefix
for file in ""$OUTPUT_DIR"/arcashla_genotype/extracted_reads.*"; do
  base=$(basename "$file")
  newname=${base#extracted_reads.}
  mv "$file" ""$OUTPUT_DIR"/arcashla_genotype/"$newname""
done

#
fetch.HLA.allele.seqs.R -g ""$OUTPUT_DIR"/arcashla_genotype/genotype.json" -o ""$OUTPUT_DIR"/arcashla_genotype/" --class-II-only


# -------------------------------
# 
# ------------------------------- 

#
zcat ""$REF_TRANSCRIPTOME"/genes/genes.gtf.gz" | \
awk -F '\t' 'BEGIN{OFS="\t"}
$1=="chr6" &&
$3=="gene" &&
$9 ~ /gene_name "HLA-(DP|DQ|DR)/ {
  match($9, /gene_name "([^"]+)"/, g)
  print $1, $4-1, $5, g[1], ".", $7
}' > ""$OUTPUT_DIR"/HLA_classII_genes.bed"

#
samtools view -h -L HLA_classII_genes.bed ""$OUTPUT_DIR"/cellranger_count/outs/possorted_genome_bam.bam" | awk '$0 ~ /^@/ || ($0 ~ /\tCB:Z:/ && $0 ~ /\tUB:Z:/)' > ""$OUTPUT_DIR"/HLA_classII_mapped_reads.sam"
samtools view -h -f 4 ""$OUTPUT_DIR"/cellranger_count/outs/possorted_genome_bam.bam" | awk '$0 ~ /^@/ || ($0 ~ /\tCB:Z:/ && $0 ~ /\tUB:Z:/)' > ""$OUTPUT_DIR"/unmapped_reads.sam"
#
samtools view -b -o ""$OUTPUT_DIR"/HLA_classII_mapped_reads.bam" ""$OUTPUT_DIR"/HLA_classII_mapped_reads.sam"
samtools view -b -o ""$OUTPUT_DIR"/unmapped_reads.bam" ""$OUTPUT_DIR"/unmapped_reads.sam"
#
samtools merge -o ""$OUTPUT_DIR"/HLA_classII_and_unmapped_reads.bam" ""$OUTPUT_DIR"/HLA_classII_mapped_reads.bam" ""$OUTPUT_DIR"/unmapped_reads.bam"
#
rm ""$OUTPUT_DIR"/HLA_classII_mapped_reads.sam" ""$OUTPUT_DIR"/unmapped_reads.sam" ""$OUTPUT_DIR"/HLA_classII_mapped_reads.bam" ""$OUTPUT_DIR"/unmapped_reads.bam"

#
samtools fastq ""$OUTPUT_DIR"/HLA_classII_and_unmapped_reads.bam" > ""$OUTPUT_DIR"/HLA_classII_and_unmapped_reads.fastq"

#
samtools view ""$OUTPUT_DIR"/HLA_classII_and_unmapped_reads.bam" | \
awk '{
  read_id=$1;
  chr=$3;
  pos=$4;
  cb="NA";
  ub="NA";
  for(i=12;i<=NF;i++){
    if($i ~ /^CB:Z:/) { cb=substr($i,6) }
    if($i ~ /^UB:Z:/) { ub=substr($i,6) }
  }
  if(cb != "NA" && ub != "NA") {
    print read_id"\t"cb"\t"ub"\t"chr"\t"pos
  }
}' > ""$OUTPUT_DIR"/overview_reads.tsv"

#
vsearch --usearch_global ""$OUTPUT_DIR"/HLA_classII_and_unmapped_reads.fastq" --db ""$OUTPUT_DIR"/arcashla_genotype/genotype.fasta" \
  --id 0.96 \
  --query_cov 1.0 \
  --maxaccepts 0 \
  --maxrejects 0 \
  --userfields query+target+id+mism+gaps+tilo+tihi \
  --userout ""$OUTPUT_DIR"/all_HLA_alignments.tsv"


# -------------------------------
# 
# ------------------------------- 


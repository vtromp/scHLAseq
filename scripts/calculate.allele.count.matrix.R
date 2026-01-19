#!/usr/bin/env Rscript

# Define the command-line options for the script using optparse
option_list <- list(
  optparse::make_option(opt_str = base::c("-r", "--reads"), action = "store", type = "character", help = "TSV file summarizing reads with cell barcode and UMI information", metavar = "TSV"),
  optparse::make_option(opt_str = base::c("-a", "--alignments"), action = "store", type = "character", help = "TSV file produced by vsearch containing read-to-HLA alignments", metavar = "TSV"),
  optparse::make_option(opt_str = base::c("-g", "--gtf"), action = "store", type = "character", help = "GTF file of reference transcriptome", metavar = "GTF"),
  optparse::make_option(opt_str = base::c("-o", "--output-dir"), action = "store", type = "character", default = base::getwd(), help = "Directory for output (default: working directory)", metavar = "DIR")
)

# Parse command-line arguments based on the defined options
opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

# Assign parsed arguments to internal variables for use in the script
read_file <- opt$`reads`
alignment_file <- opt$`alignments`
gtf_file <- opt$`gtf`
output_dir <- opt$`output-dir`

# Remove trailing slash if present
output_dir <- base::sub(pattern = "/$", replacement = "", x = output_dir)

# Read the vsearch alignment TSV file (no header) and assign column names
alignments <- utils::read.table(file = alignment_file, header = FALSE, sep = "\t")
base::colnames(alignments) <- base::c("read.id", "allele", "match", "n.mismatches", "n.gaps", "start", "end")

# Extract the HLA gene name (everything before the "*")
alignments$hla.gene <- base::sub(pattern = "\\*.*$", replacement = "", alignments$allele)

# Filter alignments to retain only reads that map unambiguously to a single HLA gene (reads with alignments to multiple different HLA genes are removed)
alignments_filtered <- alignments |>
  dplyr::group_by(read.id) |>
  dplyr::filter(dplyr::n_distinct(hla.gene) == 1) |>
  dplyr::ungroup()

# Read the TSV summarizing reads (read ID, cell barcode, UMI, reference, position)
reads <- utils::read.table(file = read_file, header = FALSE, sep = "\t")
base::colnames(reads) <- base::c("read.id", "cell.barcode", "umi", "chr", "pos")

# Merge filtered alignments with read metadata by read ID
alignments_filtered <- base::merge(x = alignments_filtered, y = reads, by = "read.id", all.x = TRUE)

# Remove reads data frame to free memory
base::rm(reads)

# Import the Gencode GTF annotation file as a data frame
gtf <- base::as.data.frame(rtracklayer::import(con = gtf_file, format = "gtf"))

# Subset the GTF to gene-level annotations on chromosome 6 only, retaining genomic coordinates and gene names
gtf <- gtf[gtf$type == "gene" & gtf$seqnames == "chr6", base::c("start", "end", "gene_name")]

# Annotate each read with the gene it overlaps based on its genomic position
alignments_filtered$gene <- base::sapply(X = alignments_filtered$pos, FUN = function(pos){
  # Unmapped reads have position 0 and are not assigned a gene
  if(pos == 0){return("")}
  # For mapped reads, identify genes whose genomic interval overlaps the read position
  if(pos != 0){
    # Identify all genes whose start–end range contains the read position
    genes <- gtf[gtf$start <= pos & gtf$end >= pos, "gene_name"]
    # Assign gene only if there is exactly one overlapping gene
    if(base::length(x = genes) == 1){return(genes)}
    # If zero or multiple genes overlap, treat as ambiguous and leave unassigned
    if(base::length(x = genes) != 1){return("")}
  }
})

# Convert placeholder values to NA for consistency
alignments_filtered$chr[alignments_filtered$chr == "*"] <- NA
alignments_filtered$pos[alignments_filtered$pos == "0"] <- NA
alignments_filtered$gene[alignments_filtered$gene == ""] <- NA

# Retain reads that either have no gene assignment or overlap an HLA gene
alignments_filtered <- alignments_filtered[base::is.na(alignments_filtered$gene) | base::grepl(pattern = "HLA-", x = alignments_filtered$gene), ]

# Remove ambiguous UMIs by retaining only cell barcode–UMI combinations that map to a single HLA gene
alignments_filtered <- alignments_filtered |>
  dplyr::group_by(cell.barcode, umi) |>
  dplyr::filter(dplyr::n_distinct(hla.gene) == 1) |>
  dplyr::ungroup()

#
alignments_filtered <- alignments_filtered[, base::c("read.id", "cell.barcode", "umi", "chr", "pos", "gene", "allele", "hla.gene", "start", "end", "n.mismatches", "n.gaps", "match")]

# Save the filtered alignments to a TSV file
write.table(x = alignments_filtered, file = base::sub(pattern = "\\.tsv$", replacement = "_filtered.tsv", alignment_file), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

#
alignments_filtered <- alignments_filtered[, base::c("cell.barcode", "umi", "hla.gene", "allele")] |>
  dplyr::group_by(cell.barcode, umi, hla.gene) |>
  dplyr::mutate(allele = base::ifelse(dplyr::n_distinct(allele) > 1, base::sub(pattern = "\\*.*$", replacement = "", x = allele), allele)) |>
  dplyr::ungroup()

# For each allele-barcode pair:
# - Count the number of unique UMIs observed
# - Pivot the data to a wide format: alleles as rows, barcodes as columns
# - Fill missing values with 0 (no UMIs observed)
mat <- alignments_filtered |>
  dplyr::group_by(allele, cell.barcode) |>
  dplyr::summarise(n_umi = dplyr::n_distinct(umi), .groups = "drop") |>
  tidyr::pivot_wider(names_from = cell.barcode, values_from = n_umi, values_fill = 0) |>
  base::as.data.frame()

# Set rownames to allele names and convert the remaining data frame to a numeric matrix
base::rownames(mat) <- mat$allele
mat <- base::as.matrix(mat[,-1])

# Convert to sparse matrix to save space for mostly-zero matrices
sparse_mat <- Matrix::Matrix(data = mat, sparse = TRUE)

# Create output directory for allele count matrix
base::dir.create(path = base::paste0(output_dir, "/allele_count_matrix/"))
# Save the matrix in Matrix Market format (widely compatible)
Matrix::writeMM(obj = sparse_mat, file = glue::glue("{output_dir}/allele_count_matrix/counts.mtx"))
# Save allele and barcode identifiers to separate text files
utils::write.table(x = base::rownames(mat), file = glue::glue("{output_dir}/allele_count_matrix/alleles.txt"), quote = FALSE, col.names = FALSE, row.names = FALSE)
utils::write.table(x = base::colnames(mat), file = glue::glue("{output_dir}/allele_count_matrix/barcodes.txt"), quote = FALSE, col.names = FALSE, row.names = FALSE)
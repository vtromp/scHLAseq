#!/usr/bin/env Rscript

# Define the command-line options for the script using optparse
option_list <- list(
  optparse::make_option(opt_str = base::c("-r", "--reads"), action = "store", type = "character", help = "TSV file summarizing reads with cell barcode and UMI information", metavar = "TSV"),
  optparse::make_option(opt_str = base::c("-a", "--alignment-file"), action = "store", type = "character", help = "TSV file produced by vsearch containing read-to-HLA alignments", metavar = "TSV"),
  optparse::make_option(opt_str = base::c("-o", "--output-dir"), action = "store", type = "character", default = base::getwd(), help = "Directory for output (default: working directory)", metavar = "DIR"),
)

# Parse command-line arguments based on the defined options
opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

# Assign parsed arguments to internal variables for use in the script
reads_overview <- opt$`reads`
alignment_file <- opt$`alignment-file`
output_dir <- opt$`output-dir`

# Remove trailing slash if present
output_dir <- base::sub(pattern = "/$", replacement = "", x = output_dir)

# Read the vsearch alignment TSV file (no header) and assign column names
alignments <- utils::read.table(file = alignment_file, header = FALSE, sep = "\t")
base::colnames(alignments) <- base::c("read.id", "allele", "match", "n.mismatches", "n.gaps", "start", "end")

# Extract the HLA gene name (everything before the "*")
alignments$gene <- base::sub(pattern = "\\*.*$", replacement = "", alignments$allele)

# Keep only relevant columns for downstream processing
alignments <- alignments[, base::c("read.id", "gene", "allele")]

# Split alignments by read ID to process reads individually
split_alignments <- base::split(x = alignments, f = alignments$read.id)

# Apply filtering rules:
# - If a read has only 1 alignment, keep it.
# - If a read has multiple alignments:
#   - Keep only if all alignments map to the same HLA gene (collapse to gene level)
#   - Remove the read if alignments map to different genes
alignments_filtered <- base::lapply(X = split_alignments, FUN = function(alignment){
  # Single alignment: keep as is
  if(base::nrow(alignment) == 1){return(alignment)}
  # Multiple alignments: check gene consistency
  if(base::nrow(alignment) > 1){
    # All alignments map to the same gene → collapse to gene
    if(base::length(x = base::unique(x = alignment$gene)) == 1){alignment$allele <- alignment$gene; return(alignment[1, ])}
    # Alignments map to different genes → discard read
    if(base::length(x = base::unique(x = alignment$gene)) > 1){return(NULL)}
  }
})

# Combine filtered alignments back into a single data frame
alignments <- base::do.call(what = base::rbind, args = alignments_filtered)

# Read the TSV summarizing reads (read ID, cell barcode, UMI, reference, position)
reads <- utils::read.table(file = reads_overview, header = FALSE, sep = "\t")
base::colnames(reads) <- base::c("read.id", "cell.barcode", "umi", "chr", "pos")

# Merge filtered alignments with read metadata by read ID
alignments <- base::merge(x = alignments, y = reads, by = "read.id", all.x = TRUE)

# Remove reads data frame to free memory
base::rm(reads)

# For each allele-barcode pair:
# - Count the number of unique UMIs observed
# - Pivot the data to a wide format: alleles as rows, barcodes as columns
# - Fill missing values with 0 (no UMIs observed)
mat <- alignments |>
  dplyr::group_by(allele, cell.barcode) |>
  dplyr::summarise(n_umi = dplyr::n_distinct(umi), .groups = "drop") |>
  tidyr::pivot_wider(names_from = cell.barcode, values_from = n_umi, values_fill = 0) |>
  base::as.data.frame()

# Set rownames to allele names and convert the remaining data frame to a numeric matrix
base::rownames(mat) <- mat$allele
mat <- base::as.matrix(mat[,-1])

# Convert to sparse matrix to save space for mostly-zero matrices
sparse_mat <- Matrix::Matrix(data = mat, sparse = TRUE)
# Save the matrix in Matrix Market format (widely compatible)
Matrix::writeMM(obj = sparse_mat, file = glue::glue("{output_dir}/allele_count_matrix.mtx"))
# Save allele and barcode identifiers to separate text files
utils::write.table(x = base::rownames(mat), file = glue::glue("{output_dir}/alleles.txt"), quote = FALSE, col.names = FALSE, row.names = FALSE)
utils::write.table(x = base::colnames(mat), file = glue::glue("{output_dir}/barcodes.txt"), quote = FALSE, col.names = FALSE, row.names = FALSE)
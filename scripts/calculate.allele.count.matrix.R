#!/usr/bin/env Rscript

# ------------------------------------------------------------------
# This script processes read-to-HLA alignment results to assign each
# barcode–UMI pair to a specific HLA gene and allele.
#
# Read-level alignments produced by vsearch are merged with read
# metadata, scored, and aggregated at the molecule level. Molecules
# are classified as confidently assigned or ambiguous based on
# alignment score support across HLA genes and alleles.
#
# The script outputs:
# - Tables of passed and dropped HLA molecules
# - A sparse allele-by-barcode UMI count matrix saved in an
#   'allele_count_matrix/' directory, with:
#     - counts.mtx containing the matrix data
#     - alleles.txt as row names
#     - barcodes.txt as column names
# ------------------------------------------------------------------

# Define the command-line options for the script using optparse
option_list <- list(
  optparse::make_option(opt_str = base::c("-r", "--reads"), action = "store", type = "character", help = "TSV file summarizing reads with cell barcode and UMI information", metavar = "TSV"),
  optparse::make_option(opt_str = base::c("-a", "--alignments"), action = "store", type = "character", help = "TSV file produced by vsearch containing all reported read-to-HLA class II alignments", metavar = "TSV"),
  optparse::make_option(opt_str = base::c("-o", "--output-dir"), action = "store", type = "character", default = base::getwd(), help = "Directory for output (default: working directory)", metavar = "DIR")
)

# Parse command-line arguments based on the defined options
opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

# Assign parsed arguments to internal variables for use in the script
read_file <- opt$`reads`
alignment_file <- opt$`alignments`
output_dir <- opt$`output-dir`

# Remove trailing slash if present
output_dir <- base::sub(pattern = "/$", replacement = "", x = output_dir)

# Read the vsearch alignment TSV file (no header) and assign column names
alignments <- utils::read.table(file = alignment_file, header = FALSE, sep = "\t")
base::colnames(alignments) <- base::c("read.id", "hla.allele", "match", "coverage", "length", "n.matches", "n.mismatches", "n.gaps", "read.start", "read.end", "allele.start", "allele.end")
# Extract the HLA gene name (everything before the "*")
alignments$hla.gene <- base::sub(pattern = "\\*.*$", replacement = "", alignments$hla.allele)
#
alignments$score <- alignments$n.matches / alignments$length * alignments$n.matches
# Reorder and sort the alignments dataframe
alignments <- alignments[base::order(alignments$read.id, -alignments$score, alignments$hla.allele), base::c("read.id", "hla.gene", "hla.allele", "score", "match", "coverage", "length", "n.matches", "n.mismatches", "n.gaps", "read.start", "read.end", "allele.start", "allele.end")]

# Read the TSV summarizing reads
reads <- utils::read.table(file = read_file, header = TRUE, sep = "\t")

# Merge alignments with read metadata by read ID
alignments <- base::merge(x = alignments, y = reads, by = "read.id", all.x = TRUE)

# Remove reads dataframe to free memory
base::rm(reads)

# Convert placeholder values to NA for consistency
alignments$chr[alignments$chr == "*"] <- NA
alignments$pos[alignments$pos == "0"] <- NA

# Reorder columns
alignments <- alignments[, base::c("read.id", "barcode", "umi", "chr", "pos", "gene", "region", "alignment.flag", "hla.gene", "hla.allele", "score", "match", "coverage", "length", "n.matches", "n.mismatches", "n.gaps", "read.start", "read.end", "allele.start", "allele.end")]
write.table(x = alignments, file = base::sub(pattern = "\\.tsv$", replacement = "_processed.tsv", x = alignment_file), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Initialize a data frame to store confidently assigned HLA molecules
passed_hla_molecules <- base::data.frame(barcode = base::character(0), umi = base::character(), hla.gene = base::character(0))
# Initialize a data frame to store molecules with ambiguous HLA assignments
dropped_hla_molecules <- base::data.frame(barcode = base::character(0), umi = base::character(), hla.genes = base::character(0))

# Identify all unique molecule identifiers defined by barcode and UMI
all_molecules <- base::unique(x = alignments[, base::c("barcode", "umi")])

# Iterate over each unique molecule
for(row in 1:base::nrow(all_molecules)){
  # Extract the barcode and UMI for the current molecule
  molecule <- all_molecules[row, ]
  # Subset all reads corresponding to the current molecule
  reads <- alignments[alignments$barcode == molecule$barcode & alignments$umi == molecule$umi, ]
  # Determine the unique HLA genes observed among reads for this molecule
  hla_genes <- base::unique(x = reads$hla.gene)
  # Check whether the molecule maps to a single HLA gene
  if(base::length(x = hla_genes) == 1){
    # Record the molecule as confidently assigned to that HLA gene
    passed_hla_molecules <- base::rbind(passed_hla_molecules, base::data.frame(barcode = molecule$barcode, umi = molecule$umi, hla.gene = hla_genes))
    # Skip further processing for this molecule
    next
  }
  # Check whether the molecule maps to multiple HLA genes
  if(base::length(x = hla_genes) > 1){
    # Compute a per-gene score by summing the best alignment score per read
    scores <- base::sapply(X = hla_genes, FUN = function(gene){base::sum(base::sapply(X = base::unique(reads[reads$hla.gene == gene, "read.id"]), FUN = function(read){base::max(reads[reads$read.id == read & reads$hla.gene == gene, "score"])}))})
    # Retain HLA genes whose scores contribute more than 60% of the total support (if present)
    if(base::sum(scores > 0.6 * base::sum(scores)) >= 1){hla_genes <- hla_genes[scores > 0.6 * base::sum(scores)]}
    # Check whether a single top-scoring HLA gene remains
    if(base::length(x = hla_genes) == 1){
      # Store the molecule as assigned to the highest-scoring HLA gene
      passed_hla_molecules <- base::rbind(passed_hla_molecules, base::data.frame(barcode = molecule$barcode, umi = molecule$umi, hla.gene = hla_genes))
    }
    # Check whether multiple HLA genes remain tied after scoring
    if(base::length(x = hla_genes) > 1){
      # Record molecules that cannot be confidently assigned to a single HLA gene
      dropped_hla_molecules <- base::rbind(dropped_hla_molecules, base::data.frame(barcode = molecule$barcode, umi = molecule$umi, hla.genes = base::paste0(base::sort(hla_genes), collapse = "/")))
    }
  }
}

# Assign a specific HLA allele to each confidently assigned molecule
passed_hla_molecules$hla.allele <- base::sapply(X = 1:base::nrow(passed_hla_molecules), function(row){
  # Subset aligned reads supporting the assigned HLA gene for this molecule
  reads <- alignments[alignments$barcode == passed_hla_molecules$barcode[row] & alignments$umi == passed_hla_molecules$umi[row] & alignments$hla.gene == passed_hla_molecules$hla.gene[row], ]
  # Identify all unique HLA alleles supported by these reads
  hla_alleles <- base::unique(x = reads$hla.allele)
  # Return the allele directly if only one allele is observed
  if(base::length(x = hla_alleles) == 1){return(hla_alleles)}
  # Resolve cases where multiple HLA alleles are observed
  if(base::length(x = hla_alleles) > 1){
    # Compute total alignment score support for each HLA allele
    scores <- base::sapply(X = hla_alleles, FUN = function(allele){base::sum(reads[reads$hla.allele == allele, "score"])})
    # Retain HLA alleles whose scores contribute more than 60% of the total support (if present)
    if(base::sum(scores > 0.6 * base::sum(scores)) >= 1){hla_alleles <- hla_alleles[scores > 0.6 * base::sum(scores)]}
    # Return the allele if a single best-supported allele remains
    if(base::length(x = hla_alleles) == 1){return(hla_alleles)}
    # Fall back to returning the HLA gene if allele assignment remains ambiguous
    if(base::length(x = hla_alleles) > 1){return(passed_hla_molecules$hla.gene[row])}
  }
})

# Write confidently assigned and dropped HLA molecules to separate TSV files
utils::write.table(x = passed_hla_molecules, file = glue::glue("{output_dir}/hla_molecule_info_passed.tsv"), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
utils::write.table(x = dropped_hla_molecules, file = glue::glue("{output_dir}/hla_molecule_info_dropped.tsv"), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

# For each allele-barcode pair:
# - Count the number of unique UMIs observed
# - Pivot the data to a wide format: alleles as rows, barcodes as columns
# - Fill missing values with 0 (no UMIs observed)
# - Order the dataframe by allele name
mat <- passed_hla_molecules |>
  dplyr::group_by(hla.allele, barcode) |>
  dplyr::summarise(n_umi = dplyr::n_distinct(umi), .groups = "drop") |>
  tidyr::pivot_wider(names_from = barcode, values_from = n_umi, values_fill = 0) |>
  dplyr::arrange(hla.allele) |>
  base::as.data.frame()

# Set rownames to allele names and convert the remaining data frame to a numeric matrix
base::rownames(mat) <- mat$hla.allele
mat <- base::as.matrix(mat[,-1])

# Convert to sparse matrix to save space for mostly-zero matrices
sparse_mat <- Matrix::Matrix(data = mat, sparse = TRUE)

# Create output directory for allele count matrix
base::dir.create(path = base::paste0(output_dir, "/allele_count_matrix/"))
# Save the matrix in Matrix Market format (widely compatible)
base::invisible(Matrix::writeMM(obj = sparse_mat, file = glue::glue("{output_dir}/allele_count_matrix/counts.mtx")))
# Save allele and barcode identifiers to separate text files
utils::write.table(x = base::rownames(mat), file = glue::glue("{output_dir}/allele_count_matrix/alleles.txt"), quote = FALSE, col.names = FALSE, row.names = FALSE)
utils::write.table(x = base::colnames(mat), file = glue::glue("{output_dir}/allele_count_matrix/barcodes.txt"), quote = FALSE, col.names = FALSE, row.names = FALSE)

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
base::colnames(alignments) <- base::c("read.id", "hla.allele", "match", "n.mismatches", "n.gaps", "start", "end")
# Extract the HLA gene name (everything before the "*")
alignments$hla.gene <- base::sub(pattern = "\\*.*$", replacement = "", alignments$hla.allele)
# Extract transcript number from the allele string and remove the transcript suffix from the allele string
alignments$transcript <- base::sub(pattern = "^HLA-.*_transcript([0-9]+)$", replacement = "\\1", x = alignments$hla.allele)
alignments$hla.allele <- base::sub(pattern = "_transcript[0-9]+$", replacement = "", x = alignments$hla.allele)
# Reorder and sort the alignments dataframe
alignments <- alignments[base::order(alignments$read.id, alignments$hla.allele, alignments$transcript), base::c("read.id", "hla.gene", "hla.allele", "transcript", "match", "n.mismatches", "n.gaps", "start", "end")]

# For each read and HLA allele pair, keep only the alignment corresponding to the transcript with the highest match score
alignments <- alignments |>
  dplyr::group_by(read.id, hla.allele) |>
  dplyr::slice_max(order_by = match, n = 1, with_ties = FALSE) |>
  dplyr::ungroup()

# Read the TSV summarizing reads (read ID, cell barcode, UMI, reference sequence / chromosome, position)
reads <- utils::read.table(file = read_file, header = FALSE, sep = "\t")
base::colnames(reads) <- base::c("read.id", "barcode", "umi", "chr", "pos")

# Merge alignments with read metadata by read ID
alignments <- base::merge(x = alignments, y = reads, by = "read.id", all.x = TRUE)

# Remove reads dataframe to free memory
base::rm(reads)

# Import the Gencode GTF annotation file as a data frame
gtf <- base::as.data.frame(rtracklayer::import(con = gtf_file, format = "gtf"))

# Subset the GTF to gene-level annotations on chromosome 6 only, retaining genomic coordinates and gene names
gtf <- gtf[gtf$type == "gene" & gtf$seqnames == "chr6", base::c("start", "end", "gene_name")]

# Annotate each read with the gene it overlaps based on its genomic position
alignments$gene <- base::sapply(X = alignments$pos, FUN = function(pos){
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

# Remove gtf dataframe to free memory
base::rm(gtf)

# Convert placeholder values to NA for consistency
alignments$chr[alignments$chr == "*"] <- NA
alignments$pos[alignments$pos == "0"] <- NA
alignments$gene[alignments$gene == ""] <- NA

# Reorder columns and save the alignments to a TSV file
alignments <- alignments[base::order(alignments$pos, alignments$read.id), base::c("read.id", "barcode", "umi", "chr", "pos", "gene", "hla.gene", "hla.allele", "transcript", "match", "n.mismatches", "n.gaps", "start", "end")]
write.table(x = alignments, file = base::sub(pattern = "\\.tsv$", replacement = "_processed.tsv", x = alignment_file), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Subset reads that were originally successfully mapped but NOT to a HLA gene and write them to a separate TSV file
alignments_overlapping <- alignments[!base::is.na(alignments$gene) & !base::grepl(pattern = "HLA-", x = alignments$gene), ]
write.table(x = alignments_overlapping, file = base::sub(pattern = "\\.tsv$", replacement = "_overlapping_with_non_HLA_genes.tsv", x = alignment_file), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
# Remove all alignments belong to reads that originally mapped to a non-HLA gene
alignments <- alignments[!alignments$read.id %in% alignments_overlapping$read.id, ]

# Subset reads that map to more than one distinct HLA gene and write them to a separate TSV file
alignments_double_mapped <- alignments |>
  dplyr::group_by(read.id) |>
  dplyr::filter(dplyr::n_distinct(hla.gene) > 1) |>
  dplyr::ungroup()
write.table(x = alignments_double_mapped, file = base::sub(pattern = "\\.tsv$", replacement = "_double_mapped.tsv", x = alignment_file), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
# Remove all alignments belonging to reads that map to multiple HLA genes
alignments <- alignments[!alignments$read.id %in% alignments_double_mapped$read.id, ]

# Collapse UMIs mapping to multiple alleles by removing the allele suffix after "*"
alignments <- alignments[, base::c("barcode", "umi", "hla.gene", "hla.allele")] |>
  dplyr::group_by(barcode, umi, hla.gene) |>
  dplyr::mutate(hla.allele = base::ifelse(dplyr::n_distinct(hla.allele) > 1, base::sub(pattern = "\\*.*$", replacement = "", x = hla.allele), hla.allele)) |>
  dplyr::ungroup()
# Remove any duplicate rows created by collapsing allele annotations
alignments <- alignments[!base::duplicated(x = alignments), ]
# Write the final molecule-level HLA annotation table to a TSV file
write.table(x = alignments, file = glue::glue("{output_dir}/molecule_info_hla.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# For each allele-barcode pair:
# - Count the number of unique UMIs observed
# - Pivot the data to a wide format: alleles as rows, barcodes as columns
# - Fill missing values with 0 (no UMIs observed)
# - Order the dataframe by allele name
mat <- alignments |>
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
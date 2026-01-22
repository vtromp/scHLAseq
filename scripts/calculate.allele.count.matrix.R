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
base::colnames(alignments) <- base::c("read.id", "hla.allele", "match", "coverage", "n.mismatches", "n.gaps", "read.start", "read.end", "allele.start", "allele.end")
# Extract the HLA gene name (everything before the "*")
alignments$hla.gene <- base::sub(pattern = "\\*.*$", replacement = "", alignments$hla.allele)
# Extract transcript number from the allele string and remove the transcript suffix from the allele string
alignments$transcript <- base::sub(pattern = "^HLA-.*_transcript([0-9]+)$", replacement = "\\1", x = alignments$hla.allele)
alignments$hla.allele <- base::sub(pattern = "_transcript[0-9]+$", replacement = "", x = alignments$hla.allele)
# Reorder and sort the alignments dataframe
alignments <- alignments[base::order(alignments$read.id, alignments$hla.allele, alignments$transcript), base::c("read.id", "hla.gene", "hla.allele", "transcript", "match", "coverage", "n.mismatches", "n.gaps", "read.start", "read.end", "allele.start", "allele.end")]

# For each read and HLA allele pair, keep only the alignment corresponding to the transcript with the highest coverage and match score
# In case of ties, the match with lowest transcript number is selected (with transcript 1, containing all exons, is selected by default)
alignments_processed <- alignments |>
  dplyr::group_by(read.id, hla.allele) |>
  dplyr::arrange(dplyr::desc(x = coverage), dplyr::desc(x = match)) |>
  dplyr::slice(1) |>
  dplyr::ungroup()

# Read the TSV summarizing reads (read ID, cell barcode, UMI, reference sequence / chromosome, position)
reads <- utils::read.table(file = read_file, header = FALSE, sep = "\t")
base::colnames(reads) <- base::c("read.id", "barcode", "umi", "chr", "pos", "type")

# Merge alignments with read metadata by read ID
alignments_processed <- base::merge(x = alignments_processed, y = reads, by = "read.id", all.x = TRUE)

# Convert placeholder values to NA for consistency
alignments_processed$chr[alignments_processed$chr == "*"] <- NA
alignments_processed$pos[alignments_processed$pos == "0"] <- NA

# Remove reads dataframe to free memory
base::rm(reads)

# Import the Gencode GTF annotation file as a data frame
gtf <- base::as.data.frame(rtracklayer::import(con = gtf_file, format = "gtf"))

# Subset the GTF to gene-level annotations on chromosome 6 only, retaining genomic coordinates and gene names
gtf <- gtf[gtf$type == "gene" & gtf$seqnames == "chr6", base::c("start", "end", "gene_name")]

# Annotate each read with the gene it overlaps based on its genomic position
alignments_processed$gene <- base::sapply(X = alignments_processed$pos, FUN = function(pos){
  # Unmapped reads have position 0 and are not assigned a gene
  if(base::is.na(pos)){return(NA)}
  # For mapped reads, identify genes whose genomic interval overlaps the read position
  if(!base::is.na(pos)){
    # Identify all genes whose start–end range contains the read position
    genes <- gtf[gtf$start <= pos & gtf$end >= pos, "gene_name"]
    # Assign gene only if there is exactly one overlapping gene
    if(base::length(x = genes) == 1){return(genes)}
    # If zero or multiple genes overlap, treat as ambiguous and leave unassigned
    if(base::length(x = genes) != 1){return(NA)}
  }
})

# Remove gtf dataframe to free memory
base::rm(gtf)

# Reorder rows and column and write to a separate (processed) alignment TSV file
alignments_processed <- alignments_processed[base::order(alignments_processed$pos, alignments_processed$read.id), base::c("read.id", "barcode", "umi", "chr", "pos", "type", "gene", "hla.gene", "hla.allele", "transcript", "match", "coverage", "n.mismatches", "n.gaps", "read.start", "read.end", "allele.start", "allele.end")]
write.table(x = alignments_processed, file = base::sub(pattern = "\\.tsv$", replacement = "_processed.tsv", x = alignment_file), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Subset reads that were originally successfully mapped but NOT to a HLA gene and write them to a separate TSV file
alignments_overlapping <- alignments_processed[!base::is.na(alignments_processed$gene) & !base::grepl(pattern = "HLA-", x = alignments_processed$gene), ]
write.table(x = alignments_overlapping, file = base::sub(pattern = "\\.tsv$", replacement = "_overlapping_with_non_HLA_genes.tsv", x = alignment_file), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
# Remove all alignments belong to reads that originally mapped to a non-HLA gene
alignments_filtered_1 <- alignments_processed[!alignments_processed$read.id %in% alignments_overlapping$read.id, ]

# For each read, retain only alignments from the HLA gene whose alignments consistently show the highest coverage and match (using the minimum value per gene as a conservative comparison)
# If multiple genes are tied on both metrics, all tied genes and their alignments are kept
alignments_filtered_2 <- alignments_filtered_1 |>
  dplyr::group_by(read.id) |>
  dplyr::mutate(row_id = dplyr::row_number()) |>
  dplyr::ungroup() |>
  dplyr::group_by(read.id) |>
  dplyr::mutate(dominated = purrr::map_lgl(.x = row_id, .f = function(row){
    this <- dplyr::cur_data_all()[row, ]
    others <- dplyr::cur_data_all() |> dplyr::filter(hla.gene != this$hla.gene)
    base::any(others$coverage >= this$coverage & others$match >= this$match & (others$coverage > this$coverage | others$match > this$match))
  })) |>
  dplyr::group_by(read.id, hla.gene) |>
  dplyr::filter(base::any(!dominated)) |>
  dplyr::select(-row_id, -dominated) |>
  dplyr::ungroup()

# Subset reads that map to more than one distinct HLA gene and write them to a separate TSV file
alignments_double_mapped <- alignments_filtered_2 |>
  dplyr::group_by(read.id) |>
  dplyr::filter(dplyr::n_distinct(hla.gene) > 1) |>
  dplyr::ungroup()
write.table(x = alignments_double_mapped, file = base::sub(pattern = "\\.tsv$", replacement = "_double_mapped.tsv", x = alignment_file), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
# Remove all alignments belonging to reads that map to multiple HLA genes
alignments_filtered_3 <- alignments_filtered_2[!alignments_filtered_2$read.id %in% alignments_double_mapped$read.id, ]

# For each barcode+umi combination, retain only alignments from the HLA gene whose alignments consistently show the highest coverage and match (using the minimum value per gene as a conservative comparison)
# If multiple genes are tied on both metrics, all tied genes and their alignments are kept
alignments_filtered_4 <- alignments_filtered_3 |>
  dplyr::group_by(barcode, umi) |>
  dplyr::mutate(row_id = dplyr::row_number()) |>
  dplyr::ungroup() |>
  dplyr::group_by(barcode, umi) |>
  dplyr::mutate(dominated = purrr::map_lgl(.x = row_id, .f = function(row){
    this <- dplyr::cur_data_all()[row, ]
    others <- dplyr::cur_data_all() |> dplyr::filter(hla.gene != this$hla.gene)
    base::any(others$coverage >= this$coverage & others$match >= this$match & (others$coverage > this$coverage | others$match > this$match))
  })) |>
  dplyr::group_by(barcode, umi, hla.gene) |>
  dplyr::filter(base::any(!dominated)) |>
  dplyr::select(-row_id, -dominated) |>
  dplyr::ungroup()

# Identify ambiguous molecules in the alignment dataframe (more than 1 HLA gene per barcode+UMI combination) and write them to a separate TSV file
alignments_ambiguous_molecules <- alignments_filtered_4 |>
  dplyr::group_by(barcode, umi) |>
  dplyr::filter(dplyr::n_distinct(hla.gene) > 1) |>
  dplyr::ungroup()
write.table(x = alignments_ambiguous_molecules, file = base::sub(pattern = "\\.tsv$", replacement = "_ambiguous_molecules.tsv", x = alignment_file), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
# Remove alignments from ambiguous molecules
alignments_filtered_5 <- alignments_filtered_4 |>
  dplyr::anti_join(y = alignments_ambiguous_molecules, by = base::c("barcode", "umi"))

# For each read, retain only alignments from the HLA allele whose alignments consistently show the highest coverage and match (using the minimum value per allele as a conservative comparison)
# If multiple alleles are tied on both metrics, all tied alleles and their alignments are kept.
alignments_filtered_6 <- alignments_filtered_5 |>
  dplyr::group_by(read.id) |>
  dplyr::mutate(row_id = dplyr::row_number()) |>
  dplyr::ungroup() |>
  dplyr::group_by(read.id) |>
  dplyr::mutate(dominated = purrr::map_lgl(row_id, function(i){
    this <- dplyr::cur_data_all()[i, ]
    others <- dplyr::cur_data_all() |> dplyr::filter(hla.allele != this$hla.allele)
    base::any(others$coverage >= this$coverage & others$match >= this$match & (others$coverage > this$coverage | others$match > this$match))
  })) |>
  dplyr::group_by(read.id, hla.allele) |>
  dplyr::filter(base::any(!dominated)) |>
  dplyr::select(-row_id, -dominated) |>
  dplyr::ungroup()

# Collapse UMIs mapping to multiple alleles by removing the allele suffix after "*"
molecule_info <- alignments_filtered_6[, base::c("barcode", "umi", "hla.gene", "hla.allele")] |>
  dplyr::group_by(barcode, umi, hla.gene) |>
  dplyr::mutate(hla.allele = base::ifelse(dplyr::n_distinct(hla.allele) > 1, base::sub(pattern = "\\*.*$", replacement = "", x = hla.allele), hla.allele)) |>
  dplyr::ungroup()
# Remove any duplicate rows created by collapsing allele annotations
molecule_info <- molecule_info[!base::duplicated(x = molecule_info), ]
# Write the final molecule-level HLA annotation table to a TSV file
write.table(x = molecule_info, file = glue::glue("{output_dir}/molecule_info_hla.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# For each allele-barcode pair:
# - Count the number of unique UMIs observed
# - Pivot the data to a wide format: alleles as rows, barcodes as columns
# - Fill missing values with 0 (no UMIs observed)
# - Order the dataframe by allele name
mat <- molecule_info |>
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
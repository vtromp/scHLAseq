#!/usr/bin/env Rscript

# Define the command-line options for the script using optparse
option_list <- list(
  optparse::make_option(opt_str = base::c("-c", "--cellranger-out"), action = "store", type = "character", help = "", metavar = "DIR"),
  optparse::make_option(opt_str = base::c("-m", "--allele-count-matrix"), action = "store", type = "character", help = "", metavar = "DIR"),
  optparse::make_option(opt_str = base::c("-p", "--project-id"), action = "store", type = "character", help = "", metavar = ""),
  optparse::make_option(opt_str = base::c("-o", "--output-dir"), action = "store", type = "character", default = base::getwd(), help = "Directory for output (default: working directory)", metavar = "DIR")
)

# Parse command-line arguments based on the defined options
opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

# Assign parsed arguments to internal variables for use in the script
cellranger_output_dir <- opt$`cellranger-out`
allele_count_matrix_dir <- opt$`allele-count-matrix`
project <- opt$`project-id`
output_dir <- opt$`output-dir`

# Remove trailing slash if present
cellranger_output_dir <- base::sub(pattern = "/$", replacement = "", x = cellranger_output_dir)
allele_count_matrix_dir <- base::sub(pattern = "/$", replacement = "", x = allele_count_matrix_dir)
output_dir <- base::sub(pattern = "/$", replacement = "", x = output_dir)

# Read filtered gene expression matrix from Cell Ranger output HDF5 file and create Seurat object using RNA assay and specified project name
data <- Seurat::Read10X_h5(filename = glue::glue("{cellranger_output_dir}/outs/filtered_feature_bc_matrix.h5"), use.names = TRUE, unique.features = TRUE)
data <- Seurat::CreateSeuratObject(counts = data, assay = "RNA", project = project)

# Extract default HLA gene counts from Seurat object (only HLA class II genes) and convert into tidy long format (gene, barcode, count)
default_counts <- base::as.data.frame(Seurat::GetAssayData(object = data, assay = "RNA", layer = "counts")[base::grep(pattern = "^HLA-(DP|DQ|DR)(A|B)[0-9]?$", base::rownames(data), value = TRUE), ]) |>
  tibble::rownames_to_column(var = "gene") |>
  tidyr::pivot_longer(cols = -gene, names_to = "barcode", values_to = "count")

# Read allele count matrix in Matrix Market format and set row and column names
allele_counts <- base::as.data.frame(Matrix::readMM(file = glue::glue("{allele_count_matrix_dir}/counts.mtx")))
base::rownames(allele_counts) <- base::readLines(con = glue::glue("{allele_count_matrix_dir}/alleles.txt"))
base::colnames(allele_counts) <- base::readLines(con = glue::glue("{allele_count_matrix_dir}/barcodes.txt"))

# Copy allele counts dataframe for collapsing by gene
allele_counts_collapsed <- allele_counts
# Loop through each unique gene name (before allele suffix)
for(gene in base::unique(base::sub(pattern = "\\*.*$", replacement = "", x = base::rownames(allele_counts)))){
  # Select all rows corresponding to the current gene (all alleles)
  gene_rows <- allele_counts_collapsed[base::grepl(pattern = gene, x = base::rownames(allele_counts_collapsed)), , drop = FALSE]
  # Sum counts across all allele rows for each barcode
  gene_sum <- base::colSums(x = gene_rows)
  # Remove all allele-specific rows for this gene
  allele_counts_collapsed <- allele_counts_collapsed[!base::grepl(pattern = gene, x = base::rownames(allele_counts_collapsed)), , drop = FALSE]
  # Add a new row with the gene name containing the summed counts
  allele_counts_collapsed[gene, ] <- gene_sum
}

# Convert collapsed allele counts into tidy long format (gene, barcode, count)
personalized_counts <- allele_counts_collapsed |>
  tibble::rownames_to_column(var = "gene") |>
  tidyr::pivot_longer(cols = -gene, names_to = "barcode", values_to = "count")

# Create a new dataframe with all possible barcode-gene combinations
hla_expression_df <- base::expand.grid(barcode = base::unique(x = default_counts$barcode), gene = base::unique(x = base::c(default_counts$gene, personalized_counts$gene)), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
# Merge default counts
hla_expression_df <- base::merge(x = hla_expression_df, y = default_counts, by = base::c("barcode", "gene"), all.x = TRUE)
base::colnames(hla_expression_df)[3] <- "default.count"
# Merge personalized counts
hla_expression_df <- base::merge(x = hla_expression_df, y = personalized_counts, by = base::c("barcode", "gene"), all.x = TRUE)
base::colnames(hla_expression_df)[4] <- "personalized.count"
# Replace all NA values with 0
hla_expression_df[base::is.na(hla_expression_df)] <- 0








#
allele_counts <- allele_counts |>
  tibble::rownames_to_column(var = "allele") |>
  tidyr::pivot_longer(cols = -allele, names_to = "barcode", values_to = "allelic.expression") |>
  dplyr::filter(barcode %in% hla_expression_df$barcode) |>
  dplyr::mutate(gene = base::sub(pattern = "\\*.*$", replacement = "", x = allele))

#
hla_expression_df <- base::merge(x = hla_expression_df, y = allele_counts, by = base::c("barcode", "gene"), all.x = TRUE)
#
hla_expression_df <- hla_expression_df[base::order(hla_expression_df$barcode, hla_expression_df$gene) , base::c("barcode", "gene", "allele", "default.expression", "personalized.expression", "allelic.expression")]

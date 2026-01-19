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

#
data <- Seurat::Read10X_h5(filename = glue::glue("{cellranger_output_dir}/outs/filtered_feature_bc_matrix.h5"), use.names = TRUE, unique.features = TRUE)
data <- Seurat::CreateSeuratObject(counts = data, assay = "RNA", project = project)

#
default_counts <- base::as.data.frame(Seurat::GetAssayData(object = data, assay = "RNA", layer = "counts")[base::grep(pattern = "^HLA-(DP|DQ|DR)(A|B)[0-9]?$", base::rownames(data), value = TRUE), ]) |>
  tibble::rownames_to_column(var = "gene") |>
  tidyr::pivot_longer(cols = -gene, names_to = "cell.barcode", values_to = "default.expression")

#
allele_counts <- base::as.data.frame(Matrix::readMM(file = glue::glue("{allele_count_matrix_dir}/counts.mtx")))
base::rownames(allele_counts) <- base::readLines(con = glue::glue("{allele_count_matrix_dir}/alleles.txt"))
base::colnames(allele_counts) <- base::readLines(con = glue::glue("{allele_count_matrix_dir}/barcodes.txt"))

#
personalized_counts <- stats::aggregate(x = allele_counts, by = base::list(gene = base::sub(pattern = "\\*.*$", replacement = "", base::rownames(allele_counts))), FUN = base::sum) |>
  tidyr::pivot_longer(cols = -gene, names_to = "cell.barcode", values_to = "personalized.expression")

#
hla_expression_df <- base::merge(x = default_counts, y = personalized_counts, by = base::c("cell.barcode", "gene"), all = TRUE)
#
hla_expression_df$default.expression[base::is.na(hla_expression_df$default.expression)] <- 0
hla_expression_df$personalized.expression[base::is.na(hla_expression_df$personalized.expression)] <- 0
#
hla_expression_df <- hla_expression_df[hla_expression_df$cell.barcode %in% default_counts$cell.barcode, ]

#
allele_counts <- allele_counts |>
  tibble::rownames_to_column(var = "allele") |>
  tidyr::pivot_longer(cols = -allele, names_to = "cell.barcode", values_to = "allelic.expression") |>
  dplyr::filter(cell.barcode %in% hla_expression_df$cell.barcode) |>
  dplyr::mutate(gene = base::sub(pattern = "\\*.*$", replacement = "", x = allele))

#
hla_expression_df <- base::merge(x = hla_expression_df, y = allele_counts, by = base::c("cell.barcode", "gene"), all.x = TRUE)
#
hla_expression_df <- hla_expression_df[base::order(hla_expression_df$cell.barcode, hla_expression_df$gene) , base::c("cell.barcode", "gene", "allele", "default.expression", "personalized.expression", "allelic.expression")]

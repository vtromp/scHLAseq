#!/usr/bin/env Rscript

# ------------------------------------------------------------------
# This script integrates recalculated HLA allele counts with Cell 
# Ranger count output to create an updated Seurat object.
#
# Default gene-level expression values of HLA class II genes 
# derived are removed from the expression matrix and replaced with 
# gene-level expression values computed from a provided allele 
# count matrix.
#
# In addition to updating the expression matrix, both the original
# Cell Ranger–derived gene-level expression values and the
# personalized scHLAseq-derived gene-level expression values are
# log-normalized and appended to cell-level metadata. Allele-level
# expression values are also log-normalized and stored in metadata,
# enabling direct subsequent analysis of allelic expression levels 
# in a single-cell manner.
#
# The script outputs a single Seurat object with personalized HLA 
# class II genes in the expression matrix and allele-resolved HLA 
# expression stored in cell-level metadata
# ------------------------------------------------------------------

# Define the command-line options for the script using optparse
option_list <- list(
  optparse::make_option(opt_str = base::c("-c", "--cellranger-out"), action = "store", type = "character", help = "Path to Cell Ranger output directory containing the outs/filtered_feature_bc_matrix.h5 file", metavar = "DIR"),
  optparse::make_option(opt_str = base::c("-m", "--allele-count-matrix"), action = "store", type = "character", help = "Path to directory containing scHLAseq allele-by-barcode count matrix (counts.mtx, alleles.txt, barcodes.txt)", metavar = "DIR"),
  optparse::make_option(opt_str = base::c("-p", "--project-id"), action = "store", type = "character", help = "Project identifier used to label the output Seurat object", metavar = "PROJECT_ID"),
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
seurat_object <- Seurat::CreateSeuratObject(counts = data, assay = "RNA", project = project)

# Extract default HLA gene counts from Seurat object (only HLA class II genes) and convert into tidy long format (gene, barcode, count)
default_counts <- base::as.data.frame(Seurat::GetAssayData(object = seurat_object, assay = "RNA", layer = "counts")[base::grep(pattern = "^HLA-(DP|DQ|DR)(A|B)[0-9]?$", base::rownames(seurat_object), value = TRUE), ]) |>
  tibble::rownames_to_column(var = "gene") |>
  tidyr::pivot_longer(cols = -gene, names_to = "barcode", values_to = "count")

# Read allele count matrix in Matrix Market format and set row and column names
allele_counts <- base::as.data.frame(Matrix::readMM(file = glue::glue("{allele_count_matrix_dir}/counts.mtx")))
base::rownames(allele_counts) <- base::readLines(con = glue::glue("{allele_count_matrix_dir}/alleles.txt"))
base::colnames(allele_counts) <- base::readLines(con = glue::glue("{allele_count_matrix_dir}/barcodes.txt"))

# Add missing Seurat barcodes to allele counts with zero values
allele_counts[base::colnames(seurat_object)[!base::colnames(seurat_object) %in% colnames(allele_counts)]] <- 0

# Create a copy of allele counts for collapsing by gene
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
scHLAseq_counts <- allele_counts_collapsed |>
  tibble::rownames_to_column(var = "gene") |>
  tidyr::pivot_longer(cols = -gene, names_to = "barcode", values_to = "count")

# Create a new dataframe with all possible barcode-gene combinations
hla_expression_df <- base::expand.grid(barcode = base::unique(x = default_counts$barcode), gene = base::unique(x = base::c(default_counts$gene, scHLAseq_counts$gene)), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
# Merge default counts
hla_expression_df <- base::merge(x = hla_expression_df, y = default_counts, by = base::c("barcode", "gene"), all.x = TRUE)
base::colnames(hla_expression_df)[3] <- "default.count"
# Merge personalized counts
hla_expression_df <- base::merge(x = hla_expression_df, y = scHLAseq_counts, by = base::c("barcode", "gene"), all.x = TRUE)
base::colnames(hla_expression_df)[4] <- "scHLAseq.count"
# Replace all NA values with 0
hla_expression_df[base::is.na(hla_expression_df)] <- 0

# Extract RNA count matrix from original Seurat object
counts <- base::as.matrix(Seurat::GetAssayData(object = seurat_object, assay = "RNA", layer = "counts"))
# Remove HLA class II genes from the count matrix
counts <- counts[!base::rownames(counts) %in% base::grep(pattern = "^HLA-(DP|DQ|DR)(A|B)[0-9]?$", base::rownames(counts), value = TRUE), ]
# Append collapsed HLA class II counts to the matrix
counts <- base::rbind(counts, allele_counts_collapsed)
# Convert combined matrix to sparse format
counts <- Matrix::Matrix(data = base::as.matrix(counts), sparse = TRUE)
# Create a new Seurat object with updated counts
seurat_object_updated <- Seurat::CreateSeuratObject(counts = counts, assay = "RNA", project = project)

# Iterate over each HLA gene
for(gene in base::unique(hla_expression_df$gene)){
  # Add default and recalculated expression values for the current gene to the metadata
  seurat_object_updated@meta.data[[base::paste0(gene, ".default")]] <- base::sapply(X = base::rownames(seurat_object_updated@meta.data), FUN = function(barcode){hla_expression_df[hla_expression_df$barcode == barcode & hla_expression_df$gene == gene, "default.count"]})
  seurat_object_updated@meta.data[[base::paste0(gene, ".scHLAseq")]] <- base::sapply(X = base::rownames(seurat_object_updated@meta.data), FUN = function(barcode){hla_expression_df[hla_expression_df$barcode == barcode & hla_expression_df$gene == gene, "scHLAseq.count"]})
  # Check if gene has allele-level counts
  if(gene %in% base::rownames(allele_counts)){
    # Iterate over each allele for the gene
    for(allele in base::grep(pattern = base::paste0(gene, "\\*"), x = base::rownames(allele_counts), value = TRUE)){
      # Add allele-level expression values to the metadata
      seurat_object_updated@meta.data[[base::paste0(allele, ".scHLAseq")]] <- base::sapply(X = base::rownames(seurat_object_updated@meta.data), FUN = function(barcode){allele_counts[allele, barcode]})
    }
  }
}

# Normalize RNA expression data
seurat_object_updated <- Seurat::NormalizeData(object = seurat_object_updated, normalization.method = "LogNormalize", scale.factor = 10000)
# Likewise, log-normalize metadata expression values
for(col in base::colnames(x = seurat_object_updated@meta.data)[-c(1:3)]){seurat_object_updated@meta.data[[col]] <- base::log1p(seurat_object_updated@meta.data[[col]] / seurat_object_updated@meta.data$nCount_RNA * 10000)}

# Save the updated Seurat object in RDS format to the output directory
base::saveRDS(object = seurat_object_updated, file = glue::glue("{output_dir}/seurat_object_updated.RDS"))

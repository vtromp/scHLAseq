#!/usr/bin/env Rscript

# ------------------------------------------------------------------
# This script reads the Cell Ranger molecule_info.h5 output file,
# extracts per-molecule cell barcode, UMI, and gene information,
# filters molecules to barcodes that passed cell filtering and to
# HLA class II genes only, decodes 2-bit encoded UMIs, and writes
# the resulting molecule-level information to a TSV file.
# ------------------------------------------------------------------

# Define the command-line options for the script using optparse
option_list <- list(
  optparse::make_option(opt_str = base::c("-c", "--cellranger-out"), action = "store", type = "character", help = "Path to the Cell Ranger output directory containing outs/molecule_info.h5", metavar = "DIR"),
  optparse::make_option(opt_str = base::c("-o", "--output-file"), action = "store", type = "character", default = "molecule_info_hla_class_II.tsv", help = "Output TSV file with barcode, decoded UMI, and HLA class II gene assignment", metavar = "TSV")
)

# Parse command-line arguments based on the defined options
opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

# Assign parsed arguments to internal variables for use in the script
cellranger_output_dir <- opt$`cellranger-out`
output_file <- opt$`output-file`

# Remove trailing slash if present
cellranger_output_dir <- base::sub(pattern = "/$", replacement = "", x = cellranger_output_dir)

# Define the full path to the Cell Ranger 'molecule_info.h5' file
cellranger_molecule_info <- glue::glue("{cellranger_output_dir}/outs/molecule_info.h5")

# Read list of all cell barcodes
barcodes <- rhdf5::h5read(file = cellranger_molecule_info, name = "barcodes")

# Read zero-based indices into the 'barcodes' vector for all molecules
# +1 converts from zero-based indexing (HDF5) to one-based indexing (R)
barcode_idx <- rhdf5::h5read(file = cellranger_molecule_info, name = "barcode_idx")+1

# Read list of integer labels distinguishing data coming from distinct 10x Genomics GEM reactions
gems <- rhdf5::h5read(file = cellranger_molecule_info, name = "gem_group")

# Read zero-based indices of barcodes that passed cell filtering
filtered_idx <- rhdf5::h5read(file = cellranger_molecule_info, name = "barcode_info/pass_filter")[1,]+1

# Read list of 2-bit encoded UMI sequences
umi <- rhdf5::h5read(file = cellranger_molecule_info, name = "umi")

# Read list of all feature/gene names
features <- rhdf5::h5read(file = cellranger_molecule_info, name = "features/name")

# Read zero-based indices into the 'features' vector for all molecules
feature_idx <- rhdf5::h5read(file = cellranger_molecule_info, name = "feature_idx")+1

# Create a data.frame linking:
#  - cell barcode
#  - 2-bit-encoded UMI
#  - gene name
molecule_info <- data.frame(barcode = barcodes[barcode_idx], gem = gems, umi = umi, gene = features[feature_idx])
# Subset molecules to those belonging to barcodes that passed filtering
molecule_info <- molecule_info[molecule_info$barcode %in% barcodes[filtered_idx], ]
# Further subset molecules to HLA class II genes only
molecule_info <- molecule_info[base::grepl(pattern = "^HLA-(DP|DQ|DR)(A|B)[0-9]?$", x = molecule_info$gene), ]

# Combine the cell barcode and GEM group into a single barcode string
molecule_info$barcode <- base::paste(molecule_info$barcode, molecule_info$gem, sep = "-")
# Drop the 'gem' column
molecule_info <- molecule_info[, base::c("barcode", "umi", "gene")]

# Define function to decode a single 2-bit encoded UMI integer into a nucleotide string
decode_umi <- function(x, umi_length){
  # Set lookup table mapping 2-bit values to nucleotides (0 = A, 1 = C, 2 = G, 3 = T)
  bases <- base::c("A", "C", "G", "T")
  # Preallocate character vector for decoded bases
  out <- base::character(umi_length)
  # Iterate over each nucleotide position
  # The least significant bits encode the 3' end
  for(i in base::seq_len(umi_length)){
    # Extract the lowest 2 bits and map to a nucleotide
    out[i] <- bases[base::bitwAnd(a = x, b = 3) + 1]
    # Right-shift by 2 bits to move to the next nucleotide
    x <- base::bitwShiftR(a = x, n =2)
  }
  # Reverse to restore 5' → 3' orientation and collapse into a string
  base::paste(base::rev(out), collapse = "")
}

# Decode all UMIs (assumes a fixed UMI length of 12)
molecule_info$umi <- base::sapply(X = molecule_info$umi, FUN = function(umi) decode_umi(x = umi, umi_length = 12))

# Write filtered molecule information to a TSV file
utils::write.table(x = molecule_info, file = output_file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
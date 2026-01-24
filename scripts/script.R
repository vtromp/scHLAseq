#!/usr/bin/env Rscript

# Define the command-line options for the script using optparse
option_list <- list(
  optparse::make_option(opt_str = base::c("-c", "--cellranger-out"), action = "store", type = "character", help = "", metavar = "DIR"),
  optparse::make_option(opt_str = base::c("-o", "--output-dir"), action = "store", type = "character", default = base::getwd(), help = "Directory for output (default: working directory)", metavar = "DIR")
)

# Parse command-line arguments based on the defined options
opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

# Assign parsed arguments to internal variables for use in the script
cellranger_output_dir <- opt$`cellranger-out`
output_dir <- opt$`output-dir`

# Remove trailing slash if present
cellranger_output_dir <- base::sub(pattern = "/$", replacement = "", x = cellranger_output_dir)
output_dir <- base::sub(pattern = "/$", replacement = "", x = output_dir)

#
param <- Rsamtools::ScanBamParam(what = base::c("qname", "seq", "qual"), tag = base::c("CB", "UB"))

#
bam <- Rsamtools::scanBam(file = glue::glue("{cellranger_output_dir}/outs/possorted_genome_bam.bam"), param = param)[[1]]
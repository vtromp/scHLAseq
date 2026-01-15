#!/usr/bin/env Rscript

# Define the command-line options for the script using optparse
option_list <- list(
  optparse::make_option(opt_str = base::c("-g", "--genotype-file"), action = "store", type = "character", help = "arcasHLA genotype JSON", metavar = "JSON"),
  optparse::make_option(opt_str = base::c("-o", "--output-dir"), action = "store", type = "character", default = base::getwd(), help = "Directory for output (default: working directory)", metavar = "DIR"),
  optparse::make_option(opt_str = "--class-I-only", action = "store_true", type = "logical", default = FALSE, help = "Include only HLA class I genes"),
  optparse::make_option(opt_str = "--class-II-only", action = "store_true", type = "logical", default = FALSE, help = "Include only HLA class II genes")
)

# Parse command-line arguments based on the defined options
opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

# Assign parsed arguments to internal variables for use in the script
genotype_file <- opt$`genotype-file`
only.HLAclassI <- opt$`class-I-only`
only.HLAclassII <- opt$`class-II-only`

# Remove trailing slash if present
output_dir <- base::sub(pattern = "/$", replacement = "", x = output_dir)

# Ensure that class I and class II filters are not used simultaneously
if(only.HLAclassI && only.HLAclassII){stop("Please specify only one of --class-I or --class-II, not both.")}

# Create a (temporary) directory to store IPD-IMGT/HLA database files
base::dir.create(path = glue::glue("{output_dir}/IPD_IMGT_HLA"))

# Download IPD-IMGT/HLA database file and allele overview
utils::download.file(url = "https://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/hla.dat.zip", destfile = glue::glue("{output_dir}/IPD_IMGT_HLA/hla.dat.zip"), quiet = TRUE)
utils::download.file(url = "https://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/Allelelist.txt", destfile = glue::glue("{output_dir}/IPD_IMGT_HLA/allele.list.txt"), quiet = TRUE)
utils::download.file(url = "https://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/Allele_status.txt", destfile = glue::glue("{output_dir}/IPD_IMGT_HLA/allele.status.txt"), quiet = TRUE)

# Unzip IPD-IMGT/HLA database archive
utils::unzip(zipfile = glue::glue("{output_dir}/IPD_IMGT_HLA/hla.dat.zip"), exdir = glue::glue("{output_dir}/IPD_IMGT_HLA"))
# Remove the archive
base::unlink(x = glue::glue("{output_dir}/IPD_IMGT_HLA/hla.dat.zip"), force = TRUE)

# Extract HLA database version from metadata in allele list file
hla_db_version <- base::readLines(con = glue::glue("{output_dir}/IPD_IMGT_HLA/allele.list.txt"))
hla_db_version <- base::grep(pattern = "version:", x = hla_db_version, value = TRUE)
hla_db_version <- base::sub(pattern = "^# version: ", replacement = "", hla_db_version)
hla_db_version <- base::sub(pattern = " ", replacement = "_", hla_db_version)

# Load the full IPD-IMGT/HLA database as raw text
hla_db <- base::readLines(con = glue::glue("{output_dir}/IPD_IMGT_HLA/hla.dat"))

# Split the HLA database into individual allele records
starts <- base::c(1, base::which(hla_db == "//") + 1)
ends <- base::c(base::which(hla_db == "//") - 1, base::length(hla_db))
hla_db <- base::mapply(FUN = function(s, e){hla_db[s:e]}, starts, ends)

# Extract accession numbers from each record
accession_numbers <- base::unlist(base::sapply(X = hla_db, FUN = function(record){
  an <- base::grep(pattern = "^AC +", x = record, value = TRUE)
  an <- base::sub(pattern = "(^AC +)(HLA\\d{5})(;$)", replacement = "\\2", x = an)
  return(an)
}))
# Assign accession numbers as names of the records
base::names(hla_db) <- accession_numbers

# Read allele list and allele status CSV files
allele_list <- read.csv(file = glue::glue("{output_dir}/IPD_IMGT_HLA/allele.list.txt"), header = TRUE, comment.char = "#")
allele_status <- read.csv(file = glue::glue("{output_dir}/IPD_IMGT_HLA/allele.status.txt"), header = TRUE, comment.char = "#")
# Merge allele list with allele status information using allele name as key
allele_list <- base::merge(x = allele_list, y = allele_status, by = "Allele")
# Retain only alleles with complete (full-length) genomic or coding sequences
allele_list <- allele_list[allele_list$Partial == "Full", ]

# Optionally filter for HLA class I or HLA class II alleles
if(only.HLAclassI){allele_list <- allele_list[base::grepl(pattern = "^(A|B|C|E|F|G|H|J|K|L|N|P|R|S|T|U|V|W|Y)\\*", x = allele_list$Allele), ]}
if(only.HLAclassII){allele_list <- allele_list[base::grepl(pattern = "^(DP|DQ|DR)(A|B)[0-9]?\\*", x = allele_list$Allele), ]}

# Extract coding sequences for each allele based on exon and UTR annotations
allele_list$seq <- base::sapply(X = allele_list$AlleleID, FUN = function(allele_id){
  # Extract the full IPD-IMGT/HLA database record for this allele
  record <- hla_db[[allele_id]]
  # Extract all annotation entries for UTRs, exons, and introns
  feature_lines <- base::grep(pattern = "FT +UTR|FT +exon|FT +intron", x = record, value = TRUE)
  # Keep only exon and UTR annotation lines
  feature_lines <- base::grep(pattern = "exon|UTR", x = feature_lines, value = TRUE)
  # Extract numeric coordinate ranges from exon and UTR annotations
  features <- base::sub(pattern = "^FT +(exon|UTR) +([0-9\\.]+)$", replacement = "\\2", x = feature_lines)
  # Convert coordinate ranges into numeric start/stop positions
  features <- base::sapply(X = features, FUN = function(feature) base::as.numeric(base::strsplit(x = feature, split = "..", fixed = TRUE)[[1]]), simplify = FALSE)
  # Extract the sequence lines (everything after 'SQ')
  seq_lines <- record[(which(grepl(pattern = "^SQ +", x = record)) + 1):base::length(record)]
  # Concatenate sequence and remove non-ACGT characters
  genomic_seq <- base::toupper(base::paste(base::gsub(pattern = "[^ACGTacgt]", replacement = "", seq_lines), collapse = ""))
  # Extract and concatenate exon/UTR regions to form the coding sequence
  coding_seq <- base::paste0(base::sapply(X = features, FUN = function(feature){base::substr(x = genomic_seq, start = feature[1], stop = feature[2])}), collapse = "")
  # Return the coding sequence
  return(coding_seq)
})

# Read HLA genotype JSON input file produced by arcasHLA
hla_genotype <- jsonlite::fromJSON(genotype_file)

# Restrict genotype object to HLA class II loci (DP, DQ, DR)
hla_genotype <- hla_genotype[base::grepl(pattern = "^DP|^DQ|^DR", x = base::names(hla_genotype))]

# Initialize an empty named character vector to store HLA allele sequences
hla_genotype_seqs <- base::c()

# Iterate over each HLA gene in the genotype file
for(gene in base::names(hla_genotype)){
  # Iterate over all unique alleles for this HLA gene
  for(allele in base::unique(hla_genotype[[gene]])){
    # Find the first matching entry in allele_list whose allele field contains the current allele string (exact substring match)
    allele_info <- allele_list[base::grepl(pattern = allele, x = allele_list$Allele, fixed = TRUE), ][1,]
    # Store the sequence in the output vector, using the allele name as the key
    hla_genotype_seqs[allele_info$Allele] <- allele_info$seq
  }
}

# Convert extracted sequences into a DNAStringSet object
hla_genotype_seqs <- Biostrings::DNAStringSet(x = hla_genotype_seqs)
# Write all HLA allele sequences to a FASTA file
Biostrings::writeXStringSet(hla_genotype_seqs, filepath = glue::glue("{output_dir}/HLA_genotype.fasta"))

# Remove the temporary IPD-IMGT/HLA database directory
base::unlink(x = glue::glue("{output_dir}/IPD_IMGT_HLA/"), recursive = TRUE, force = TRUE)
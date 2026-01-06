#!/usr/bin/env Rscript

# Define the command-line options for the script using optparse
option_list <- list(
  optparse::make_option(opt_str = base::c("-g", "--genotype"), type = "character", default = NULL, help = "Path to HLA genotype JSON file", metavar = "FILE"),
  optparse::make_option(opt_str = base::c("-o", "--output-dir"), type = "character", default = base::getwd(), help = "Directory for output (default: working directory)", metavar = "DIR")
)

# Parse command-line arguments based on the defined options
opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

# Validate that the path to the genotype file was provided
if(base::is.null(opt$genotype)){stop("Error: You must provide a genotype file using -g or --genotype.")}

# Assign parsed arguments to internal variables for use in the script
genotype_file <- opt$`genotype`
output_dir    <- opt$`output-dir`

# Remove trailing slash if present
output_dir <- base::sub(pattern = "/$", replacement = "", x = output_dir)


####################################################################################################
# 1. Create output directories and define reference versions
####################################################################################################

# Print status message
base::message(base::paste0(base::format(x = base::Sys.time(), "%Y-%m-%d %H:%M:%S "), "Starting personalized reference build. Creating temporary directories to store Gencode GRCh38 (release 44) human reference genome and the IPD-IMGT/HLA database..."))

# Create directory to store Gencode reference files within the output directory
gencode_dir <- "Gencode"
base::dir.create(path = glue::glue("{output_dir}/{gencode_dir}"))

# Create directory to store IPD-IMGT/HLA database file within the output directory
IPD_IMGT_HLA_dir <- "IPD_IMGT_HLA"
base::dir.create(path = glue::glue("{output_dir}/{IPD_IMGT_HLA_dir}"))

# Define genome version and Gencode release number
genome <- "GRCh38"
release <- "44"


####################################################################################################
# 2. Download FASTA and GTF file of primary assembly from Genome Reference Consortium Human Build 38 (GRCH38), and download the IPD-IMGT/HLA database file
####################################################################################################

# Print status message
base::message(base::paste0(base::format(x = base::Sys.time(), "%Y-%m-%d %H:%M:%S "), "Downloading Gencode GRCh38 (release 44) genome FASTA and GTF files, and the IPD-IMGT/HLA database file..."))

# Construct URLs for downloading Gencode FASTA genome file, GTF annotation file, and IPD-IMGT/HLA database file
fasta_url <- glue::glue("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{release}/GRCh38.primary_assembly.genome.fa.gz")
gtf_url <- glue::glue("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{release}/gencode.v{release}.primary_assembly.annotation.gtf.gz")
hla_db_url <- "https://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/hla.dat.zip"
hla_alleles_url <- "https://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/Allelelist.txt"

# Extend the timeout limit to 5 minutes to allow large genome files to download successfully
base::options(timeout = 300)

# Download the Gencode FASTA genome file
utils::download.file(url = fasta_url, destfile = glue::glue("{output_dir}/{gencode_dir}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"), quiet = TRUE)
# Download the Gencode GTF annotation file
utils::download.file(url = gtf_url, destfile = glue::glue("{output_dir}/{gencode_dir}/Homo_sapiens.GRCh38.dna.primary_assembly.gtf.gz"), quiet = TRUE)
# Download the IPD-IMGT/HLA database archive and allele list
utils::download.file(url = hla_db_url, destfile = glue::glue("{output_dir}/{IPD_IMGT_HLA_dir}/hla.dat.zip"), quiet = TRUE)
utils::download.file(url = hla_alleles_url, destfile = glue::glue("{output_dir}/{IPD_IMGT_HLA_dir}/allelelist.txt"), quiet = TRUE)

# Unzip IPD-IMGT/HLA database archive
utils::unzip(zipfile = glue::glue("{output_dir}/{IPD_IMGT_HLA_dir}/hla.dat.zip"), exdir = glue::glue("{output_dir}/{IPD_IMGT_HLA_dir}"))


####################################################################################################
# 3. Import and preprocess reference genome and GTF annotation
#    Preprocessing follows the "Build Notes for Reference Packages" guidelines from 10x Genomics (https://www.10xgenomics.com/support/software/cell-ranger/downloads/cr-ref-build-steps)
####################################################################################################

# Print status message
base::message(base::paste0(base::format(x = base::Sys.time(), "%Y-%m-%d %H:%M:%S "), "Loading Gencode reference FASTA and GTF file and applying preprocessing following 10x Genomics guidelines..."))

# Import the Gencode FASTA file
fasta <- Biostrings::readDNAStringSet(filepath = glue::glue("{output_dir}/{gencode_dir}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"), format = "fasta")
# Import the Gencode GTF annotation file as a data frame
gtf <- base::as.data.frame(rtracklayer::import(glue::glue("{output_dir}/{gencode_dir}/Homo_sapiens.GRCh38.dna.primary_assembly.gtf.gz"), format = "gtf"))

# Loop through gene, transcript, and exon identifiers to extract version numbers
for(col in base::c("gene", "transcript", "exon")){
  # For each feature type, set up the ID and version column names
  id_col <- glue::glue("{col}_id")
  version_col <- glue::glue("{col}_version")
  # Extract the version number (digits after the final dot) into the version column
  gtf[[version_col]] <- base::sub(pattern = "(^.*\\.)(\\d{1,}$)", replacement = "\\2", gtf[[id_col]])
  # If the regex does not match (no version present), replace with empty string
  gtf[[version_col]][gtf[[version_col]] == gtf[[id_col]]] <- "" 
  # Remove version suffix from the ID (strip everything after the final dot)
  gtf[[id_col]] <- base::sub(pattern = "\\..*$", replacement = "", gtf[[id_col]])
}
# Reorder the GTF columns to a standardized structure
gtf <- gtf[, base::c("seqnames", "start", "end", "width", "strand", "source", "type", "score", "phase", "gene_id", "gene_version", "gene_type", "gene_name", "level", "tag",  "transcript_id", "transcript_version", "transcript_type", "transcript_name", "exon_number", "exon_id", "exon_version", "transcript_support_level", "havana_transcript", "hgnc_id", "havana_gene", "ont", "protein_id", "ccdsid", "artif_dupl")]

# Define the set of allowed gene and transcript biotypes for downstream filtering
allowed_biotypes <- base::c("protein_coding", "protein_coding_LoF", "lncRNA", "IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_LV_gene", "IG_V_gene", "IG_V_pseudogene", "IG_J_pseudogene", "IG_C_pseudogene", "TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_V_gene", "TR_V_pseudogene", "TR_J_pseudogene")
# Extract gene IDs that have at least one associated transcript with:
# - an allowed 'gene_type'
# - an allowed 'transcript_type'
# - no 'readthrough_transcript' tag 
allowed_gene_ids <- gtf |>
  dplyr::filter(type == "transcript") |>
  dplyr::filter(gene_type %in% allowed_biotypes) |>
  dplyr::filter(transcript_type %in% allowed_biotypes) |>
  dplyr::filter(!base::grepl(pattern = "readthrough_transcript", tag)) |>
  dplyr::pull(gene_id) |>
  base::unique()
# Subset the GTF annotation data frame to retain only entries belonging to the allowed genes
gtf_filtered <- gtf[gtf$gene_id %in% allowed_gene_ids, ]
# Filter chrY annotations: exclude all PAR-Y genes (including XGY2, ENSG00000290840), but retain SRY and ENSG00000286130, which lie within intronic or overlapping regions of XGY2/RPS4Y1
gtf_filtered <- gtf_filtered |>
  dplyr::filter(seqnames != "chrY" | (seqnames == "chrY" & start >= 2752083 & start < 56887903 & gene_id != "ENSG00000290840"))


####################################################################################################
# 4. Process HLA genotype JSON file and infer expected locus configuration
####################################################################################################

# Print status message
base::message(base::paste0(base::format(x = base::Sys.time(), "%Y-%m-%d %H:%M:%S "), "Reading HLA genotype JSON file and inferring expected DRB paralog and pseudogene configuration..."))

# Read HLA genotype JSON input file produced by arcasHLA
hla_genotype <- jsonlite::fromJSON(genotype_file)

# Restrict genotype object to HLA class II loci (DP, DQ, DR)
hla_genotype <- hla_genotype[base::grepl(pattern = "^DP|^DQ|^DR", x = base::names(hla_genotype))]

# Extract the DRB1 allele calls from the genotype object (key determinants of DRB haplotype structure)
DRB1_alleles <- hla_genotype[["DRB1"]]
# Extract the 2-digit serological group from the allele names
allele_groups <- base::unique(base::sapply(X = DRB1_alleles, FUN = function(allele){base::sub(pattern = "DRB1\\*(\\d{2}):.*$", replacement = "\\1", x = allele)}))

# Define expected DRB paralog and pseudogene composition for each DRB1 allele group, based on known HLA-DR haplotype architectures
group_dictionary <- base::list(`01` = base::c("DRB1", "DRB6", "DRB9"),
                               `03` = base::c("DRB1", "DRB2", "DRB3", "DRB9"),
                               `04` = base::c("DRB1", "DRB4", "DRB7", "DRB8", "DRB9"),
                               `07` = base::c("DRB1", "DRB4", "DRB7", "DRB8", "DRB9"),
                               `08` = base::c("DRB1", "DRB9"),
                               `09` = base::c("DRB1", "DRB4", "DRB7", "DRB8", "DRB9"),
                               `10` = base::c("DRB1", "DRB6", "DRB9"),
                               `11` = base::c("DRB1", "DRB2", "DRB3", "DRB9"),
                               `12` = base::c("DRB1", "DRB2", "DRB3", "DRB9"),
                               `13` = base::c("DRB1", "DRB2", "DRB3", "DRB9"),
                               `14` = base::c("DRB1", "DRB2", "DRB3", "DRB9"),
                               `15` = base::c("DRB1", "DRB5", "DRB6", "DRB9"),
                               `16` = base::c("DRB1", "DRB5", "DRB6", "DRB9"))

# Infer the complete set of DRB genes expected for the individual by combining haplotype configurations from all observed DRB1 allele groups
DRB_genes <- base::unique(base::unlist(group_dictionary[allele_groups]))

# Define the full expected set of HLA class II genes, including DP, DQ, DR, with DRB paralogs inferred from DRB1 haplotype structure
all_HLAclassII_genes <- base::c("DPA1", "DPA2", "DPB1", "DPB2", "DQA1", "DQA2", "DQB1", "DQB2", "DRA", DRB_genes)

# Define canonical reference alleles to use when genotype calls are missing by selecting, for each locus, the first allele in the IPD-IMGT/HLA database with a complete genomic sequence
# NOTE: No full genomic sequences were available for DRB2, DRB6, DRB7, DRB8, or DRB9, and these loci are therefore excluded from this selection
canonical_reference_alleles <- base::list(DPA1 = "DPA1*01:03:01",
                                          DPA2 = "DPA2*01:01:01",
                                          DPB1 = "DPB1*01:01:01",
                                          DPB2 = "DPB2*01:01:01",
                                          DQA1 = "DQA1*01:01:01",
                                          DQA2 = "DQA2*01:01:01",
                                          DQB1 = "DQB1*02:01:01",
                                          DQB2 = "DQB2*01:01:01",
                                          DRA = "DRA*01:01:01",
                                          DRB1 = "DRB1*01:01:01",
                                          DRB3 = "DRB3*01:01:02",
                                          DRB4 = "DRB4*01:01:01",
                                          DRB5 = "DRB5*01:01:01")

# Identify HLA class II genes expected based on haplotype inference but missing from the genotype call set
missing_HLAclassII_genes <- all_HLAclassII_genes[!all_HLAclassII_genes %in% base::names(hla_genotype)]

# For HLA class II genes that are expected based on haplotype inference but missing from the genotype call set, add a fallback allele using the predefined canonical reference allele (when available) to ensure complete locus representation
for(gene in missing_HLAclassII_genes){if(gene %in% base::names(canonical_reference_alleles)){hla_genotype[[gene]] <- canonical_reference_alleles[[gene]]}}


####################################################################################################
# 5. Import HLA database, parse allele sequences, and build personalized chrHLA
####################################################################################################

# Print status message
base::message(base::paste0(base::format(x = base::Sys.time(), "%Y-%m-%d %H:%M:%S "), "Building chrHLA by concatenating IPD-IMGT/HLA allele sequences and writing its feature annotations, while masking the original HLA regions on chr6 and removing their annotations..."))

# Extract HLA database version from metadata in allele list file
hla_db_version <- base::readLines(con = glue::glue("{output_dir}/{IPD_IMGT_HLA_dir}/allelelist.txt"))
hla_db_version <- base::grep(pattern = "version:", x = hla_db_version, value = TRUE)
hla_db_version <- base::sub(pattern = "^# version: ", replacement = "", hla_db_version)
hla_db_version <- base::sub(pattern = " ", replacement = "_", hla_db_version)

# Load the full IPD-IMGT/HLA database as raw text
hla_db <- base::readLines(con = glue::glue("{output_dir}/{IPD_IMGT_HLA_dir}/hla.dat"))

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

# Read allele list file and rename columns
allele_list <- read.csv(file = glue::glue("{output_dir}/{IPD_IMGT_HLA_dir}/allelelist.txt"), header = TRUE, comment.char = "#")
base::colnames(allele_list) <- base::c("allele.id", "allele.name")

# Subset HLA database to retain only alleles present in the current allele list
hla_db <- hla_db[allele_list$allele.id]

# Identify full-length HLA alleles (contain UTRs and introns in expected pattern)
allele_list$full.length <- base::unlist(base::sapply(X = hla_db, FUN = function(record){
  annotation_lines <- base::grep(pattern = "^FT +(UTR|exon|intron)", x = record, value = TRUE)
  annotations <- base::gsub(pattern = "^FT +(UTR|exon|intron) +[0-9\\.]+$", replacement = "\\1", x = annotation_lines)
  if(annotations[1] == "UTR" && annotations[base::length(annotations)] == "UTR" && base::sum(annotations == "intron") >= base::sum(annotations == "exon") - 1){return(TRUE)}else{return(FALSE)}
}))

# Initialize container for concatenated HLA allele sequences
chr_HLA <- ""
# Create a working copy of the reference FASTA to be personalized
fasta_personalized <- fasta
# Initialize empty GTF structure to store coordinates for all inserted HLA alleles
gtf_HLA <- base::data.frame(seqnames = base::character(0), start = base::numeric(0), end = base::numeric(0), width = base::numeric(0), strand = base::character(0), phase = base::numeric(0), 
                            source = base::character(0), type = base::character(0),
                            gene_id = base::character(0), gene_type = base::character(0), gene_name = base::character(0), 
                            transcript_id = base::character(0), transcript_type = base::character(0), transcript_name = base::character(0),
                            exon_number = base::numeric(0), exon_id = base::character(0))

# Mask HLA class II gene loci in the reference FASTA and remove their annotations from the GTF to remove native HLA class II genes from the reference, allowing them to be reintroduced later as personalized, allele-resolved loci (e.g. on a new chromosome 'chrHLA') without sequence or annotation duplication
for(gene in base::c("HLA-DPA1", "HLA-DPA2", "HLA-DPB1", "HLA-DPB2",
                    "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DQB2",
                    "HLA-DRA", "HLA-DRB1", "HLA-DRB2", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5", "HLA-DRB6", "HLA-DRB7", "HLA-DRB8", "HLA-DRB9")){
  # Proceed only if the gene is present in the filtered GTF annotation
  if(gene %in% gtf_filtered$gene_name){
    # Identify the chromosome/contig containing this HLA gene in the personalized FASTA
    mask_seqname <- base::grep(pattern = base::unique(gtf_filtered[gtf_filtered$gene_name == gene, "seqnames"]), x = base::names(fasta_personalized), value = TRUE)
    # Determine the genomic start coordinate of the gene to be masked
    mask_start <- base::min(gtf_filtered[gtf_filtered$gene_name == gene, "start"])
    # Determine the genomic end coordinate of the gene to be masked
    mask_end <- base::max(gtf_filtered[gtf_filtered$gene_name == gene, "end"])
    # Replace the gene region in the personalized FASTA with Ns to mask the HLA sequence of the default sequence
    fasta_personalized[[mask_seqname]] <- Biostrings::replaceLetterAt(x = fasta_personalized[[mask_seqname]], at = mask_start:mask_end, letter = base::paste(base::rep(x = "N", mask_end - mask_start + 1), collapse = ""))
    # Remove the masked gene entries from the filtered GTF annotation
    gtf_filtered <- gtf_filtered[gtf_filtered$gene_name != gene, ]
  }
}

# Running genomic coordinate for placing each allele consecutively in chrHLA
start_pos <- 0

# Iterate over each HLA gene in the genotype file
for(gene in base::names(hla_genotype)){
  # Iterate over all unique alleles for this HLA gene
  for(allele in base::unique(hla_genotype[[gene]])){
    # If a full-length sequence is available...
    if(allele_list[base::grepl(pattern = allele, x = allele_list$allele.name, fixed = TRUE) , "full.length"][1]){
      # Retrieve (first) allele ID corresponding to the allele entry from the allel list
      allele_id <- allele_list[base::grepl(pattern = allele, x = allele_list$allele.name, fixed = TRUE) , "allele.id"][1]
    } 
    # If not...
    if(!allele_list[base::grepl(pattern = allele, x = allele_list$allele.name, fixed = TRUE) , "full.length"][1]){
      # Drop last field from the allele name
      lower_resolution <- base::paste(head(base::strsplit(x = allele, split = ":")[[1]], -1), collapse = ":")
      # Retrieve (first) allele ID corresponding to the lower-resolution allele from the allele list
      allele_id <- allele_list[base::grepl(pattern = lower_resolution, x = allele_list$allele.name, fixed = TRUE) , "allele.id"][1]
    }
    # Extract the full IMGT/HLA database record for this allele
    record <- hla_db[[allele_id]]
    # Extract the sequence lines (everything after 'SQ')
    seq_lines <- record[(which(grepl(pattern = "^SQ +", x = record)) + 1):base::length(record)]
    # Concatenate sequence and remove non-ACGT characters
    seq <- base::toupper(base::paste(base::gsub(pattern = "[^ACGTacgt]", replacement = "", seq_lines), collapse = ""))
    # Append allele sequence to chrHLA
    chr_HLA <- base::paste0(chr_HLA, seq)
    # Extract all annotation entries for UTRs, exons, and introns
    annotation_lines <- base::grep(pattern = "FT +UTR|FT +exon|FT +intron", x = record, value = TRUE)
    # Extract exon annotation lines and convert them into numeric coordinate ranges
    exon_annotations_lines <- base::grep(pattern = "exon", x = annotation_lines, value = TRUE)
    exon_annotations <- base::sub(pattern = "^FT +exon +([0-9\\.]+)$", replacement = "\\1", x = exon_annotations_lines)
    exon_annotations <- base::sapply(X = exon_annotations, FUN = function(exon) base::as.numeric(base::strsplit(x = exon, split = "..", fixed = TRUE)[[1]]), simplify = FALSE)
    base::names(exon_annotations) <- base::paste0("exon", 1:base::length(exon_annotations))
    # Extract CDS annotation block and convert joined segments into numeric coordinate ranges
    CDS_annotation_lines <- record[base::which(base::grepl(pattern = "^FT +CDS", x = record)):(base::which(base::grepl(pattern = "^FT +/codon_start=", x = record)) - 1)]
    CDS_annotation_line <- base::paste0(base::gsub(pattern = "^FT +CDS +join\\(|^FT +|\\)$", replacement = "", x = CDS_annotation_lines), collapse = "")
    CDS_annotations <- base::strsplit(x = CDS_annotation_line, split = ",", fixed = TRUE)[[1]]
    CDS_annotations <- base::sapply(X = CDS_annotations, FUN = function(CDS) base::as.numeric(base::strsplit(x = CDS, split = "..", fixed = TRUE)[[1]]), simplify = FALSE)
    base::names(CDS_annotations) <- sapply(CDS_annotations, function(CDS){base::names(exon_annotations)[sapply(exon_annotations, function(exon){CDS[1] >= exon[1] && CDS[2] <= exon[2]})]})
    # Extract UTR annotation lines convert them into numeric coordinate ranges
    UTR_annotation_lines <- base::grep(pattern = "UTR", x = annotation_lines, value = TRUE)
    UTR_annotations <- base::sub(pattern = "^FT +UTR +([0-9\\.]+)$", replacement = "\\1", x = UTR_annotation_lines)
    UTR_annotations <- base::sapply(X = UTR_annotations, FUN = function(UTR) base::as.numeric(base::strsplit(x = UTR, split = "..", fixed = TRUE)[[1]]), simplify = FALSE)
    base::names(UTR_annotations) <- base::rep(x = "UTR", base::length(UTR_annotations))
    # Extract the offset at which the first complete codon of a coding feature can be found
    codon_start_line <- base::grep(pattern = "codon_start", x = record, value = TRUE)
    codon_start <- base::as.numeric(base::sub(pattern = "^FT +/codon_start=(\\d{1,})$", replacement = "\\1", x = codon_start_line))
    # Build gene-level GTF entry
    gene_row <- base::data.frame(seqnames = "chrHLA", 
                                 start = base::min(base::as.numeric(base::unlist(UTR_annotations))), 
                                 end = base::max(base::as.numeric(base::unlist(UTR_annotations))), 
                                 width = NA, 
                                 strand = "+",
                                 phase = NA,
                                 source = hla_db_version, 
                                 type = "gene",
                                 gene_id = base::paste0("HLA-", gene), 
                                 gene_type = "protein_coding", 
                                 gene_name = base::paste0("HLA-", gene), 
                                 transcript_id = NA, 
                                 transcript_type = NA, 
                                 transcript_name = NA,
                                 exon_number = NA, 
                                 exon_id = NA)
    # Build transcript-level GTF entry
    transcript_row <- base::data.frame(seqnames = "chrHLA", 
                                       start = base::min(base::as.numeric(base::unlist(UTR_annotations))), 
                                       end = base::max(base::as.numeric(base::unlist(UTR_annotations))), 
                                       width = NA, 
                                       strand = "+",
                                       phase = NA,
                                       source = hla_db_version, 
                                       type = "transcript",
                                       gene_id = base::paste0("HLA-", gene), 
                                       gene_type = "protein_coding", 
                                       gene_name = base::paste0("HLA-", gene), 
                                       transcript_id = allele_id, 
                                       transcript_type = "protein_coding", 
                                       transcript_name = base::paste0("HLA-", base::gsub(pattern = "\\*|\\:", replacement = "_", allele)),
                                       exon_number = NA, 
                                       exon_id = NA)
    # Build CDS rows for all CDS intervals
    CDS_rows <- base::data.frame(seqnames = "chrHLA", 
                                 start = base::sapply(X = CDS_annotations, FUN = function(CDS) CDS[1]), 
                                 end = base::sapply(X = CDS_annotations, FUN = function(CDS) CDS[2]), 
                                 width = NA, 
                                 strand = "+",
                                 phase = base::c(0, base::sapply(2:base::length(CDS_annotations), function(x){(CDS_annotations[[x-1]][2] - CDS_annotations[[1]][1] + 1) %% 3})),
                                 source = hla_db_version, 
                                 type = "CDS",
                                 gene_id = base::paste0("HLA-", gene), 
                                 gene_type = "protein_coding", 
                                 gene_name = base::paste0("HLA-", gene), 
                                 transcript_id = allele_id, 
                                 transcript_type = "protein_coding", 
                                 transcript_name = base::paste0("HLA-", base::gsub(pattern = "\\*|\\:", replacement = "_", allele)),
                                 exon_number = base::sapply(X = base::names(CDS_annotations), FUN = function(CDS) base::as.numeric(base::gsub(pattern = "^exon", replacement = "", x = CDS))), 
                                 exon_id = base::paste(allele_id, base::names(CDS_annotations), sep = "_"))
    # Build exon rows for all exons
    exon_rows <- base::data.frame(seqnames = "chrHLA", 
                                  start = base::sapply(X = exon_annotations, FUN = function(exon) exon[1]), 
                                  end = base::sapply(X = exon_annotations, FUN = function(exon) exon[2]), 
                                  width = NA, 
                                  strand = "+",
                                  phase = NA,
                                  source = hla_db_version, 
                                  type = "exon",
                                  gene_id = base::paste0("HLA-", gene), 
                                  gene_type = "protein_coding", 
                                  gene_name = base::paste0("HLA-", gene), 
                                  transcript_id = allele_id, 
                                  transcript_type = "protein_coding", 
                                  transcript_name = base::paste0("HLA-", base::gsub(pattern = "\\*|\\:", replacement = "_", allele)),
                                  exon_number = base::sapply(X = base::names(exon_annotations), FUN = function(exon) base::as.numeric(base::gsub(pattern = "^exon", replacement = "", x = exon))), 
                                  exon_id = base::paste(allele_id, base::names(exon_annotations), sep = "_"))
    # Build start and stop codon row
    start_codon_row <- base::data.frame(seqnames = "chrHLA", 
                                        start = CDS_annotations[[1]][1], 
                                        end = CDS_annotations[[1]][1] + 2, 
                                        width = NA, 
                                        strand = "+",
                                        phase = 0,
                                        source = hla_db_version, 
                                        type = "start_codon",
                                        gene_id = base::paste0("HLA-", gene), 
                                        gene_type = "protein_coding", 
                                        gene_name = base::paste0("HLA-", gene), 
                                        transcript_id = allele_id, 
                                        transcript_type = "protein_coding", 
                                        transcript_name = base::paste0("HLA-", base::gsub(pattern = "\\*|\\:", replacement = "_", allele)),
                                        exon_number = base::as.numeric(base::gsub(pattern = "^exon", replacement = "", x = base::names(CDS_annotations)[1])), 
                                        exon_id = base::paste(allele_id, base::names(CDS_annotations)[1], sep = "_"))
    stop_codon_row <- base::data.frame(seqnames = "chrHLA", 
                                       start = CDS_annotations[[base::length(CDS_annotations)]][2] - 2, 
                                       end = CDS_annotations[[base::length(CDS_annotations)]][2], 
                                       width = NA, 
                                       strand = "+",
                                       phase = 0,
                                       source = hla_db_version, 
                                       type = "stop_codon",
                                       gene_id = base::paste0("HLA-", gene), 
                                       gene_type = "protein_coding", 
                                       gene_name = base::paste0("HLA-", gene), 
                                       transcript_id = allele_id, 
                                       transcript_type = "protein_coding", 
                                       transcript_name = base::paste0("HLA-", base::gsub(pattern = "\\*|\\:", replacement = "_", allele)),
                                       exon_number = base::as.numeric(base::gsub(pattern = "^exon", replacement = "", x = base::names(CDS_annotations)[base::length(CDS_annotations)])), 
                                       exon_id = base::paste(allele_id, base::names(CDS_annotations)[base::length(CDS_annotations)], sep = "_"))
    # Build UTR rows
    UTR_rows <- base::data.frame(seqnames = "chrHLA", 
                                 start = base::sapply(X = UTR_annotations, FUN = function(UTR) UTR[1]), 
                                 end = base::sapply(X = UTR_annotations, FUN = function(UTR) UTR[2]), 
                                 width = NA, 
                                 strand = "+",
                                 phase = 0,
                                 source = hla_db_version, 
                                 type = "UTR",
                                 gene_id = base::paste0("HLA-", gene), 
                                 gene_type = "protein_coding", 
                                 gene_name = base::paste0("HLA-", gene), 
                                 transcript_id = allele_id, 
                                 transcript_type = "protein_coding", 
                                 transcript_name = base::paste0("HLA-", base::gsub(pattern = "\\*|\\:", replacement = "_", allele)),
                                 exon_number = NA, 
                                 exon_id = NA)
    # Combine all annotation rows for this allele
    annotations <- base::rbind(gene_row, transcript_row, CDS_rows, exon_rows, start_codon_row, stop_codon_row, UTR_rows)
    rownames(annotations) <- NULL
    # Compute feature widths
    annotations$width <- annotations$end - annotations$start + 1
    # Offset all coordinates by current chrHLA position
    annotations$start <- annotations$start + start_pos
    annotations$end <- annotations$end + start_pos
    # Append allele annotations to chrHLA GTF
    gtf_HLA <- rbind(gtf_HLA, annotations)
    # Update running coordinate for next allele sequence
    start_pos <- start_pos + base::nchar(seq)
  }
}

# For each HLA gene, collapse multiple gene-level entries (one per allele) into a single gene annotation spanning all concatenated alleles
for(gene in base::unique(gtf_HLA$gene_name)){
  # Extract all gene-level rows for this gene
  gene_annotations <- gtf_HLA[gtf_HLA$gene_id == gene & gtf_HLA$type == "gene", ]
  # Remove the original (allele-specific) gene rows from the GTF
  gtf_HLA <- gtf_HLA[!(gtf_HLA$gene_id == gene & gtf_HLA$type == "gene"), ]
  # Recompute gene start and end to cover the full span of all alleles
  gene_annotations$start <- base::min(gene_annotations$start)
  gene_annotations$end <- base::max(gene_annotations$end)
  # Update gene width based on the new coordinates
  gene_annotations$width <- base::max(gene_annotations$end) - base::min(gene_annotations$start) + 1
  # Reinsert a single, non-duplicated gene entry into the GTF data frame
  gtf_HLA <- base::rbind(gtf_HLA, dplyr::distinct(gene_annotations, .keep_all = FALSE))
}

# Convert feature type to an ordered factor
gtf_HLA$type <- base::factor(x = gtf_HLA$type, levels = base::c("gene", "transcript", "CDS", "exon", "start_codon", "stop_codon", "UTR"), ordered = TRUE)
# Sort the GTF by genomic start coordinate, then by feature type order
gtf_HLA <- gtf_HLA[base::order(gtf_HLA$start, gtf_HLA$type), ]


####################################################################################################
# 6. Write personalized reference FASTA and GTF file
####################################################################################################

# Print status message
base::message(base::paste0(base::format(x = base::Sys.time(), "%Y-%m-%d %H:%M:%S "), "Saving personalized FASTA and GTF file with chrHLA..."))

# Add newly constructed chrHLA as a new chromosome to personalized FASTA
fasta_personalized[["chrHLA HLA"]] <- Biostrings::DNAString(x = chr_HLA)
# Write personalized GRCh38 genome FASTA file
Biostrings::writeXStringSet(x = fasta_personalized, filepath = glue::glue("{output_dir}/Homo_sapiens.GRCh38.dna.primary_assembly.personalized.fa"), format = "fasta", compress = FALSE)

# Add any missing GTF columns to HLA annotations
gtf_HLA[base::colnames(gtf_filtered)[!base::colnames(gtf_filtered) %in% base::colnames(gtf_HLA)]] <- NA
# Combine GRCh38 filtered GTF with chrHLA GTF entries
gtf_personalized <- base::rbind(gtf_filtered, gtf_HLA)
# Write personalized GTF file
rtracklayer::export(object = gtf_personalized, con = glue::glue("{output_dir}/Homo_sapiens.GRCh38.dna.primary_assembly.personalized.gtf"), format = "gtf")


####################################################################################################
# 7. Finish and cleanup source directories
####################################################################################################

# Print status message
base::message(base::paste0(base::format(x = base::Sys.time(), "%Y-%m-%d %H:%M:%S "), "Removing temporary Gencode and IPD-IMGT/HLA directories. Personalized reference build complete."))

# Delete source directories and all their contents used for downloading and processing Gencode and IPD-IMGT/HLA data files (this only keeps the final personalized FASTA and GTF outputs)
base::unlink(x = glue::glue("{output_dir}/{gencode_dir}"), recursive = TRUE, force = TRUE)
base::unlink(x = glue::glue("{output_dir}/{IPD_IMGT_HLA_dir}"), recursive = TRUE, force = TRUE)
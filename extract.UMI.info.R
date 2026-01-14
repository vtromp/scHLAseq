# 
allele_list <- "analysis/data/Allelelist.txt"

#
default_molecule_info_file <- "analysis/data/default_reference/DRneg/molecule_info.h5"
personalized_molecule_info_file <- "analysis/data/personalized_reference/DRneg/molecule_info.h5"

# Define a function to convert a 2-bit-encoded number to an UMI sequence
decompress.UMI <- function(UMI){
  # Convert compressed UMI to binary string, zero-padded to 24 digits
  binary_str <- base::paste0(base::rev(x = base::as.integer(base::intToBits(UMI)[1:24])), collapse = "")
  # Split binary string into chunks of 2
  nuc_list <- substring(binary_str, seq(1, nchar(binary_str), 2), seq(2, nchar(binary_str), 2))
  # Map binary pairs to nucleotides
  UMI_seq <- base::sapply(X = nuc_list, FUN = function(i){if(i == "00"){return("A")}else if(i == "01"){return("C")}else if(i == "10"){return("G")}else{return("T")}})
  # Collapse into a single string and return
  return(base::paste0(UMI_seq, collapse = ""))
}

# Read zero-based indices into the barcode list, indicating the cell-barcode assigned to this transcript
default_barcode_indices <- rhdf5::h5read(file = default_molecule_info_file, name = "barcode_idx")
personalized_barcode_indices <- rhdf5::h5read(file = personalized_molecule_info_file, name = "barcode_idx")
# Read list of barcodes
default_barcodes <- rhdf5::h5read(file = default_molecule_info_file, name = "barcodes")
personalized_barcodes <- rhdf5::h5read(file = personalized_molecule_info_file, name = "barcodes")
# Read zero-based indices into the gene list, indicating the the gene to which each transcript is assigned
default_feature_indices <- rhdf5::h5read(file = default_molecule_info_file, name = "feature_idx")
personalized_feature_indices <- rhdf5::h5read(file = personalized_molecule_info_file, name = "feature_idx")
# Read list of gene names and gene IDs
default_feature_names <- rhdf5::h5read(file = default_molecule_info_file, name = "features/name")
personalized_feature_names <- rhdf5::h5read(file = personalized_molecule_info_file, name = "features/name")
default_feature_ids <- rhdf5::h5read(file = default_molecule_info_file, name = "features/id")
personalized_feature_ids <- rhdf5::h5read(file = personalized_molecule_info_file, name = "features/id")
# Read 2-bit encoded UMI sequence of each transcript
default_umis <- rhdf5::h5read(file = default_molecule_info_file, name = "umi")
personalized_umis <- rhdf5::h5read(file = personalized_molecule_info_file, name = "umi")

# Create dataframes combining barcode, UMI, gene name, and gene ID
df_default <- base::data.frame(barcode = default_barcodes[default_barcode_indices + 1],
                               umi = default_umis, 
                               gene.name = default_feature_names[default_feature_indices + 1],
                               gene.id = default_feature_ids[default_feature_indices + 1])
df_personalized <- base::data.frame(barcode = personalized_barcodes[personalized_barcode_indices + 1],
                                    umi = personalized_umis,
                                    gene.name = personalized_feature_names[personalized_feature_indices + 1],
                                    gene.id = personalized_feature_ids[personalized_feature_indices + 1])

# Remove intermediate objects to free memory
rm(default_barcode_indices)
rm(personalized_barcode_indices)
rm(default_barcodes)
rm(personalized_barcodes)
rm(default_feature_indices)
rm(personalized_feature_indices)
rm(default_feature_names)
rm(personalized_feature_names)
rm(default_feature_ids)
rm(personalized_feature_ids)
rm(default_umis)
rm(personalized_umis)

# Read allele list file and rename columns
allele_list <- read.csv(file = allele_list, header = TRUE, comment.char = "#")
base::colnames(allele_list) <- base::c("allele.id", "allele.name")

# Replace HLA gene names in the personalized data frame with the corresponding allele names from the allele list (which only apply to rows where 'gene.id' starts with "HLA")
df_personalized[base::grepl(pattern = "^HLA", x = df_personalized$gene.id), "gene.name"] <- base::sapply(X = base::grep(pattern = "^HLA", x = df_personalized$gene.id, value = TRUE), FUN = function(id){base::paste0("HLA-", allele_list[allele_list$allele.id == id, "allele.name"])})

# Create a combined dataframe of HLA-related transcripts from both datasets
df_combined <- base::rbind(df_default[base::grepl(pattern = "^HLA-", x = df_default$gene.name), c("barcode", "umi")], df_personalized[base::grepl(pattern = "^HLA-", x = df_personalized$gene.name), c("barcode", "umi")])
# Remove duplicate rows (same barcode + UMI)
df_combined <- df_combined[!base::duplicated(x = df_combined), ]

# Merge gene names from the default reference into the combined dataframe by matching barcode and UMI and rename the merged 'gene.name' column to "default.reference"
df_combined <- base::merge(x = df_combined, y = df_default[, base::c("barcode", "umi", "gene.name")], by = base::c("barcode", "umi"), all.x = TRUE)
base::colnames(df_combined)[base::colnames(df_combined) == "gene.name"] <- "default.reference"
# Merge gene names from the personalized reference into the combined dataframe by matching barcode and UMI and rename the merged 'gene.name' column to "personalized.reference"
df_combined <- base::merge(x = df_combined, y = df_personalized[, base::c("barcode", "umi", "gene.name")], by = base::c("barcode", "umi"), all.x = TRUE)
base::colnames(df_combined)[base::colnames(df_combined) == "gene.name"] <- "personalized.reference"

# Decompress the 2-bit-encoded UMIs into nucleotide sequences
df_combined$umi <- base::sapply(X = df_combined$umi, FUN = decompress.UMI)

# Replace NA values with "unmapped" for transcripts that were not assigned to a gene
df_combined$default.reference[base::is.na(df_combined$default.reference)] <- "unmapped"
df_combined$personalized.reference[base::is.na(df_combined$personalized.reference)] <- "unmapped"
# Assign "other" to gene names that are not classical HLA genes
df_combined$default.reference[!base::grepl(pattern = "^HLA-|unmapped", x = df_combined$default.reference)] <- "other"
df_combined$personalized.reference[!base::grepl(pattern = "^HLA-|unmapped", x = df_combined$personalized.reference)] <- "other"

write.csv(df_combined, "DRneg.csv", row.names = FALSE)
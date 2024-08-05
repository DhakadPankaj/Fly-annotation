
library(cubar)
suppressPackageStartupMessages(library(Biostrings))
# Create a directory to store files
output_dir <- "CUB_output"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Read the input sequences
# Get the file path as an input argument
file_path <- commandArgs(trailingOnly = TRUE)[1]
ref_path <- commandArgs(trailingOnly = TRUE)[2]


# Read the sequences from the files
sequences <- readDNAStringSet(file_path)
ref_sequences <- readDNAStringSet(ref_path)

seq_qc <- check_cds(sequences, min_len = 30, check_len = TRUE, check_start = TRUE, check_stop = TRUE, check_istop = TRUE, 
rm_start = TRUE, rm_stop = TRUE, start_codons = c("ATG"))

ref_qc <- check_cds(ref_sequences, min_len = 30, check_len = TRUE, check_start = TRUE, check_stop = TRUE, check_istop = TRUE, 
rm_start = TRUE, rm_stop = TRUE, start_codons = c("ATG"))

# Count the codons
cf_all <- count_codons(seq_qc)
cf_ref <- count_codons(ref_qc)
rscu_heg <- est_rscu(cf_ref)

cai <- get_cai(cf_all)
gc <- get_gc(cf_all)
gc3 <- get_gc3s(cf_all)
gc4d <- get_gc4d(cf_all)
enc <- get_enc(cf_all)

file_name <- basename(file_path)
result_file <- file.path(output_dir, gsub(".fa", "_results.tsv", file_name))

# Write all the results to a single file
result_data <- data.frame(gene = names(gc), GC = gc, GC3 = gc3, GC4d = gc4d, ENC = enc, CAI = cai)
write.table(result_data, result_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, append = FALSE)

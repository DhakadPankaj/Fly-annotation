
library(cubar)
library(parallel)
suppressPackageStartupMessages(library(Biostrings))
# Create a directory to store files
output_dir <- "CUB_output"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

species_list <- readLines("species_list.txt") 

process_species <- function(species) {
  # Read the input sequences
  # Get the file path as an input argument
  file_path <- paste0("CDS_trimmed/", species, ".fa")
  #ref_path <- paste0("HEG_genes/", species, "_highly_expressed.fa")

# Read the sequences from the files
sequences <- readDNAStringSet(file_path)
#ref_sequences <- readDNAStringSet(ref_path)

#seq_qc <- check_cds(sequences, min_len = 30, check_len = TRUE, check_start = TRUE, check_stop = TRUE, check_istop = TRUE, 
#rm_start = TRUE, rm_stop = TRUE, start_codons = c("ATG"))

#ref_qc <- check_cds(ref_sequences, min_len = 30, check_len = TRUE, check_start = TRUE, check_stop = TRUE, check_istop = TRUE, 
#rm_start = TRUE, rm_stop = TRUE, start_codons = c("ATG"))

seq_qc <- sequences
#ref_qc <- ref_sequences

# Count the codons
cf_all <- count_codons(seq_qc)
#cf_ref <- count_codons(ref_qc)

# Revome stop codons
#cf_all <- as.data.frame(cf_all) %>% select(-TAA, -TAG, -TGA)

# save codon frequencies
cf_all_file <- file.path(output_dir, paste0(species, "_cf_all.tsv"))
#write.table(cf_all, cf_all_file, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

cf_ref_file <- file.path(output_dir, paste0(species, "_cf_ref.tsv"))
#write.table(cf_ref, cf_ref_file, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

#rscu_heg <- est_rscu(cf_ref)
# save RSCU of HEG genes
#rscu_heg_file <- file.path(output_dir, paste0(species, "_rscu_heg.tsv"))
#write.table(rscu_heg, rscu_heg_file, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

#cai <- get_cai(cf_all, rscu_heg)
gc <- get_gc(cf_all)
gc3 <- get_gc3s(cf_all)
gc4d <- get_gc4d(cf_all)
enc <- get_enc(cf_all)
#optimal_codons <- est_optimal_codons(cf_ref)
#fop <- get_fop(seq_qc)

file_name <- basename(file_path)
result_file <- file.path(output_dir, gsub(".fa", "_results.tsv", file_name))

# Write all the results to a single file
result_data <- data.frame(gene = names(gc), GC = gc, GC3 = gc3, GC4d = gc4d, ENC = enc)
#write.table(result_data, result_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, append = FALSE)

# Return the mean CAI value and mean ENC value for the species
results <- data.frame(species = species, mean_gc = mean(gc)*100, mean_gc3 = mean(gc3)*100, mean_gc4d = mean(gc4d)*100, mean_enc = mean(enc))

return(results)

}

# Process for each species in parallel in "Species_list.txt"
num_cores <- 40
results <- mclapply(species_list, process_species, mc.cores = num_cores)

results_df <- do.call(rbind, results)

#mean_fop <- sapply(results, function(x) x[3])
# Export the mean CAI and ENC values for each species
write.table(results_df, file = "GC_ENC.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

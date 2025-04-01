library(cubar)
suppressPackageStartupMessages(library(Biostrings))
# Create a directory to store files
output_dir <- "CUB_output"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

#CDS_path <- commandArgs(trailingOnly = TRUE)[1]

species_list_path <- "species_list.txt"
species_list <- readLines(species_list_path)


library(parallel)

# Function to process each species
process_species <- function(species) {
    tryCatch({
    CDS_path = "HOGs_seq/"
    species_CDS_path <- paste0(CDS_path, species, ".fa")
    sequences <- readDNAStringSet(species_CDS_path)
    
    seq_qc <- check_cds(sequences, min_len = 30, check_len = TRUE, check_start = TRUE, check_stop = TRUE, check_istop = TRUE, 
                                            rm_start = TRUE, rm_stop = TRUE, start_codons = c("ATG"))
    
    cf_all <- count_codons(seq_qc)
    
    rscu_all <- est_rscu(cf_all)

    # Convert the RSCU values to a data frame and add the species name
    rscu_df <- as.data.frame(rscu_all)
    rscu_df$species <- species
    
    return(rscu_df)
     }, error = function(e) {
        message(sprintf("Error processing species %s: %s", species, e$message))
        return(NULL)
    })
}

# Run the process in parallel
num_cores <- 40
results <- mclapply(species_list, process_species, mc.cores = num_cores)

# Combine results into a single dataframe
combined_rscu_df <- do.call(rbind, results)

# Reshape the dataframe to have species names as rows and codons as columns
library(reshape2)
final_rscu_df <- dcast(combined_rscu_df, species ~ codon, value.var = "RSCU")

# Write the combined results to a file
write.table(combined_rscu_df, file = file.path(output_dir, "rscu_all_combined.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(final_rscu_df, file = file.path(output_dir, "rscu_all_species.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
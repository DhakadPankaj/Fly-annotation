# Load the required packages
library(seqinr)
library(coRdon)
library(parallel)

# Create a directory to store files
output_dir <- "CUB_output"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

species_list <- readLines("species_list.txt") 

process_species <- function(species) {
    tryCatch({
  # Read the input sequences
  # Get the file path as an input argument
  file_path <- paste0("CDS_trimmed/", species, ".fa")
  ref_path <- paste0("HEG_genes/", species, "_highly_expressed.fa")
  # Read the input sequences
  sequences <- readSet(file=file_path)
  
  # Calculate codon usage frequencies
  codon_table <- codonTable(sequences)

  
  # Calculate codon adaptation index (CAI)
  ref_seq <- readSet(file=ref_path) # highly expressed genes in Dmel
  ref_table <- codonTable(ref_seq)

  cai <- CAI(codon_table, subsets = list(ref_table), ribosomal = FALSE, id_or_name2 = "1", alt.init = FALSE, stop.rm = TRUE,
  len.threshold = 300)

  fop <- Fop(codon_table, subsets = list(ref_table), ribosomal = FALSE, id_or_name2 = "1", alt.init = FALSE, stop.rm = TRUE,
    len.threshold = 300)

    fop <- as.data.frame(fop)
    #rownames(fop) <- getID(codon_table)

  enc <- ENC(codon_table, len.threshold = 300)
    enc <- as.data.frame(enc)
    rownames(enc) <- getID(codon_table)
  cai <- as.data.frame(cai)
    rownames(cai) <- getID(codon_table)
    
    # merge the two data frames
    df <- merge(cai, enc, by = "row.names")
    colnames(df) <- c("gene", "CAI", "ENC")
    fop$gene <- getID(codon_table)
    colnames(fop) <- c("FOP", "gene")
    df <- merge(df, fop, by = "gene")
    
    # Export the results
    #write.table(df, file.path(output_dir, paste0(species, "_cai_enc.tsv")), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

    # Return the mean CAI value and mean ENC value for the species
    return(c(mean(df$CAI), mean(df$ENC), mean(df$FOP)))

    }, error = function(e) {
        message(sprintf("Error processing species %s: %s", species, e$message))
        return(NULL)
    })
}


# Process for each species in parallel in "Species_list.txt"
num_cores <- 40
results <- mclapply(species_list, process_species, mc.cores = num_cores)

# Combine the results into a single data frame
mean_cai <- sapply(results, function(x) x[1])
mean_enc <- sapply(results, function(x) x[2])
mean_fop <- sapply(results, function(x) x[3])

# Export the mean CAI and ENC values for each species
write.table(data.frame(species = species_list, mean_cai = mean_cai, mean_enc = mean_enc, mean_fop = mean_fop), file.path(output_dir, "mean_cai_enc_fop_coRdon.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

#cai_output <- paste0(output_dir, "/", basename(file_path), ".cai")
#write.table(cai, cai_output, sep = "\t", row.names = TRUE, col.names = FALSE, quote = FALSE)

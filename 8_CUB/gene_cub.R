library(data.table)
library(parallel)
suppressPackageStartupMessages(library(dplyr))

num_cores <- 20

hog <- read.table("../Results_MSA/HOG_filtered.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

species_list <- readLines("species_list.txt") 

process_species <- function(species) {
    gc_file <- paste0("CUB_output/", species, "_gc.tsv")
    gc3_file <- paste0("CUB_output/", species, "_gc3.tsv")
    gc4d_file <- paste0("CUB_output/", species, "_gc4d.tsv")
    enc_file <- paste0("CUB_output/", species, "_enc.tsv")
    cai_file <- paste0("CUB_output/", species, "_cai.tsv")
    gc <- read.table(gc_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    gc3 <- read.table(gc3_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    gc4d <- read.table(gc4d_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    enc <- read.table(enc_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    cai <- read.table(cai_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    hog[[species]] <- strsplit(hog[[species]], ",\\s*")
    hog <- tidyr::unnest(hog, species)
    hog_subset <- hog %>% select(all_of(c("HOG", species)))
    # Convert the data frames to data tables for faster manipulation
    gc <- as.data.table(gc)
    gc3 <- as.data.table(gc3)
    gc4d <- as.data.table(gc4d)
    enc <- as.data.table(enc)
    cai <- as.data.table(cai)
    hog_subset <- as.data.table(hog_subset)
    output <- system2("bioawk", args = c("-c", "fastx", "'{ print $name, length($seq) }'", paste0("/data/home/s2215768/fly_annotation/complement_annotations/transcripts/", species, ".fa")), stdout = TRUE)
    output_df <- read.table(text = output, header = FALSE, stringsAsFactors = FALSE)
    cds_len <- as.data.table(output_df)
    # Set the first column as the key for each data table
    setkey(gc, V1)
    setkey(gc3, V1)
    setkey(gc4d, V1)
    setkey(enc, V1)
    setkey(cai, V1)
    setkey(cds_len, V1)
    # Join the data tables
    result_table <- gc[gc3][gc4d][enc][cai][cds_len]
    setnames(result_table, old = c("V1", "V2", "i.V2", "i.V2.1", "i.V2.2", "i.V2.3", "i.V2.4"), new = c("Gene", "GC", "GC3", "GC4D", "ENC", "CAI", "CDS_length"))
    final_table <- merge(hog_subset, result_table, by.x = species, by.y = "Gene", all.x = TRUE)
    final_table$Species <- species
    # Write the result table to a file
    result_file <- paste0("CUB_output/", species, "_result_table.tsv")
    fwrite(final_table, file = result_file, sep = "\t", quote = FALSE, row.names = FALSE)
    # save hog_subset
    #hog_file <- paste0("CUB_output/", species, "_hog_subset.tsv")
    #fwrite(hog_subset, file = hog_file, sep = "\t", quote = FALSE, row.names = FALSE)
}

results <- mclapply(species_list, process_species, mc.cores = num_cores)
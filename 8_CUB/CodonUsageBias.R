# Load the required packages
library(seqinr)
library(coRdon)

# Create a directory to store files
output_dir <- "CUB_output"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Read the input sequences
# Get the file path as an input argument
file_path <- commandArgs(trailingOnly = TRUE)[1]

# Read the input sequences
sequences <- readSet(file=file_path)

# Calculate codon usage frequencies
codon_table <- codonTable(sequences)
#codon_counts <- codonCounts(codon_table)

enc <- ENC(codon_table, len.threshold = 80)
enc <- as.data.frame(enc)
rownames(enc) <- getID(codon_table)

# Perform relative synonymous codon usage (RSCU) analysis
rscu <- lapply(read.fasta(file = file_path), uco, index = "rscu", as.data.frame = FALSE)
# Transform the RSCU table
rscu <- as.data.frame(rscu)
# Create a data frame from rscu
rscu_df <- as.data.frame(t(rscu))
#colnames(rscu_df) <- rownames(rscu)
rownames(rscu_df) <- NULL
# Create a new column as the first column and get the values from colnames(rscu)
rscu_df <- cbind(colnames(rscu), rscu_df)
colnames(rscu_df) <- c("gene", rownames(rscu))


# Calculate GC content for each sequence
gc_content <- sapply(sequences, function(x) {
    seqs <- tolower(paste(x, collapse = ""))
    seqs <- s2c(seqs) # Convert sequence to character vector
    gc_content <- GC(seqs)
    gc1 <- ifelse(length(seqs) >= 1, GC1(seqs), NA)
    gc2 <- ifelse(length(seqs) >= 2, GC2(seqs), NA)
    gc3 <- ifelse(length(seqs) >= 3, GC3(seqs), NA)
    return(c(gc_content, gc1, gc2, gc3))
})


# Convert to data frame
gc_content_df <- as.data.frame(t(gc_content))
colnames(gc_content_df) <- c("gc_content", "gc1", "gc2", "gc3")


# Get the file name from file_path
file_name <- basename(file_path)
# Output the results
write.table(rscu_df, file.path(output_dir, gsub(".fa", "_rscu.tsv", file_name)), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(enc, file.path(output_dir, gsub(".fa", "_enc.tsv", file_name)), sep = "\t", row.names = TRUE, col.names = FALSE, quote = FALSE)
write.table(gc_content_df, file.path(output_dir, gsub(".fa", "_gc.tsv", file_name)), sep = "\t", row.names = TRUE, col.names = FALSE, quote = FALSE)

# Load the required packages
library(seqinr)
library(coRdon)

# Create a directory to store files
output_dir <- "CUB_output"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Read the input sequences
# Get the file path as an input argument
file_path <- commandArgs(trailingOnly = TRUE)[1]
ref_path <- commandArgs(trailingOnly = TRUE)[2]

# Read the input sequences
sequences <- readSet(file=file_path)

# Calculate codon usage frequencies
codon_table <- codonTable(sequences)

# Calculate codon adaptation index (CAI)
ref_seq <- readSet(file=ref_path) # highly expressed genes in Dmel
ref_table <- codonTable(ref_seq)

cai <- CAI(codon_table, subsets = list(ref_table), ribosomal = FALSE, id_or_name2 = "11", alt.init = TRUE, stop.rm = TRUE)

cai <- as.data.frame(cai)
rownames(cai) <- getID(codon_table)

cai_output <- paste0(output_dir, "/", basename(file_path), ".cai")
write.table(cai, cai_output, sep = "\t", row.names = TRUE, col.names = FALSE, quote = FALSE)

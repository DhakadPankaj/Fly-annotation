# Make sure "Nes.R" is present in the current working directory.


# Functions to calculate the strength of selected codon usage, S = 4*Ne*s
# (2*Ne*s haploids) for amino acids encoded by only two codons.

# The function takes a codon frequency table

library(cubar)
library(parallel)
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(dplyr))

# Function to calculate the strength of selected codon usage, S = 4*Ne*s

process_species <- function(species) {
    tryCatch({
        # Run codonM in bash to get the codon frequency table (bash command line)
        # Ex. Nes/codonM CDS_trimmed/DROSOPHILA_MELANOGASTER.fa CDS_trimmed/DROSOPHILA_MELANOGASTER.m
        #system(paste0("Nes/codonM CDS_trimmed/", species, ".fa CDS_trimmed/", species, ".m"))

        #m <- matrix(scan(paste0("CDS_trimmed/", species, ".m")), ncol=61, byrow=T)
        #codons <- scan("Nes/codons", what="char")
        #colnames(m) <- codons
        
        # Codon frequency table
        sequences <- readDNAStringSet(paste0("CDS_trimmed/", species, ".fa"))
        cf_all <- count_codons(sequences)
        
        # Remove stop codons
        cf_all <- as.data.frame(cf_all) %>% select(-TAA, -TAG, -TGA)
     
        gene_exp <- read.table(paste0("gene_expression/", species, "_genes_HOGs_FPKM.txt"), header=TRUE, sep="\t")

        # Combine gene expression data with codon frequency table
        cf_all$Gene <- rownames(cf_all)
        m.cat <- merge(cf_all, gene_exp[,c("Gene", "Category")], by = "Gene")
        m.cat <- arrange(m.cat, Category)

        source("Nes.R")   
        # Bootstrap over codon frequency table (m)
        boot_results <- replicate(10000, {
            # Sample with replacement from the codon frequency table
            m.sample <- m.cat %>% select(-Gene) %>% group_by(Category)  %>% sample_n(size= n(), replace = TRUE) %>% ungroup()
            # Aggregate (sum) the sampled values according to the Category column
            m.aggregated <- m.sample %>%
                group_by(Category) %>%
                summarize(across(everything(), sum, na.rm = TRUE))

            GC_codons <- match(c("TTC", "TAC", "TGC", "CAC", "CAG", "AAC", "AAG", "GAC", "GAG"), colnames(as.data.frame(m.aggregated %>% select(-Category))))
            AT_codons <-  match(c("TTT", "TAT", "TGT", "CAT", "CAA", "AAT", "AAA", "GAT", "GAA"), colnames(as.data.frame(m.aggregated %>% select(-Category))))  
            Nes(as.matrix(m.aggregated %>% select(-Category)), GC_codons, AT_codons, ref=1, mean=T)
        })

        # Mean S value (From actual data)
        m.cat <- m.cat %>% select(-Gene) %>% group_by(Category) %>% summarize(across(everything(), sum, na.rm = TRUE))
        GC_codons <- match(c("TTC", "TAC", "TGC", "CAC", "CAG", "AAC", "AAG", "GAC", "GAG"), colnames(as.data.frame(m.cat %>% select(-Category))))
        AT_codons <-  match(c("TTT", "TAT", "TGT", "CAT", "CAA", "AAT", "AAA", "GAT", "GAA"), colnames(as.data.frame(m.cat %>% select(-Category))))
        S <- Nes(as.matrix(m.cat %>% select(-Category)), GC_codons, AT_codons, ref=1, mean=T)
        S.mean <- S["Mean", 20]
        # CI values for S
        CI <- quantile(boot_results["Mean", , ][20,], c(0.025, 0.975))
        
        # Store species, S.mean, CI[1], CI[2] as a dataframe, with round up to 3 decimal places
        results <- data.frame(species=species, S=round(S.mean, 3), CI1=round(as.numeric(CI[1]), 3), CI2=round(as.numeric(CI[2]), 3))

        # Get Nes and CI values for all 9 optimal codons
        #results <- lapply(1:length(op), function(i) {
        #    codon_pair <- c(codons[op[i]], codons[nop[i]])
        #    res <- Nes.boot(m.cat[40, codon_pair], m.cat[1, codon_pair], N=10000, alpha=0.05)
            # Store codons[op[i]], codons[nop[i]], S, CI[1], CI[2] as a dataframe
        #    return(data.frame(species=species,op=codons[op[i]], nop=codons[nop[i]], S=res$S[[1]], CI1=res$CI[[1]], CI2=res$CI[[2]]))
        #})

        return(results)

     }, error = function(e) {
        message(sprintf("Error processing species %s: %s", species, e$message))
        return(NULL)
    })
}   

# Run the process in parallel
num_cores <- 60
species_list_path <- "species_list.txt"
species_list <- readLines(species_list_path)

results <- mclapply(species_list, process_species, mc.cores = num_cores)

# Combine results into a single dataframe
results_df <- do.call(rbind, results)

# Write the combined results to a file
write.table(results_df, file = "S_values.tsv", sep = "\t", quote = FALSE, row.names = FALSE)


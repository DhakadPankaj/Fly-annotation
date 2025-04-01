#!/bin/bash

# Function to calculate RSCU from cusp output
calculate_rscu() {
    species=$1
    cusp_output="rscu/${species}.rscu"
    rscu_output="rscu/${species}_rscu.txt"

    # Parse the cusp output and calculate RSCU
    awk '
    BEGIN {
        FS = "\t";
        OFS = "\t";
        # Initialize the number of synonymous codons for each amino acid
        aa_synonyms["A"] = 4;
        aa_synonyms["C"] = 2;
        aa_synonyms["D"] = 2;
        aa_synonyms["E"] = 2;
        aa_synonyms["F"] = 2;
        aa_synonyms["G"] = 4;
        aa_synonyms["H"] = 2;
        aa_synonyms["I"] = 3;
        aa_synonyms["K"] = 2;
        aa_synonyms["L"] = 6;
        aa_synonyms["M"] = 1;
        aa_synonyms["N"] = 2;
        aa_synonyms["P"] = 4;
        aa_synonyms["Q"] = 2;
        aa_synonyms["R"] = 6;
        aa_synonyms["S"] = 6;
        aa_synonyms["T"] = 4;
        aa_synonyms["V"] = 4;
        aa_synonyms["W"] = 1;
        aa_synonyms["Y"] = 2;
    }
    {
        if ($1 !~ /^#/) {
            codon = $1;
            aa = $2;
            number = $5;
            aa_count[aa] += number;
            codon_count[codon] = number;
            codon_aa[codon] = aa;
        }
    }
    END {
        for (codon in codon_count) {
            aa = codon_aa[codon];
            if (aa_count[aa] > 0) {
                rscu = codon_count[codon] / (aa_count[aa] / aa_synonyms[aa]);
            } else {
                rscu = "NA";
            }
            print codon, aa, rscu;
        }
    }
    ' $cusp_output > $rscu_output
}

# Ensure the output directory exists
mkdir -p rscu

calculate_rscu ZAPRIONUS_VITTIGER

# Calculate RSCU for each species
#for species in $(cut -f1 ../genome_info_selected.tsv); do
#    calculate_rscu $species
#done
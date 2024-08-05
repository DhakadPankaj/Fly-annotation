#!/bin/bash

# This script calculates the Codon Usage Bias (CUB) for each species in a directory of CDS fasta files. 
# It uses the CodonUsageBias.R script to calculate GC content, ENC, and CAI for each species.
# The results are saved in the CUB_results directory.

threads=30
CDS_dir=$1

# Check if input argument is provided
if [ -z "$CDS_dir" ]; then
    echo "Input argument is missing. Please provide the CDS directory."
    exit 1
fi
# List of species
basename -a $CDS_dir/*.fa | sed 's/.fa//' > species_list.txt

# Function to run CodonUsageBias.R on a species
run_codon_usage_bias() {
    species=$1
    CDS_dir=$2
    if [ ! -s CUB_output/${species}_gc.tsv ]; then
        echo "Calculating CUB for ${species}"
       #Rscript CodonUsageBias.R ${CDS_dir}/${species}.fa
       Rscript cubar_CAI.R ${CDS_dir}/${species}.fa Highly_exp_genes/${species}.heg.fa
    fi
}

# Export the function so that it can be used by parallel
export -f run_codon_usage_bias

# Run CodonUsageBias.R on each species in parallel
parallel -j $threads run_codon_usage_bias {1} {2} :::: species_list.txt ::: $CDS_dir


get_GC_noncoding() {
    species=$1
    assembly=$2
    genome_dir="/data/home/s2215768/fly_annotation/data/CAT_genomes/selected_genomes"
    if [ ! -s CUB_output/${species}_gc_noncoding.tsv ]; then
        echo "Calculating GC % non-codingp part and whole genome for ${species}"
        bash get_GC_noncoding.sh $species $assembly $genome_dir
    fi
}

export -f get_GC_noncoding
parallel -j $threads --colsep "\t" get_GC_noncoding {1} {2} :::: genome_info_selected.tsv


# Combine the results into a single file
echo -e "Species\tGC\tGC1\tGC2\tGC3\tGC_genomic\tGC_nonCoding\tENC\tENC_range\tCAI\tCAI_range" > CUB_results.txt
get_CUB_results() {
    species=$1
    tmp_dir=$2
    GC=$(awk '{print $2}' CUB_output/${species}_gc.tsv | awk '{sum+=$1} END{printf "%.2f", (sum/NR)*100}')
    #GC1=$(awk '{print $3}' CUB_output/${species}_gc.tsv | awk '{sum+=$1} END{printf "%.2f", (sum/NR)*100}')
    #GC2=$(awk '{print $4}' CUB_output/${species}_gc.tsv | awk '{sum+=$1} END{printf "%.2f", (sum/NR)*100}')
    GC3=$(awk '{print $2}' CUB_output/${species}_gc3.tsv | awk '{sum+=$1} END{printf "%.2f", (sum/NR)*100}')
    GC4d=$(awk '{print $2}' CUB_output/${species}_gc4d.tsv | awk '{sum+=$1} END{printf "%.2f", (sum/NR)*100}')
    GC_genomic=$(awk 'NR>1 {print $2}' CUB_output/${species}_gc_whole_genome.tsv | awk '{sum+=$1} END{printf "%.2f", (sum/NR)*100}')
    GC_nonCoding=$(awk 'NR>1 {print $2}' CUB_output/${species}_gc_noncoding.tsv | awk '{sum+=$1} END{printf "%.2f", (sum/NR)*100}')
    Avg_ENC=$(awk '{print $2}' CUB_output/${species}_enc.tsv | awk '{sum+=$1} END{printf "%.2f", sum/NR}')
    ENC_range=$(awk '{print $2}' CUB_output/${species}_enc.tsv | sort -n | awk 'NR==1{min=$1} END{printf "%.2f-%.2f", min, $1}')
    Avg_CAI=$(awk '{print $2}' CUB_output/${species}_cai.tsv | awk '{sum+=$1} END{printf "%.2f", (sum/NR)}')
    CAI_range=$(awk '{print $2}' CUB_output/${species}_cai.tsv | sort -n | awk 'NR==1{min=$1} END{printf "%.2f-%.2f", min, $1}')
    echo -e "${species}\t${GC}\t${GC3}\t${GC4d}\t${GC_genomic}\t${GC_nonCoding}\t${Avg_ENC}\t${ENC_range}\t${Avg_CAI}\t${CAI_range}" > $tmp_dir/${species}_cub.txt
}

tmp_dir=$(mktemp -d)
export -f get_CUB_results
parallel -j $threads --colsep "\t" get_CUB_results :::: species_list.txt ::: $tmp_dir
cat $tmp_dir/*_cub.txt >> CUB_results.txt
rm -rf $tmp_dir

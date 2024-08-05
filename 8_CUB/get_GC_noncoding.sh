#!/bin/bash

# This script calculates the GC content of non-coding regions and the whole genome for a given species.
# It uses the geecee tool to calculate the GC content of the sequences.

# Usage: bash get_GC_noncoding.sh <species> <genome_dir>

# Check if input arguments are provided
species=$1
assembly=$2
genome_dir=$3

if [ -z "$species" ] || [ -z "$genome_dir" ]; then
    echo "Input arguments are missing. Please provide the species name and the genome directory."
    exit 1
fi

#sp=$(echo "$species" | awk '{print toupper(substr($0, 1, 1)) tolower(substr($0, 2))}')
sp=$(echo "$species" | awk '{print toupper($0)}')
# Find genome file in genome_dir
assembly_file=${genome_dir}/${assembly}

if [ ! -s "$assembly_file" ]; then
    echo "Genome file not found in the specified directory."
    exit 1
fi

gff3_file="/data/home/s2215768/fly_annotation/complement_annotations/clean_annotations/${sp}.gff"

tmp_dir=$(mktemp -d)

# Convert GFF3 to BED format
awk 'BEGIN{OFS="\t"} $3=="CDS" {print $1,$4-1,$5}' $gff3_file | sort -k1,1 -k2,2n >  $tmp_dir/${sp}_annotations.bed

rm -rf ${assembly_file}.fai
samtools faidx $assembly_file
cut -f1,2 ${assembly_file}.fai|sort -k1,1 -k2,2n > $tmp_dir/${sp}_genome.bed

# Get non-coding regions
bedtools complement -i $tmp_dir/${sp}_annotations.bed -g $tmp_dir/${sp}_genome.bed > $tmp_dir/${sp}_noncoding_regions.bed

# Extract non-coding sequences
bedtools getfasta -fi $assembly_file -bed $tmp_dir/${sp}_noncoding_regions.bed -fo $tmp_dir/${sp}_noncoding_sequences.fa

# Calculate GC content of non-coding sequences
geecee -sequence $tmp_dir/${sp}_noncoding_sequences.fa -outfile CUB_output/${sp}_gc_noncoding.tsv

# Calculate GC content of the whole genome
geecee -sequence $assembly_file -outfile CUB_output/${sp}_gc_whole_genome.tsv


# Clean up
rm -rf $tmp_dir

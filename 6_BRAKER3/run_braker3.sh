#!/bin/bash
# This script is used to run the BRAKER3 gene annotation pipeline for multiple genomes.
# It requires the following inputs:
# - A genome_info file containing the species name and their corresponding assemblies.
# - A protein sequence file named "Diptera_noDupl.faa".
# - RNA-seq data in the form of compressed fastq files.

# Note: This script assumes that the necessary software and data files are already set up and available in the specified directories.

source ~/miniconda3/etc/profile.d/conda.sh
conda activate braker3

genome_info=$1

if [ -z "$genome_info" ]; then
    echo "Error: genome_info not provided."
    exit 1
fi


if [ ! -s Diptera_noDupl.faa ]; then
    echo "Error: protein seq not provided."
    exit 1
fi


rna_seq(){
sp=$1
rna_dir="/data/home/s2215768/fly_annotation/data/rnaseq_all"

if [ -s ${rna_dir}/${sp}.tar.gz ] && [ ! -s bams/${sp}_1.fastq ];then
tar xzf ${rna_dir}/${sp}.tar.gz -C bams/
zcat bams/${sp}/${sp}_rnaseq_norm_R1.fastq.gz > bams/${sp}_1.fastq
zcat bams/${sp}/${sp}_rnaseq_norm_R2.fastq.gz > bams/${sp}_2.fastq
rm -rf bams/${sp}/
fi
}

export -f rna_seq
#parallel -j20 --colsep '\t' rna_seq {1} :::: $genome_info


run_braker(){
genome=$2
species=`echo $1|tr "[a-z]" "[A-Z]"`
genomes_dir="/data/home/s2215768/fly_annotation/data/CAT_genomes/selected_genomes"
sp=$1
mkdir -p $species

if [ ! -s ${species}/braker.gff3 ];then

if [ -s bams/${sp}_1.fastq ] && [ -s bams/${sp}_2.fastq ];then
echo "Run braker (For $species with RNAseq)"
braker.pl --species=fly --genome=${genomes_dir}/${genome} --prot_seq=Diptera_noDupl.faa --AUGUSTUS_ab_initio --threads=25 \
 --useexisting --gff3 --workingdir=$species --augustus_args="--species=fly" --rnaseq_sets_ids=$sp --rnaseq_sets_dir=bams/
cat "bams/${sp}_1.fastq" | pigz -p 25 > "bams/${sp}_1.fastq.gz"
cat "bams/${sp}_2.fastq" | pigz -p 25 > "bams/${sp}_2.fastq.gz"
rm -rf bams/${sp}_1.fastq bams/${sp}_2.fastq
echo "Braker run ended $species"

else
echo "Run braker (For $species without RNAseq)"
braker.pl --species=fly --genome=${genomes_dir}/${genome} --prot_seq=Diptera_noDupl.faa --AUGUSTUS_ab_initio --threads=25 \
 --useexisting --gff3 --workingdir=$species --augustus_args="--species=fly"
echo "Braker run ended $species"
fi

fi
}

export -f run_braker

parallel -j3 --colsep '\t' run_braker {1} {2} :::: $genome_info

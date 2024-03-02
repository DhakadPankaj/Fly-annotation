#! /bin/bash

# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate annotation

# Create directory for selected genomes
mkdir -p /data/home/s2215768/fly_annotation/data/CAT_genomes/selected_genomes

# Function to get species list
get_list(){

    sp=$1
    path=$3

    # Get taxon ID
    tax_id=`echo $sp|taxonkit name2taxid |cut -f2`

    # Find assembly file
    assembly=`find ${path}/$(echo ${sp}|tr " " "_").${2}*.rm.fna -type f`

    if [ -z ${assembly} ];then
        echo "$sp assembly not available"
    else
        # Copy assembly file to selected genomes directory
        spcy=`echo ${sp}|tr " " "_"`
        genome=`basename -a ${assembly}`
        echo -e "${spcy}\t${tax_id}\t${genome}" > tmp/${spcy}.tmp
    fi
}

export -f get_list

# Create temporary directory
mkdir -p tmp/

# Process nanopore files
for file in only_nanopore.tmp nanopore_better.tsv; do
    echo $file
    path="/data/home/s2215768/fly_annotation/data/CAT_genomes/masked_genomes"
    parallel -j 20 --colsep "\t" get_list :::: <(cut -f1 $file|grep -v "species") ::: "nanopore" ::: $path
done

# Process refseq files
for file in only_ncbi.tmp Selected_refseq.tsv; do
    echo $file
    path="/data/home/s2215768/fly_annotation/data/CAT_genomes/refseq_genomes/masked_genomes"
    parallel -j 20 --colsep "\t" get_list :::: <(cut -f1 $file|grep -v "species") ::: "GC" ::: $path
done

# Create species list file
>species_list.tsv
cat tmp/*.tmp >> species_list.tsv

# Copy species list file to selected genomes directory
cp species_list.tsv ../../CAT_genomes/selected_genomes/

# Remove temporary directory
rm -r tmp/

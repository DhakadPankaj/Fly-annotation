#!/bin/bash

# This script is used to complement the CAT annotation with the BRAKER annotation


tail -n +2 genome_info.tsv | awk -F'\t' -v OFS='\t' '{gsub(/ /,""); print $1, $3, $10}' | grep -Fwvf <(cut -f1 outgroup_species.tsv) > genome_info_selected.tsv

run_CATBraker(){
  species=`echo $1|tr [a-z] [A-Z]`
  genome=$2
  Refseq=$3
  genome_path="/data/home/s2215768/fly_annotation/data/CAT_genomes/selected_genomes"
  if [ ! -s /data/home/s2215768/fly_annotation/complement_annotations/longest_isoforms/${species}.fa ]; then
    set -x
    set -e
    cd $species
    rm -rf ${species}_longest_isoform.gff
    agat_sp_keep_longest_isoform.pl -gff ${species}.gff -o ${species}_longest_isoform.gff > ../logs/${species}_longest_isoform.log 2>&1
    rm -rf *report.txt *.log
    gffread -y - -C -H -V -J -g $genome_path/$genome ${species}_longest_isoform.gff | fasta_formatter | \
     sed "/^>/ s/mrna-/${species}_/" > ../longest_isoforms/${species}.fa
    cd ..
    set +x
  fi
}
mkdir -p longest_isoforms
export -f run_CATBraker
parallel -j 20 -k --colsep '\t' run_CATBraker {1} {2} {3} :::: genome_info_selected.tsv

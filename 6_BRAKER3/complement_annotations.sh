#!/bin/bash

# This script is used to complement the CAT annotation with the BRAKER annotation


tail -n +2 genome_info.tsv | awk -F'\t' -v OFS='\t' '{gsub(/ /,""); print $1, $3, $10}' | grep -Fwvf <(cut -f1 outgroup_species.tsv) > genome_info_selected.tsv

run_CATBraker(){
  species=`echo $1|tr [a-z] [A-Z]`
  genome=$2
  Refseq=$3
  if [ ! -s /data/home/s2215768/fly_annotation/complement_annotations/${species}/${species}.gff ]; then
    set -x
    echo "Running CATBraker for $species"
    bash /data/home/s2215768/fly_annotation/complement_annotations/CATBraker.sh $species $genome $Refseq > logs/${species}_CATBraker.log 2>&1
    set +x
  fi


}
mkdir -p logs
export -f run_CATBraker
parallel -j 20 -k --colsep '\t' run_CATBraker {1} {2} {3} :::: genome_info_selected.tsv


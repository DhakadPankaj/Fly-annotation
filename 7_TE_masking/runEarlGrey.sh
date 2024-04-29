#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh
conda activate earlgrey


genome_info=$1

if [ -z "$genome_info" ]; then
    echo "Error: genome_info not provided."
    exit 1
fi
repeatmasker(){
    genome=$2
    species=`echo $1|tr "[a-z]" "[A-Z]"`
    genomes_dir="/data/home/s2215768/fly_annotation/data/CAT_genomes/selected_genomes"
    if [ ! -s TE_annotation/${species}_EarlGrey/${species}_summaryFiles/${species}.filteredRepeats.gff ];then
    echo "Start": $(date)
    set -x
    mv eddie/${species}_EarlGrey TE_annotation/
    cd TE_annotation/${species}_EarlGrey/${species}_RepeatModeler
    rec_dir=$(find . -maxdepth 1 -type d -name "RM*" -print)

    count=$(echo "$rec_dir" | wc -l)

    if [ "$count" -gt 1 ]; then
       echo -e "Multiple RepeatModeler directories found. Keep only one recovery (RM*) directory in ${species}_EarlGrey/${species}_RepeatModeler/"
       exit 1
    else
       echo -e "RepeatModeler directory found: ${rec_dir}"
    fi
    RepeatModeler -engine ncbi -threads 15 -database ../${species}_Database/${species} -recoverDir ${rec_dir}/
    cd ../../
    earlGrey -g ${genomes_dir}/${genome} -s $species -o ./ -t 15
    set +x
    echo "End": $(date)
    else
    echo "${species} Repeatmasked genome already present... EXITING"
    fi
}

export -f repeatmasker
mkdir -p TE_annotation/
parallel -j 5 -k --colsep '\t' repeatmasker {1} {2} :::: $genome_info

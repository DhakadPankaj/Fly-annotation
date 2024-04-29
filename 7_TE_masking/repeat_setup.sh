#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh
conda activate earlgrey

mkdir -p Repeatmasked_genomes/


genome_info=$1

if [ -z "$genome_info" ]; then
    echo "Error: genome_info not provided."
    exit 1
fi
repeatmasker(){
    genome=$2
    species=`echo $1|tr "[a-z]" "[A-Z]"`
    genomes_dir="/data/home/s2215768/fly_annotation/data/CAT_genomes/selected_genomes"
    if [ ! -s Repeatmasked_genomes/${species}/masked/${genome}.out.gff ];then
    mkdir -p Repeatmasked_genomes/$species
    cd Repeatmasked_genomes/$species
    echo "Start": $(date)
    set -x
    BuildDatabase -name ${species}_db ${genomes_dir}/${genome}
    RepeatModeler -database ${species}_db -threads 22 -LTRStruct -quick 1>repeatmodeler.log 2>&1

    cd-hit-est -i ${species}_db-families.fa -o ${species}_reduced.fa -d 0 -aS 0.8 -c 0.9 -G 0 -g 1 -b 500

    # mask reference genome using above TE library
    RepeatMasker -pa 22 -gff -lib ${species}_reduced.fa -xsmall -dir masked ${genomes_dir}/${genome} 1>repeatmasker.log 2>&1
    cd ..
    set +x
    echo "End": $(date)
    else
    echo "${species} Repeatmasked genome already present... EXITING"
    fi
}

export -f repeatmasker
parallel -j 3 -k --colsep '\t' repeatmasker {1} {2} :::: <(head -25 $genome_info)


genome_size(){
    genome_path="/data/home/s2215768/fly_annotation/data/CAT_genomes/selected_genomes"
    size=`bioawk -c fastx '{ sum += length($seq) } END { print sum }' $genome_path/$2`
    echo -e "$1\t$2\t$3\t$size" >> genome_sizes.tmp
}
export -f genome_size
#> genome_sizes.tmp
#parallel -j 20 -k --colsep '\t' genome_size {1} {2} {3} :::: $1
#sort -k4,4nr genome_sizes.tmp > genome_sizes_sorted.tmp
#rm genome_sizes.tmp
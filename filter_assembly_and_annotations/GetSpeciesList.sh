#! /bin/bash

source ~/miniconda3/etc/profile.d/conda.sh
conda activate annotation

mkdir -p /data/home/s2215768/fly_annotation/data/CAT_genomes/selected_genomes

get_list(){

sp=$1
path=$3
tax_id=`echo $sp|taxonkit name2taxid |cut -f2`
assembly=`find ${path}/$(echo ${sp}|tr " " "_").${2}*.rm.fna -type f`

if [ -z ${assembly} ];then
echo "$sp assembly not available"

else
#cp $assembly /data/home/s2215768/fly_annotation/data/CAT_genomes/selected_genomes/
spcy=`echo ${sp}|tr " " "_"`
genome=`basename -a ${assembly}`

echo -e "${spcy}\t${tax_id}\t${genome}" > tmp/${spcy}.tmp

fi

}

export -f get_list
mkdir -p tmp/
for file in only_nanopore.tmp nanopore_better.tsv;do
echo $file
path="/data/home/s2215768/fly_annotation/data/CAT_genomes/masked_genomes"
parallel -j 20 --colsep "\t" get_list :::: <(cut -f1 $file|grep -v "species") ::: "nanopore" ::: $path
done

for file in only_ncbi.tmp Selected_refseq.tsv;do
echo $file
path="/data/home/s2215768/fly_annotation/data/CAT_genomes/refseq_genomes/masked_genomes"
parallel -j 20 --colsep "\t" get_list :::: <(cut -f1 $file|grep -v "species") ::: "GC" ::: $path
done

>species_list.tsv
cat tmp/*.tmp >> species_list.tsv
cp species_list.tsv ../../CAT_genomes/selected_genomes/
rm -r tmp/

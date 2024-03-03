#! /bin/bash

# Script Usage:
# bash gff_clean.sh
#
# Required Files:
# 1. sp_list_busco.tsv: This file contains information about the species and genome, and their BUSCO scores. It is used to generate the genome_info.tsv file.
#
# Usage Explanation:
# When you run this script using the bash command, it will perform the following steps:
# 1. Activate the conda environment named annotation.
# 2. Download annotated Refseq gff based on the information provided in the refseq_annotated.tsv file. The downloaded gff will be saved as *.genomic.gff.gz files.
# 3. Clean the downloaded GFF files by removing specific entries that are not required. The cleaned GFF files will be saved as *.cleaned.gff3 files.
# 4. Generate a genome_info.tsv file that contains information about the downloaded assemblies, including taxon IDs, assembly names, genome BUSCO scores, N50 values, and more.
#
# Please make sure to have the refseq_annotated.tsv and sp_list_busco.tsv files in the same directory as the script before running it.

# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate annotation

# Download annotated Refseq assemblies
# datasets download genome taxon 7214 --annotated --reference --exclude-seq --exclude-genomic-cds --exclude-protein --exclude-rna --dehydrated \
# --filename annotatedAssemblies.dataset.zip

# Define fields for annotation
field=$(echo "organism-name,assminfo-accession,annotinfo-name,annotinfo-release-date,annotinfo-source,\
annotinfo-featcount-gene-total,annotinfo-featcount-gene-protein-coding,annotinfo-busco-complete,annotinfo-busco-missing,assmstats-contig-n50,assminfo-level")

# Generate TSV file with selected annotations
dataformat tsv genome --fields $field --package annotatedAssemblies.dataset.zip|dos2unix |awk -F"\t" '($5=="FlyBase"||$5=="NCBI RefSeq")'|awk -F"\t" '$10>=1000000' > refseq_annotated.tsv

assembly_path="/data/home/s2215768/fly_annotation/data/CAT_genomes/selected_genomes"

# Function to download annotation
download_annotation() {
    if [[ $2 = GCF_* ]]; then
        esearch -db assembly -query $2 </dev/null | esummary| grep 'FtpPath type="RefSeq"'|sed 's/<FtpPath type="RefSeq">//'|sed 's_</FtpPath>__'|sed 's_\s__g' \
        | while read -r url ; do
            fname=$(echo $url | grep -o 'GCF_.*' | sed 's/$/_genomic.gff.gz/') ;
            wget -O "${1}.${2}_genomic.gff.gz"  "$url/$fname" ;
            gunzip ${1}.${2}_genomic.gff.gz
        done;
    else
        esearch -db assembly -query $2 </dev/null | esummary | xtract -pattern DocumentSummary -element FtpPath_GenBank \
        | while read -r url ; do
            fname=$(echo $url | grep -o 'GCA_.*' | sed 's/$/_genomic.gff.gz/') ;
            wget -O "${1}.${2}_genomic.gff.gz" "$url/$fname" ;
            gunzip ${1}.${2}_genomic.gff.gz
        done ;
    fi
}

# Function to clean gff file
gff_clean(){
    >${1}.out
    SP=`echo $1|tr "[a-z]" "[A-Z]"`
    cat $2 |grep -v "gene=mod(mdg4)" | grep -v "rna:NR_156799.1" | grep -v "rna-NR_133496.1" |grep -v "RNaseP:RNA" | grep -v "RNaseMRP:RNA" > ${SP}.tmp.gff
    gt gff3 -force -tidy  -sort  -retainids  -checkids  -o ${SP}.sorted.gff ${SP}.tmp.gff && rm ${SP}.tmp.gff
    grep -v "^##s" ${SP}.sorted.gff > ${SP}.tmp.gff
    rm ${SP}.sorted.gff
    convert_ncbi_gff3 ${SP}.tmp.gff ${SP}.tmp2.gff && rm ${SP}.tmp.gff
    grep -v "gene_biotype=tRNA" ${SP}.tmp2.gff > ${SP}.cleaned.gff3 && rm ${SP}.tmp2.gff
    gt gff3validator ${SP}.cleaned.gff3 >> ${1}.out
    validate_gff3 ${SP}.cleaned.gff3 >> ${1}.out
}

# Function to process reference annotations
ref_annotations(){
    sp=$1
    acc=$2
    assembly_path="/data/home/s2215768/fly_annotation/data/CAT_genomes/selected_genomes"

    if [ -s ${sp}.${acc}_genomic.gff ];then
        echo $sp gff file already downloaded
    else
        echo $sp gff downloading
        download_annotation $sp $acc
    fi
    SP=`echo $1|tr "[a-z]" "[A-Z]"`
    if [ -s ${SP}.cleaned.gff3 ];then
        echo $sp clean gff3 already exists
    else
        echo $sp gff3 cleaning
        gff_clean $sp ${sp}.${acc}_genomic.gff
    fi
}

# Export functions for parallel execution
export -f ref_annotations
export -f download_annotation
export -f gff_clean

# Process annotations in parallel
cat refseq_annotated.tsv|tr " " "_" > refseq_annotated.tmp
parallel -j 5  --colsep "\t" ref_annotations :::: refseq_annotated.tmp
rm refseq_annotated.tmp

#########################
### Create genome_info.tsv file (needed to generate cactus guide tree)
#########################

# Function to generate genome_info.tsv file
genome_info(){
    sp=$1
    SP=`echo $sp|tr "[a-z]" "[A-Z]"`
    assembly_path="/data/home/s2215768/fly_annotation/data/CAT_genomes/selected_genomes"
    assembly=$3
    id=$(echo $sp|tr "_" " "|taxonkit name2taxid |cut -f2)
    genome_busco=$4
    N50=`assembly-stats -s -u ${assembly_path}/${assembly}|cut -f9`
    if [ $(echo "$N50 >= 10000000"|bc) == 1 ];then
        cactus_name=`echo "*${SP}"`
        cactus_ref="Yes"
    else
        cactus_name="${SP}"
        cactus_ref="No"
    fi

    if [ -s ${SP}.cleaned.gff3 ];then
        annot_info=`grep "^${sp}$(printf "\t")" <(cat refseq_annotated.tsv|tr " " "_")`
        annot_busco=`grep "^${sp}$(printf "\t")" <(cat refseq_annotated.tsv|tr " " "_")|cut -f8`
        genes=`grep "^${sp}$(printf "\t")" <(cat refseq_annotated.tsv|tr " " "_")|cut -f7`
        annot="Yes"
    else
        annot_busco="NA"
        genes="NA"
        annot="No"
    fi
    bampath="/data/home/s2215768/fly_annotation/data/rnaseq_all/bams"
    if [ -s $bampath/${SP}_rnaseq.bam ];then
        RNA="Yes"
    else
        RNA="No"
    fi

    echo -e "${sp}\t${id}\t${assembly}\t${genome_busco}\t${N50}\t${cactus_ref}\t${cactus_name}\t${SP}\t${RNA}\t${annot}\t${annot_busco}\t${genes}" > tmp/${sp}.tmp
}

# Create temporary directory
mkdir -p tmp
export -f genome_info

# Process genome_info in parallel
parallel -j 30  --colsep "\t" genome_info :::: <(cat sp_list_busco.tsv|tr " " "_")
echo -e "species_name\taccession\tassembly\tgenome_busco\tN50\tcactus_ref\tcactus_name\thal_name\tRNA_seq\tAnnotated\tannot_busco\tprotein_coding_genes" > genome_info.tsv
cat tmp/*.tmp >> genome_info.tsv
rm -r tmp/

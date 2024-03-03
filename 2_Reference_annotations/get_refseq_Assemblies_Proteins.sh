# /bin/bash

# Script Usage:
# bash get_refseqprotein.sh
#
# Required Files:
# 1. refseq_annotated.tsv: This file contains information about the annotated Refseq assemblies.
# Usage Explanation:
# When you run this script using the bash command, it will perform the following steps:
# 1. Download protein sequences and genomic assemblies based on the information provided in the refseq_annotated.tsv file. The downloaded files will be saved as *.protein.faa.gz and *.genome.fna.gz files respectively.
# 2. The fasta_header_clean.sh script is run on the downloaded genomic assemblies to clean the fasta headers.
#
# Please make sure to have the refseq_annotated.tsv file in the same directory as the script before running it.


source ~/miniconda3/etc/profile.d/conda.sh
conda activate annotation


download_protein() {
    SP=`echo $1|tr "[a-z]" "[A-Z]"`
    if [[ $2 = GCF_* ]]
    then
    esearch -db assembly -query $2 </dev/null | esummary| grep 'FtpPath type="RefSeq"'|sed 's/<FtpPath type="RefSeq">//'|sed 's_</FtpPath>__'|sed 's_\s__g' \
     |while read -r url ; do
    fname=$(echo $url | grep -o 'GCF_.*' | sed 's/$/_protein.faa.gz/') ;
    wget -O "${SP}.${2}_protein.faa.gz"  "$url/$fname" ;
    gunzip ${SP}.${2}_protein.faa.gz
    done;
    else
    esearch -db assembly -query $2 </dev/null | esummary | xtract -pattern DocumentSummary -element FtpPath_GenBank \
     | while read -r url ; do
    fname=$(echo $url | grep -o 'GCA_.*' | sed 's/$/_protein.faa.gz/') ;
    wget -O "${SP}.${2}_protein.faa.gz" "$url/$fname" ;
    gunzip ${SP}.${2}_protein.faa.gz
    done ;
    fi
}

export -f download_protein
cat refseq_annotated.tsv|tr " " "_" > refseq_annotated.tmp
parallel -j 5  --colsep "\t" download_protein :::: refseq_annotated.tmp


download_assembly() {
    SP=`echo $1|tr "[a-z]" "[A-Z]"`
    if [[ $2 = GCF_* ]]
    then
    esearch -db assembly -query $2 </dev/null | esummary| grep 'FtpPath type="RefSeq"'|sed 's/<FtpPath type="RefSeq">//'|sed 's_</FtpPath>__'|sed 's_\s__g' \
     |while read -r url ; do
    fname=$(echo $url | grep -o 'GCF_.*' | sed 's/$/_genomic.fna.gz/') ;
    wget -O "${SP}.genome.fna.gz"  "$url/$fname" ;
    gunzip ${SP}.genome.fna.gz
    fasta_header_clean.sh ${SP}.genome.fna
    done;
    else
    esearch -db assembly -query $2 </dev/null | esummary | xtract -pattern DocumentSummary -element FtpPath_GenBank \
     | while read -r url ; do
    fname=$(echo $url | grep -o 'GCA_.*' | sed 's/$/_genomic.fna.gz/') ;
    wget -O "${SP}.genome.fna.gz" "$url/$fname" ;
    gunzip ${SP}.genome.fna.gz
    fasta_header_clean.sh ${SP}.genome.fna
    done ;
    fi
}

export -f download_assembly
cat refseq_annotated.tsv|tr " " "_" > refseq_annotated.tmp
parallel -j 5  --colsep "\t" download_assembly :::: refseq_annotated.tmp

rm refseq_annotated.tmp
#! /bin/bash

# Usage: ./genome_setup.sh
# This script performs genome masking using RepeatMasker for all *_genomic.fna files in the current directory.
# It requires the following dependencies:
# - conda (with the 'repeatmasker' environment activated)
# - gmes_linux_64_4 (probuild)
# - RepeatMasker

source ~/miniconda3/etc/profile.d/conda.sh
conda activate repeatmasker

threads=8

genome_masking() {
    fasta=$1

    PROBUILD="/data/home/s2215768/fly_annotation/tools/gmes_linux_64_4/probuild"
    RepeatMasker_path="/data/home/s2215768/fly_annotation/tools/RepeatMasker/RepeatMasker"

    # Fasta header cleaner
    cat ${fasta} | sed -E 's/(>[0-9A-Z\._]+).+/\1/' > $(echo ${fasta} | sed -E 's/(\.f.a)/\.renamed\1/')
    cat ${fasta} | grep ">" | sed -E 's/>([0-9A-Z\._]+)(.+)/\1\t\2/' > $(echo ${fasta} | sed -E 's/\.(f.a)/\.\1info/')

    #all uppercase
    $PROBUILD --reformat_fasta --in $fasta --out ${fasta%_genomic.fna}_UP.fna --uppercase 1 --letters_per_line 60 --original

    mkdir -p repeatmasker_data/${fasta%_genomic.fna}.rmdata

    $RepeatMasker_path --species drosophila_flies_genus -pa 10 -xsmall ${fasta%_genomic.fna}_UP.fna

    cp ${fasta%_genomic.fna}_UP.fna.masked repeatmasker_data/${fasta%_genomic.fna}.rmdata/
    mv ${fasta%_genomic.fna}_UP.fna.tbl repeatmasker_data/${fasta%_genomic.fna}.rmdata/
    mv ${fasta%_genomic.fna}_UP.fna.out repeatmasker_data/${fasta%_genomic.fna}.rmdata/
    mv ${fasta%_genomic.fna}_UP.fna.cat.gz repeatmasker_data/${fasta%_genomic.fna}.rmdata/

    mv ${fasta%_genomic.fna}_UP.fna.masked masked_genomes/${fasta%_genomic.fna}.rm.fna

    cd repeatmasker_data/
    tar cf - ${fasta%_genomic.fna}.rmdata | pigz -p 8 > ${fasta%_genomic.fna}.rmdata.tar.gz
    cd ..

    rm ${fasta%_genomic.fna}_UP.fna
}

main() {
    mkdir -p repeatmasker_data/
    mkdir -p masked_genomes/

    export -f genome_masking

    if [ -z "$(ls -1 *_genomic.fna 2>/dev/null)" ]; then
        echo "No *_genomic.fna files found."
        exit 1
    else
        parallel -j $threads genome_masking {} ::: $(ls -1 *_genomic.fna)
    fi
}

main

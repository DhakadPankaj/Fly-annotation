#!/bin/bash

# This script is used to download and compare genome assemblies from NCBI for a specific species.

# Function to download genome from NCBI
# Usage: download_assembly_ftp species_name assembly_acc
download_assembly_ftp() {
    refseq_dir=$3
    if [[ $2 = GCF_* ]]; then
        esearch -db assembly -query $2 </dev/null | esummary | grep 'FtpPath type="RefSeq"' | sed 's/<FtpPath type="RefSeq">//' | sed 's_</FtpPath>__' | sed 's_\s__g' \
        | while read -r url ; do
            fname=$(echo $url | grep -o 'GCF_.*' | sed 's/$/_genomic.fna.gz/') ;
            wget -O "${1}.${2}_genomic.fna.gz"  "$url/$fname" ;
            gunzip ${1}.${2}_genomic.fna.gz && mv ${1}.${2}_genomic.fna ${refseq_dir}/;
        done;
    else
        esearch -db assembly -query $2 </dev/null | esummary | xtract -pattern DocumentSummary -element FtpPath_GenBank \
        | while read -r url ; do
            fname=$(echo $url | grep -o 'GCA_.*' | sed 's/$/_genomic.fna.gz/') ;
            wget -O "${1}.${2}_genomic.fna.gz" "$url/$fname" ;
            gunzip ${1}.${2}_genomic.fna.gz && mv ${1}.${2}_genomic.fna ${refseq_dir}/;
        done ;
    fi
}

# Create directories
mkdir -p ~/fly_annotation/data/CAT_genomes
nanopore_dir="/data/home/s2215768/fly_annotation/data/CAT_genomes/masked_genomes"
refseq_dir="/data/home/s2215768/fly_annotation/data/CAT_genomes/refseq_genomes"

# Define fields for genome metadata
fields="organism-name,tax-id,assminfo-accession,assminfo-name,assminfo-biosample-accession,assminfo-level,assminfo-sequencing-tech,assmstats-total-number-of-chromosomes,assmstats-number-of-contigs,assmstats-number-of-scaffolds,assmstats-contig-n50,assmstats-scaffold-n50,annotinfo-busco-complete,annotinfo-busco-duplicated,annotinfo-busco-fragmented,annotinfo-busco-missing,annotinfo-busco-totalcount,annotinfo-featcount-gene-total"

# Download reference assemblies metadata
datasets download genome taxon 7214 --reference --exclude-seq --exclude-genomic-cds --exclude-protein --exclude-rna --dehydrated --filename refonly.dataset.zip

# Download all assemblies metadata
datasets download genome taxon 7214 --exclude-seq --exclude-genomic-cds --exclude-protein --exclude-rna --dehydrated --filename allAssemblies.dataset.zip

# Create nanopore genome table
if [[ -s nanopore_genomes.tsv ]]; then
    echo "Nanopore genome table already present. Skipping..."
else
    >nanopore_genomes.tsv
    mkdir -p tmp
    export -f get_stats
    parallel -j 50 get_stats ::: $(basename -a ${nanopore_dir}/*.nanopore.rm.fna)
    cat tmp/*tsv >> nanopore_genomes.tsv
    rm -r tmp/
fi

# Species with only nanopore genomes
grep -Fxvf <(cut -f1,2 all.tmp|sort -u) <(cut -f1,2 nanopore_genomes.tsv|sort) > only_nanopore.tmp

# Species with only NCBI genomes
grep -Fxvf <(cut -f1,2 nanopore_genomes.tsv|sort) <(cut -f1,2 all.tmp|sort -u) > only_ncbi.tmp

# Species with both NCBI & nanopore genomes
grep -Fxf <(cut -f1,2 all.tmp|sort -u) <(cut -f1,2 nanopore_genomes.tsv|sort) > both.tmp

# Download Refseq genomes
download_only_refseq() {
    acc=$(grep -w "${1}$(printf '\t')${2}" ref_assemblies.tmp | cut -f3)
    sp=$(echo $1 | tr " " "_")
    N50=$(grep -w "${1}$(printf '\t')${2}" ref_assemblies.tmp | cut -f12)

    if [ $N50 -lt 50000 ]; then
        rm -f ${3}/${sp}.${acc}_genomic.fna
        echo "$sp $N50" >> refseq_N50_lt50k.txt
    else
        download_assembly_ftp $sp $acc $3
    fi
}

export -f download_only_refseq
export -f download_assembly_ftp
parallel -j 10 --colsep "\t" download_only_refseq :::: <(grep -Fxf <(cut -f1,2 ref_assemblies.tmp) <(cat only_ncbi.tmp )) ::: ${refseq_dir}

# Download Genebank genomes
download_only_genebank() {
    acc=$(grep -w "${1}$(printf '\t')${2}" all.tmp | sort -t$'\t' -k6,6 -k12,12nr -k9,9n | head -1 | cut -f3)
    sp=$(echo $1 | tr " " "_")
    N50=$(grep -w "${1}$(printf '\t')${2}" all.tmp | sort -t$'\t' -k6,6 -k12,12nr -k9,9n | head -1 | cut -f12)

    if [ $N50 -lt 50000 ]; then
        rm -f ${3}/${sp}.${acc}_genomic.fna
        echo "$sp $N50" >> genebank_N50_lt50k.txt
    else
        download_assembly_ftp $sp $acc $3
    fi
}

export -f download_only_genebank
export -f download_assembly_ftp
parallel -j 10 --colsep "\t" download_only_genebank :::: <(grep -Fxvf <(grep -Fxf <(cut -f1,2 ref_assemblies.tmp) <(cat only_ncbi.tmp )) only_ncbi.tmp) ::: ${refseq_dir}

# Download NCBI genomes and compare with nanopore
export -f download_only_genebank
export -f download_assembly_ftp
export -f download_only_refseq
mkdir -p compare_dir

parallel -j 10 --colsep "\t" download_only_refseq :::: <(grep -Fxf <(cut -f1,2 ref_assemblies.tmp) both.tmp) ::: "compare_dir"
parallel -j 5 --colsep "\t" download_only_genebank :::: <(grep -Fxvf <(cut -f1,2 ref_assemblies.tmp) both.tmp) ::: "compare_dir"

# Compare with nanopore
compare() {
    nanopore_dir="/data/home/s2215768/fly_annotation/data/CAT_genomes/masked_genomes"
    refseq_dir="/data/home/s2215768/fly_annotation/data/CAT_genomes/refseq_genomes"
    sp=$(echo $1 | tr " " "_")
    ncbi=$(find compare_dir/${sp}.GC*_genomic.fna -type f)
    nanopore=$(find ${nanopore_dir}/${sp}.nanopore.rm.fna -type f)

    if [ -z $ncbi ]; then
        echo "$sp NCBI genome: N50 lt50k"
        echo -e "${1}\t${nanopre}" >> Selected_genomes.tsv
    else
        num=$(bioawk -c fastx '{ print length($seq) }' $ncbi | sort -k1nr | awk '$1<1000000' | paste -sd+ - | bc)
        tot=$(bioawk -c fastx '{ print length($seq) }' $ncbi | sort -k1nr | paste -sd+ - | bc)

        ncbi_frag=$(awk -vn=$num -vt=$tot 'BEGIN{printf("%.2f\n",n*100/t)}')
        ncbi_N50=$(assembly-stats -s -u $ncbi | cut -f9)
        ncbi_contig=$(assembly-stats -s -u $ncbi | cut -f3)

        num1=$(bioawk -c fastx '{ print length($seq) }' $nanopore | sort -k1nr | awk '$1<1000000' | paste -sd+ - | bc)
        tot1=$(bioawk -c fastx '{ print length($seq) }' $nanopore | sort -k1nr | paste -sd+ - | bc)

        nano_frag=$(awk -vn=$num1 -vt=$tot1 'BEGIN{printf("%.2f\n",n*100/t)}')
        nano_N50=$(assembly-stats -s -u $nanopore | cut -f9)
        nano_contig=$(assembly-stats -s -u $nanopore | cut -f3)

        if [ $(echo "$ncbi_frag $nano_frag" | awk '{print ($1 <= $2)}') == 1 ] || [ $(echo "$ncbi_contig $nano_contig" | awk '{print ($1 <= $2)}') == 1 ]; then
            echo -e "${sp}\t${ncbi_frag}\t${nano_frag}\t${ncbi_contig}\t${nano_contig}\t$(basename -a $ncbi)\t${ncbi_N50}\t${nano_N50}" >> Selected_refseq.tsv
            echo -e "${1}\t${ncbi}" >> Selected_genomes.tsv
        else
            echo "For ${sp} Nanopore is better"
            echo -e "${sp}\t${ncbi_frag}\t${nano_frag}\t${ncbi_contig}\t${nano_contig}\t$(basename -a $nanopore)\t${ncbi_N50}\t${nano_N50}" >> nanopore_better.tsv
            echo -e "${1}\t${nanopre}" >> Selected_genomes.tsv
        fi
    fi
}

export -f compare
>nanopore_better.tsv
>Selected_refseq.tsv
>Selected_genomes.tsv
parallel -j 10 --colsep "\t" compare :::: both.tmp

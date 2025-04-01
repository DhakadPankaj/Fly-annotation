#! /bin/bash


# This script is used to get the sequences from aligned HOGs for each species.
mkdir -p CDS_trimmed
get_all_seq(){
    align_dir="/data/home/s2215768/fly_annotation/complement_annotations/Results_MSA/og_alignments"
    HOG=$1
    # Check if alignment file exists
    if [ -s ${align_dir}/${HOG}/${HOG}_final_align_NT.aln ]; then
        awk '/^>/ {print; next} {gsub("-", ""); print}' ${align_dir}/${HOG}/${HOG}_final_align_NT.aln > CDS_trimmed/${HOG}_final_align_NT.tmp
    fi
}

if [ ! -s all_HOGs_NT.aln ]; then
export -f get_all_seq
parallel -j 50 get_all_seq ::: $(tail -n +2 ~/fly_annotation/complement_annotations/Results_MSA/HOG_filtered.tsv|cut -f 1)
cat CDS_trimmed/*_final_align_NT.tmp > all_HOGs_NT.aln
fi

get_seq_sp(){
    SP=$1
    col=$(head -1 ~/fly_annotation/complement_annotations/Results_MSA/HOG_filtered.tsv|tr "\t" "\n"|grep -n -w "${SP}"|cut -d ":" -f1)
    cut -f1,$col ../Results_MSA/HOG_filtered.tsv | awk -F'\t' '$2 != ""'|tail -n +2| awk -F'\t' '{split($2, genes, ", ");for (i in genes) {print $1, genes[i];}}' OFS='\t' > tmp/${SP}_genes_HOGs.txt

    # get sequences
    seqtk subseq all_HOGs_NT.aln <(cut -f2 tmp/${SP}_genes_HOGs.txt) > CDS_trimmed/${SP}.tmp.fa 

    # get genes with FPKM values from "HOGs_categories_20.txt"
    join -t $'\t' -1 1 -2 1 -o 2.3,1.2,2.2 <(sort -k1 tmp/${SP}_genes_HOGs.txt) <(sort -k1 HOGs_categories_20.txt)|sort -k1n|grep -f <(grep ">" CDS_trimmed/${SP}.tmp.fa|sed 's/^>//') > gene_expression/${SP}_genes_HOGs_FPKM.txt

    # Reorder sequences according to fpkm
    samtools faidx CDS_trimmed/${SP}.tmp.fa 
    samtools faidx CDS_trimmed/${SP}.tmp.fa $(cut -f2 gene_expression/${SP}_genes_HOGs_FPKM.txt) > CDS_trimmed/${SP}.fa

    # Add header in genes_HOGs_FPKM.txt
    sed -i '1i Category\tGene\tExpression' gene_expression/${SP}_genes_HOGs_FPKM.txt

    rm CDS_trimmed/${SP}.tmp.fa*
    rm tmp/${SP}_genes_HOGs.txt
    
    # get highly expressed genes
    #grep ">" Highly_exp_genes/${SP}.heg.fa|sed 's/>//g' > tmp/${SP}_highly_expressed.genes
    #seqtk subseq all_HOGs_NT.aln tmp/${SP}_highly_expressed.genes > HEG_genes/${SP}_highly_expressed.fa

}

export -f get_seq_sp
mkdir -p tmp
mkdir -p HEG_genes
mkdir -p gene_expression
#seqkit grep -v -n -f <(grep ">" ../HEG_genes/CACOXENUS_INDAGATOR_highly_expressed.fa|sed 's/^>//') ../CDS_trimmed/CACOXENUS_INDAGATOR.fa


parallel -j 50 get_seq_sp ::: $(cut -f2- ~/fly_annotation/complement_annotations/Results_MSA/HOG_filtered.tsv|head -1|tr "\t" "\n")
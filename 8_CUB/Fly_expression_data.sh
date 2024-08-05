#!/bin/bash

# This script is used to extract the top 1000 genes with the highest expression levels from the gene_rpkm_matrix_fb_2024_01.tsv file.
#sconda agat
threads=10

if [ ! -s gene_rpkm_matrix_fb_2024_01.tsv ]; then
    echo "gene_rpkm_matrix_fb_2024_01.tsv file not found."
    exit 1
    else
    # Get the top 1000 genes with the highest expression levels for adult whole body
    # RNA-Seq_Profile_FlyAtlas2_Adult_Female_Whole_(FBlc0003634): 144      RNA-Seq_Profile_FlyAtlas2_Adult_Male_Whole_(FBlc0003649):159
    awk 'BEGIN{FS=OFS="\t"} $4=="protein_coding_gene" { sum = $144 + $159; print $1, $2, sum }' gene_rpkm_matrix_fb_2024_01.tsv|\
    sort -k3,3nr | head -n 1000 > top_1000_genes.txt
    #awk 'BEGIN{FS=OFS="\t"} $4=="protein_coding_gene" { sum = 0; for (i = 5; i <= NF; i++) sum += $i; print $1, $2, sum }' gene_rpkm_matrix_fb_2024_01.tsv | \
  #sort -k3,3nr | head -n 1000 > top_1000_genes.txt
fi


if [ ! -s NP_FBid.tsv ]; then
    agat_sp_extract_attributes.pl -m --gff \
    ../clean_annotations/DROSOPHILA_MELANOGASTER.gff --att ID,Parent,dbxref,gene_name -p mRNA \
    -o DROSOPHILA_MELANOGASTER.gffattributes.tsv >> /dev/null 2>&1

    cat DROSOPHILA_MELANOGASTER.gffattributes.tsv |sort -u |\ 
    awk '{FS="\t"; OFS="\t"} { if(match($3, /FLYBASE:FBgn[0-9]+/)) { print $1,$2,substr($3, RSTART, RLENGTH),$4 } }'|sed 's_FLYBASE:__' > NP_FBid.tsv
fi


#join -t $'\t' -1 3 -2 1 -o 1.1,2.1,2.2,2.3 <(sort -k3 NP_FBid.tsv) <(sort -k1 top_1000_genes.txt)|sort -k4,4nr > top_1000_genes_with_FBid.txt



# Get the HOG id associated with each gene
get_hogid(){
    NP=$1
    FB=$2
    gene=$3
    rpkm=$4
    tmp_dir=$5
    if [ `grep -w "${NP}" ../Results_MSA/HOG_filtered.tsv|cut -f1 || [[ $? == 1 ]]` ];then
        OG=`grep -w "${NP}" ../Results_MSA/HOG_filtered.tsv|cut -f1`
        echo -e "${OG}\t${rpkm}\t${gene}\t${FB}\t${NP}" > "${tmp_dir}/${gene}_${FB}.tmp"
    fi
}

export -f get_hogid
if [ ! -s top_1000_genes_with_FBid_HOGid.txt ]; then
parallel -j $threads --colsep '\t' get_hogid {1} {2} {3} {4} {5} :::: $PWD/top_1000_genes_with_FBid.txt ::: $tmp_dir
cat $tmp_dir/*.tmp |sort -k2,2nr > top_1000_genes_with_FBid_HOGid.txt
rm -rf $tmp_dir/
fi


grep -Fwf <(awk '$5>=300' ../Results_MSA/Anc_OG.tsv |cut -f2) top_1000_genes_with_FBid_HOGid.txt |head -100 > top_100_genes_with_FBid_HOGid.txt

mkdir -p Highly_exp_genes
get_100_genes(){
    sp=$1
    tmp_dir=$(mktemp -d)
    col=$(head -1 ../Results_MSA/HOG_filtered.tsv|tr '\t' '\n'|grep -n -w "${sp}"|cut -d ":" -f1)
    grep -Fwf <(cut -f1 top_100_genes_with_FBid_HOGid.txt) ../Results_MSA/HOG_filtered.tsv |cut -f$col|sed 's/, /\n/g'|sed '/^$/d' > ${tmp_dir}/genes.tmp
    seqtk subseq CDS/${sp}.fa ${tmp_dir}/genes.tmp |seqtk seq -L 241 > Highly_exp_genes/${sp}.heg.fa
    rm -rf ${tmp_dir}
}
export -f get_100_genes
parallel -j $threads --colsep '\t' get_100_genes :::: <(head -1 ../Results_MSA/HOG_filtered.tsv|cut -f2-|tr '\t' '\n')

# Get the top 50 HOGs with the highest expression levels for each species
: <<'END_COMMENT'
tmp_dir=$(mktemp -d)
# Top ~50 HOG with the highest expression for each species
species=$1
col=$(head -1 HOGs.tsv |awk -F'\t' -v sp="$species" '{ for (i=1; i<=NF; i++) { if ($i == sp) { print i; break; } } }')

awk -F'\t' -v sp=$col '$sp != ""' HOGs.tsv|cut -f1,$col > $tmp_dir/HOGs_${species}.tsv

awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$0;next} $1 in a{print a[$1],$0}' top_1000_genes_with_FBid_HOGid.txt $tmp_dir/HOGs_${species}.tsv | \
sort -k2,2nr | cut -f1-5,7 | awk -F'\t' '{split($NF, arr, ", "); print $1, $2, $3, $4, $5, arr[int(rand()*length(arr))+1]}' 
END_COMMENT

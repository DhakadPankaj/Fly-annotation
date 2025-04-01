# This script is used to categorize the genes based on their expression levels.
# Dmel genes are first categorised first and then the orthologous genes in the other species are categorised based on the Dmel categories.

# Get FPKM values for Dmel genes
awk 'BEGIN{FS=OFS="\t"} $4=="protein_coding_gene" {
    sum = 0;
    for (i = 5; i <= NF; i++) {
        sum += $i;
    }
    print $1, $2, sum
}' gene_rpkm_matrix_fb_2024_01.tsv | sort -k3,3nr > Dmel_genes_FPKM.txt

# get HOGs and Dmel genes
col=$(head -1 ../Results_MSA/HOG_filtered.tsv | tr '\t' '\n' | grep -n -w "DROSOPHILA_MELANOGASTER" | cut -d ":" -f1)

cut -f1,$col ../Results_MSA/HOG_filtered.tsv | awk -F'\t' '$2 != ""'|tail -n +2| awk -F'\t' '{
    split($2, genes, ", ");
    for (i in genes) {
        print $1, genes[i];
    }
}' OFS='\t' > Dmel_genes_HOGs.txt

if [ ! -s NP_FBid.tsv ]; then
    agat_sp_extract_attributes.pl -m --gff \
    ../clean_annotations/DROSOPHILA_MELANOGASTER.gff --att ID,Parent,dbxref,gene_name -p mRNA \
    -o DROSOPHILA_MELANOGASTER.gffattributes.tsv >> /dev/null 2>&1

    cat DROSOPHILA_MELANOGASTER.gffattributes.tsv |sort -u |\ 
    awk '{FS="\t"; OFS="\t"} { if(match($3, /FLYBASE:FBgn[0-9]+/)) { print $1,$2,substr($3, RSTART, RLENGTH),$4 } }'|sed 's_FLYBASE:__' > NP_FBid.tsv
fi

join -t $'\t' -1 3 -2 1 -o 1.1,2.1,2.2,2.3 <(sort -k3 NP_FBid.tsv) <(sort -k1 Dmel_genes_FPKM.txt) > Dmel_genes_FBid_FPKM.txt

# Get the HOG id associated with each gene
join -t $'\t' -1 2 -2 1 -o 1.1,2.4 <(sort -k2 Dmel_genes_HOGs.txt) <(sort -k1 Dmel_genes_FBid_FPKM.txt)|sort -k1 |awk 'BEGIN{FS=OFS="\t"} {count[$1]++; sum[$1] += $2;} END {for (id in count) {print id, sum[id] / count[id];}}' |sort -k2n > HOGs_categories.txt
 
# Assign categories of HOGs based on their FPKM values (40 categories)
total=$(wc -l HOGs_categories.txt | cut -d " " -f1)
cat HOGs_categories.txt | awk -v total=$total 'BEGIN{FS=OFS="\t"} {print $1, $2, int((NR - 1) / total * 40) + 1;}' > HOGs_categories_40.txt




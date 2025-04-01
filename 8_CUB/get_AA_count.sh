#! /bin/bash

head -1 ../Results_MSA/HOG_filtered.tsv > HOGs.txt
awk -v FS='\t' -v OFS='\t' -v num1=305 -v num2=30 '{non_empty=0; for(i=2; i<=NF; i++) if($i >= "1") non_empty++; if(non_empty < num1 && non_empty >= num2) print $0}' <(tail -n +2 ../Results_MSA/HOG_stat_filtered.tsv) | cut -f1 > hog_mt29sp.txt
grep -Fwf hog_mt29sp.txt ../Results_MSA/HOG_filtered.tsv >> HOGs.txt

# extract genes from single copy HOGs

get_seq() {
    species=$(echo $1|tr [a-z] [A-Z] )
    if [ ! -s HOGs_seq/${species}.fa ]; then
    col=$(grep -w -n "${species}" <( head -1 ../Results_MSA/HOG_filtered.tsv|tr "\t" "\n")|cut -d: -f1)
    echo $species $col
    tail -n +2 HOGs.txt | cut -f$col| sed '/^$/d'| tr ", " "\n"|sed '/^$/d' > tmp/${species}.tmp
    seqtk subseq ../Results_MSA/All_CDS.fasta tmp/${species}.tmp > HOGs_seq/${species}.fa
    fi
    rm -rf tmp/${species}.tmp
    cusp -sequence HOGs_seq/${species}.fa -outfile rscu/${species}.rscu

    #get the total number of Amino acids
    total_AA=$(sed -e '/^#/d' -e '/^$/d' rscu/${species}.rscu |awk '{print $2"\t"$5}'|awk '{sum[$1] += $2} END {for (i in sum) print i"\t"sum[i]}' |\
    sort -k1,1|grep -w -v "*"|cut -f2|sed '/^$/d'|awk '{sum += $1} END {print sum}')
    echo "$1 $total_AA"
    #get the proportion of each Amino acids
    AA_prop=$(sed -e '/^#/d' -e '/^$/d' rscu/${species}.rscu |awk '{print $2"\t"$5}'|awk '{sum[$1] += $2} END {for (i in sum) print i"\t"sum[i]}' |\
    sort -k1,1|grep -w -v "*"|cut -f2|sed '/^$/d'|awk -v total=$total_AA '{for (i = 1; i <= NF; i++) printf "%.6f\t", $i / total}')
    # store in AA Prop table
    echo -e "${1}\t${AA_prop}" > tmp/${species}_AA_count.txt
    
}

export -f get_seq

mkdir -p tmp
mkdir -p HOGs_seq
mkdir -p rscu

echo -e "Species\tA\tC\tD\tE\tF\tG\tH\tI\tK\tL\tM\tN\tP\tQ\tR\tS\tT\tV\tW\tY" > AA_count_mt29spp.txt
parallel -j 30 --colsep "\t" get_seq {1} :::: ../genome_info_selected.tsv

cat tmp/*_AA_count.txt >> AA_count_mt29spp.txt
#rm -rf tmp


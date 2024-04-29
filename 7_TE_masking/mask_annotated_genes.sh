#!/bin/bash

mask_genes(){
    species=`echo $1|tr [a-z] [A-Z]`
    genome=$2
    genome_path="/data/home/s2215768/fly_annotation/data/CAT_genomes/selected_genomes"
    gff="/data/home/s2215768/fly_annotation/complement_annotations/${species}/${species}_longest_isoform.gff"
    blastdb="/data/home/s2215768/fly_annotation/TE_annotations/filter_TEs/dmel-all-CDS-r6.56.fasta"
    if [ -s TE_annotation/${species}_EarlGrey/${species}_summaryFiles/${species}.filteredRepeats.gff ];then
    repeats="/data/home/s2215768/fly_annotation/TE_annotations/TE_annotation/${species}_EarlGrey/${species}_summaryFiles/${species}.filteredRepeats.gff"
    strained="/data/home/s2215768/fly_annotation/TE_annotations/TE_annotation/${species}_EarlGrey/${species}_summaryFiles/${species}-families.fa.strained"
    echo "Masking $species genes with $repeats"
    mkdir -p temp4/$species
    cd temp4/$species

    grep -P -v "^>.*#Unknown" $strained|grep "^>"|sed 's/>//' > known_TEs.txt
    seqtk subseq $strained known_TEs.txt|awk '/^>/ {sub(/#.*/, ""); print toupper($0); next} 1' > known_TEs.fa
    # Remove TEs that are genuine protein hits.
    blastn -task blastn -query known_TEs.fa -db $blastdb \
    -outfmt '6 qseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
    -max_target_seqs 25 -culling_limit 2 -num_threads 4 -evalue 1e-2 -out ${species}_blast.out

    cut -f1 ${species}_blast.out|sort -u|sed 's/^/ID=/' > genuine_proteins.txt

# Exon-repeat overlaps filtering criterias
#: <<'END'
    bedtools intersect -a <(awk '$3=="exon"' $gff) -b <(awk '$3!="Unknown"' $repeats|grep -v "SHORTTE=T"|grep -Fwvf genuine_proteins.txt) -wao|cut -f1-5,9,19|tr ';' '\t'|\
    awk -v FS="\t" -v OFS="\t" '{a[$1 FS $2 FS $3 FS $4 FS $5 FS $6 FS $7]+=$8} END{for (i in a) print i, a[i]}'| \
    sort -k6,6|awk -v FS="\t" -v OFS="\t" '{print $5-$4+1,$7,$8}'| \
    awk -v FS="\t" -v OFS="\t" '{a[$2]+=$1; b[$2]+=$3} END{for (i in a) print i, a[i], b[i]}'| \
    awk -v FS="\t" -v OFS="\t" '{print $1,$3/$2}'|awk '$2>0.9'|cut -f1|sed 's/Parent=//' > kill_list_0.9.txt

    bedtools intersect -a <(awk '$3=="exon"' $gff) -b <(awk '$3!="Unknown"' $repeats|grep -v "SHORTTE=T"|grep -Fwvf genuine_proteins.txt) -wao|cut -f1-5,9,19|tr ';' '\t'|\
    awk -v FS="\t" -v OFS="\t" '{a[$1 FS $2 FS $3 FS $4 FS $5 FS $6 FS $7]+=$8} END{for (i in a) print i, a[i]}'| \
    sort -k6,6|awk -v FS="\t" -v OFS="\t" '{print $5-$4+1,$7,$8}'| \
    awk -v FS="\t" -v OFS="\t" '{a[$2]+=$1; b[$2]+=$3} END{for (i in a) print i, a[i], b[i]}'| \
    awk -v FS="\t" -v OFS="\t" '{print $1,$3/$2}'|awk '$2>0.95'|cut -f1|sed 's/Parent=//' > kill_list_0.95.txt

    bedtools intersect -a <(awk '$3=="exon"' $gff) -b <(awk '$3!="Unknown"' $repeats|grep -v "SHORTTE=T"|grep -Fwvf genuine_proteins.txt) -wao|cut -f1-5,9,19|tr ';' '\t'|\
    awk -v FS="\t" -v OFS="\t" '{a[$1 FS $2 FS $3 FS $4 FS $5 FS $6 FS $7]+=$8} END{for (i in a) print i, a[i]}'| \
    sort -k6,6|awk -v FS="\t" -v OFS="\t" '{print $5-$4+1,$7,$8}'| \
    awk -v FS="\t" -v OFS="\t" '{a[$2]+=$1; b[$2]+=$3} END{for (i in a) print i, a[i], b[i]}'| \
    awk -v FS="\t" -v OFS="\t" '{print $1,$3/$2}'|awk '$2>0.5'|cut -f1|sed 's/Parent=//' > kill_list_0.5.txt

    bedtools intersect -a <(awk '$3=="exon"' $gff) -b <(awk '$3!="Unknown"' $repeats|grep -v "SHORTTE=T"|grep -Fwvf genuine_proteins.txt) -wao|cut -f1-5,9,19|tr ';' '\t'|\
    awk -v FS="\t" -v OFS="\t" '{a[$1 FS $2 FS $3 FS $4 FS $5 FS $6 FS $7]+=$8} END{for (i in a) print i, a[i]}'| \
    sort -k6,6|awk -v FS="\t" -v OFS="\t" '{print $5-$4+1,$7,$8}'| \
    awk -v FS="\t" -v OFS="\t" '{a[$2]+=$1; b[$2]+=$3} END{for (i in a) print i, a[i], b[i]}'| \
    awk -v FS="\t" -v OFS="\t" '{print $1,$3/$2}'|awk '$2>0.3'|cut -f1|sed 's/Parent=//' > kill_list_0.3.txt
#END

# gene-repeat overlaps filtering criterias
: <<'END'
    bedtools intersect -a <(awk '$3=="gene"' $gff) -b <(awk '$3!="Unknown"' $repeats|grep -v "SHORTTE=T"|grep -Fwvf genuine_proteins.txt) -wa -f 0.92 | sort -u | grep -o 'ID=[^;]*' | cut -d= -f2 > kill_list_0.92.txt
    bedtools intersect -a <(awk '$3=="gene"' $gff) -b <(awk '$3!="Unknown"' $repeats|grep -v "SHORTTE=T"|grep -Fwvf genuine_proteins.txt) -wa -f 0.98 | sort -u | grep -o 'ID=[^;]*' | cut -d= -f2 > kill_list_0.98.txt
    bedtools intersect -a <(awk '$3=="gene"' $gff) -b <(awk '$3!="Unknown"' $repeats|grep -v "SHORTTE=T"|grep -Fwvf genuine_proteins.txt) -wa -f 0.95 | sort -u | grep -o 'ID=[^;]*' | cut -d= -f2 > kill_list_0.95.txt
    bedtools intersect -a <(awk '$3=="gene"' $gff) -b <(awk '$3!="Unknown"' $repeats|grep -v "SHORTTE=T"|grep -Fwvf genuine_proteins.txt) -wa -f 1 | sort -u | grep -o 'ID=[^;]*' | cut -d= -f2 > kill_list_1.txt
END
    agat_sp_filter_feature_from_kill_list.pl --gff $gff --kill_list kill_list_0.3.txt --output ${species}_clean_0.3.gff > log.tmp
    agat_sp_filter_feature_from_kill_list.pl --gff $gff --kill_list kill_list_0.95.txt --output ${species}_clean_0.95.gff > log.tmp1
    agat_sp_filter_feature_from_kill_list.pl --gff $gff --kill_list kill_list_0.5.txt --output ${species}_clean_0.5.gff > log.tmp2
    agat_sp_filter_feature_from_kill_list.pl --gff $gff --kill_list kill_list_0.9.txt --output ${species}_clean_0.9.gff > log.tmp3

    gffread -x - -C -H -V -J -g $genome_path/$genome ${species}_clean_0.3.gff | fasta_formatter > ${species}_mRNA_0.3.fasta
    gffread -x - -C -H -V -J -g $genome_path/$genome ${species}_clean_0.5.gff | fasta_formatter > ${species}_mRNA_0.5.fasta
    gffread -x - -C -H -V -J -g $genome_path/$genome ${species}_clean_0.95.gff | fasta_formatter > ${species}_mRNA_0.95.fasta
    gffread -x - -C -H -V -J -g $genome_path/$genome ${species}_clean_0.9.gff | fasta_formatter > ${species}_mRNA_0.9.fasta

    data_95=`bioawk -c fastx '{ print $name, length($seq) }' ${species}_mRNA_0.95.fasta|cut -f2|../../data_summary.sh |cut -f2`
    data_9=`bioawk -c fastx '{ print $name, length($seq) }' ${species}_mRNA_0.9.fasta|cut -f2|../../data_summary.sh |cut -f2`
    data_05=`bioawk -c fastx '{ print $name, length($seq) }' ${species}_mRNA_0.5.fasta|cut -f2|../../data_summary.sh |cut -f2`
    data_03=`bioawk -c fastx '{ print $name, length($seq) }' ${species}_mRNA_0.3.fasta|cut -f2|../../data_summary.sh |cut -f2`
    comp_gene=$(grep -w "${species}" /data/home/s2215768/fly_annotation/complement_annotations/gene_number.tsv|cut -f2)
    echo -e "${species}\t${comp_gene}\t${data_95}\t${data_9}\t${data_05}\t${data_03}" >> ../../summary4.txt

    rm -rf *report.txt  *.agat.log log.tmp*
    cd ../../
    fi

}

export -f mask_genes
>summary4.txt
parallel -j 5 -k --colsep '\t' mask_genes {1} {2} :::: $1

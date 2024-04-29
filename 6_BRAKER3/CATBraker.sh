#!/bin/bash

# This script is used to complement the Braker annotation with CAT annotation

set -e
set -o pipefail
set -u

genome_path="/data/home/s2215768/fly_annotation/data/CAT_genomes/selected_genomes"
species=$1
genome=$2
Refseq=$3
braker_gff="/data/home/s2215768/fly_annotation/braker_annotation/${species}/braker.gff3"
cat_gff="/data/home/s2215768/fly_annotation/CAT_annotation/CAT_gffs/${species}.gff3"

# Check if all the required arguments are provided
if [ $# -lt 3 ]; then
    echo "Usage: $0 <species> <genome> <Refseq>"
    exit 1
fi

# Check if the gff files exist
if [ ! -s "$braker_gff" ]; then
    echo "Braker gff file not found: $braker_gff"
    exit 1
fi

if [ ! -s "$cat_gff" ]; then
    echo "CAT gff file not found: $cat_gff"
    exit 1
fi

# Create a temporary directory

mkdir -p $species
mkdir -p ${species}/tmp
temp_dir="$(pwd)/${species}/tmp"

set -x
cd $species
agat_sp_manage_attributes.pl -gff $cat_gff  -att all_attributes -p level3 -o $temp_dir/cat_clean.gff > $temp_dir/output.txt 2>&1
# Filter out overlapping genes
agat_sp_fix_overlaping_genes.pl -f $temp_dir/cat_clean.gff -o $temp_dir/cat_gff_nooverlap.gff > $temp_dir/output.txt 2>&1
agat_sp_fix_overlaping_genes.pl -f $braker_gff -o $temp_dir/braker_gff_nooverlap.gff > $temp_dir/output.txt 2>&1
rm -rf $temp_dir/*report.txt $temp_dir/cat_clean.gff 
# Clean and sort the gff files using genometools
agat_sp_fix_cds_phases.pl --gff $temp_dir/braker_gff_nooverlap.gff --fasta $genome_path/$genome -o $temp_dir/braker_gff.sorted > $temp_dir/output.txt 2>&1
agat_sp_fix_cds_phases.pl --gff $temp_dir/cat_gff_nooverlap.gff --fasta $genome_path/$genome -o $temp_dir/cat_gff.sorted > $temp_dir/output.txt 2>&1
rm -rf $temp_dir/*report.txt 
gffread -y - -C -H -V -J -g $genome_path/$genome $temp_dir/cat_gff.sorted | fasta_formatter > $temp_dir/cat_mRNA.fasta
grep ">" $temp_dir/cat_mRNA.fasta | sed 's/>//' | sort > $temp_dir/cat_keep.tmp

gffread -y - -C -H -V -J -g $genome_path/$genome $temp_dir/braker_gff.sorted | fasta_formatter > $temp_dir/braker_mRNA.fasta
grep ">" $temp_dir/braker_mRNA.fasta | sed 's/>//' | sed 's/\.t[0-9]+*//' > $temp_dir/braker_keep.tmp

# bioawk to get length of CDSs
bioawk -c fastx '{print $name, length($seq)}' $temp_dir/cat_mRNA.fasta > $temp_dir/cat_mRNA_length.tmp
bioawk -c fastx '{print $name, length($seq)}' $temp_dir/braker_mRNA.fasta > $temp_dir/braker_mRNA_length.tmp

# Cat removed transcripts
gffread -y - -C -g $genome_path/$genome $temp_dir/cat_gff.sorted | fasta_formatter | grep ">" | sed 's/>//' | sort > $temp_dir/cat_all.tmp
grep -Fwvf $temp_dir/cat_keep.tmp $temp_dir/cat_all.tmp > $temp_dir/cat_remove.tmp

# Braker gff filter transcript IDs that are too short (cds<150) and have score < 0.3
awk '$3=="mRNA"' $temp_dir/braker_gff.sorted | cut -f4-6,9 | awk -v FS='\t' -v OFS='\t' '{print $2-$1+1, $0}' | \
sort -k1,1n | awk '$1 > 150 && ($4 == "." || $4 > 0.3)' | sort -t$'\t' -k1,1n -k4,4n | \
awk -F';' -v OFS=';' '{sub(/\.t[0-9]+/, "", $1); print $0}' | grep -o 'ID=[^;]*' | cut -d= -f2 >> $temp_dir/braker_keep.tmp

sort -u $temp_dir/braker_keep.tmp > $temp_dir/braker_keep.tmp.tmp && mv $temp_dir/braker_keep.tmp.tmp $temp_dir/braker_keep.tmp

# Agat to filter the gff files
agat_sp_filter_feature_from_kill_list.pl --gff $temp_dir/cat_gff.sorted --kill_list $temp_dir/cat_remove.tmp --output $temp_dir/cat_gff.clean > $temp_dir/output.txt 2>&1
rm $temp_dir/cat_keep.tmp $temp_dir/cat_all.tmp  
agat_sp_filter_feature_from_keep_list.pl --gff $temp_dir/braker_gff.sorted --keep_list $temp_dir/braker_keep.tmp --output $temp_dir/braker_gff.clean > $temp_dir/output.txt 2>&1 
rm $temp_dir/braker_keep.tmp $temp_dir/braker_gff.sorted
rm -rf $temp_dir/*report.txt 
# Get gff of removed CAT transcripts 
agat_sp_filter_feature_from_keep_list.pl --gff $temp_dir/cat_gff.sorted --keep_list $temp_dir/cat_remove.tmp --output $temp_dir/cat_gff.trash > $temp_dir/output.txt 2>&1
rm $temp_dir/cat_remove.tmp $temp_dir/cat_gff.sorted
rm -rf $temp_dir/*report.txt 
awk -v FS='\t' -v OFS='\t' '{if($3=="transcript") $3="mRNA"; print $0}' $temp_dir/cat_gff.clean > $temp_dir/cat_gff.modified && rm $temp_dir/cat_gff.clean
# Recover transcripts that were not found by braker
awk -v FS='\t' -v OFS='\t' '{if($3=="transcript") $3="mRNA"; print $0}' $temp_dir/cat_gff.trash > $temp_dir/cat_gff.trash.modified && rm $temp_dir/cat_gff.trash


if [ $Refseq == "No" ];then 
# compare annotations (FOR CAT annotations only)
agat_sp_compare_two_annotations.pl -gff1 $temp_dir/braker_gff.clean -gff2 $temp_dir/cat_gff.modified -o $temp_dir/compare_original > $temp_dir/output.txt 2>&1
rm -rf $temp_dir/*report.txt
# compare CDS lengths if CDS overlapa between Braker and CAT
cat $temp_dir/compare_original/gene@mrna@cds_1_1_id_list.txt |tr -d ' '|tr "|" "\t" |grep -v '^$' | sort -u > $temp_dir/compare_geneID.tmp
compare_CDS(){
braker_id=$1
cat_id=`echo $2|tr [a-z] [A-Z]`
temp_dir=$3

cat_len=$(grep -Fwf <(awk -F'\t' '$3=="mRNA" && $9 ~ /Parent='"$cat_id"'/' $temp_dir/cat_gff.modified |grep -o 'ID=[^;]*' | cut -d= -f2) $temp_dir/cat_mRNA_length.tmp | sort -k2,2nr|head -1|cut -f2)
braker_len=$(grep -Fw "$braker_id" $temp_dir/braker_mRNA_length.tmp | sort -k2,2nr|head -1|cut -f2)

if [ -z "$cat_len" ]; then
    cat_len=0
fi

if [ -z "$braker_len" ]; then
    braker_len=0
fi

if [ $braker_len -le $cat_len ]; then
    echo -e "$braker_id" >> $temp_dir/braker_geneID.tmp
else
    echo -e "$cat_id" >> $temp_dir/cat_geneID.tmp
fi

}
export -f compare_CDS
>$temp_dir/braker_geneID.tmp
>$temp_dir/cat_geneID.tmp
parallel -j 4 -k --colsep '\t' compare_CDS {1} {2} {3} ::::  $temp_dir/compare_geneID.tmp ::: $temp_dir

agat_sp_filter_feature_from_kill_list.pl --gff $temp_dir/cat_gff.modified --kill_list $temp_dir/cat_geneID.tmp --output $temp_dir/cat_final.gff > $temp_dir/output.txt 2>&1 
rm $temp_dir/cat_gff.modified 

agat_sp_filter_feature_from_kill_list.pl --gff $temp_dir/braker_gff.clean --kill_list $temp_dir/braker_geneID.tmp --output $temp_dir/braker_final.gff > $temp_dir/output.txt 2>&1 
rm $temp_dir/braker_gff.clean
rm -rf $temp_dir/*report.txt
cp -r $temp_dir/compare_original ./

# Complemet CAT annotation with Braker annotation
if [ -s $temp_dir/cat_gff.trash.modified ]; then
agat_sp_complement_annotations.pl --ref $temp_dir/cat_final.gff --add $temp_dir/braker_final.gff --add $temp_dir/cat_gff.trash.modified --out $temp_dir/complement.gff > $temp_dir/output.txt 2>&1
else
agat_sp_complement_annotations.pl --ref $temp_dir/cat_final.gff --add $temp_dir/braker_final.gff --out $temp_dir/complement.gff > $temp_dir/output.txt 2>&1
fi
rm -rf $temp_dir/*report.txt
else

if [ -s $temp_dir/cat_gff.trash.modified ]; then
agat_sp_complement_annotations.pl --ref $temp_dir/cat_gff.modified --add $temp_dir/braker_gff.clean --add $temp_dir/cat_gff.trash.modified --out $temp_dir/complement.gff > $temp_dir/output.txt 2>&1
else
agat_sp_complement_annotations.pl --ref $temp_dir/cat_gff.modified --add $temp_dir/braker_gff.clean --out $temp_dir/complement.gff > $temp_dir/output.txt 2>&1
fi
rm -rf $temp_dir/*report.txt
fi

# Add IDs to the GFF file
rm -rf ${species}/${species}.gff
agat_sp_manage_IDs.pl --gff $temp_dir/complement.gff --type_dependent --ensembl --collective -o ${species}.gff > $temp_dir/output.txt 2>&1 

rm -r $temp_dir
rm *.log
#echo "### CATBraker for $species completed successfully ###" >> agat.run.out
cd ..
set +x

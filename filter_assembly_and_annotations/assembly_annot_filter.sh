#! /bin/bash

#NCBI dataset to download genome metadata (dehydrated)
#input: file containing list of species name
#file=$1
#oldIFS=$IFS
#IFS=$'\n'
#for sp in `cat $file`
# do
 #taxid=`echo $sp |taxonkit name2taxid`
 #species=`echo $sp |tr ' ' '_'`
 #datasets download genome taxon $taxid --reference --assembly-source refseq --exclude-seq  --exclude-genomic-cds --exclude-protein --exclude-rna --dehydrated --filename ${species}.dataset.zip
 #if [[ -f ${species}.dataset.zip ]]
 #then
 #dataformat tsv genome --fields organism-name,tax-id,assminfo-accession,assminfo-name,assminfo-level,assminfo-sequencing-tech,assmstats-total-number-of-chromosomes,assmstats-number-of-contigs,assmstats-number-of-scaffolds,assmstats-contig-n50,assmstats-scaffold-n50,annotinfo-busco-complete,annotinfo-busco-duplicated,annotinfo-busco-fragmented,annotinfo-busco-missing,annotinfo-busco-totalcount,annotinfo-featcount-gene-total --package ${species}.dataset.zip | tail -n +2 >> tmp_dataset.txt
 #rm ${species}.dataset.zip
 #fi
# unzip ${species}.dataset.zip -d ${species}_dataset
# datasets rehydrate --directory ${species}_dataset
#done
#IFS=$oldIFS


#function to download genome from NCBI
#Usage: download_assembly_ftp species_name assembly_acc
download_assembly_ftp() {
    if [[ $2 = GCF_* ]]
    then
    esearch -db assembly -query $2 </dev/null | esummary| grep 'FtpPath type="RefSeq"'|sed 's/<FtpPath type="RefSeq">//'|sed 's_</FtpPath>__'|sed 's_\s__g' \
     |while read -r url ; do
    fname=$(echo $url | grep -o 'GCF_.*' | sed 's/$/_genomic.fna.gz/') ;
    wget -O "${1}.${2}_genomic.fna.gz"  "$url/$fname" ;
    done;
    else
    esearch -db assembly -query $2 </dev/null | esummary | xtract -pattern DocumentSummary -element FtpPath_GenBank \
     | while read -r url ; do
    fname=$(echo $url | grep -o 'GCA_.*' | sed 's/$/_genomic.fna.gz/') ;
    wget -O "${1}.${2}_genomic.fna.gz" "$url/$fname" ;
    done ;
    fi
}




mkdir -p  ~/fly_annotation/data/CAT_genomes
mkdir -p  ~/fly_annotation/data/buscos_all


fields=`echo "organism-name,tax-id,assminfo-accession,assminfo-name,assminfo-biosample-accession,assminfo-level,assminfo-sequencing-tech,assmstats-total-number-of-chromosomes,assmstats-number-of-contigs,\
assmstats-number-of-scaffolds,assmstats-contig-n50,assmstats-scaffold-n50,annotinfo-busco-complete,annotinfo-busco-duplicated,annotinfo-busco-fragmented,annotinfo-busco-missing,\
annotinfo-busco-totalcount,annotinfo-featcount-gene-total"`

#reference assemblies metadata
datasets download genome taxon 7214 --reference --exclude-seq --exclude-genomic-cds --exclude-protein --exclude-rna --dehydrated --filename refonly.dataset.zip

#All assemblies
datasets download genome taxon 7214  --exclude-seq --exclude-genomic-cds --exclude-protein --exclude-rna --dehydrated --filename allAssemblies.dataset.zip

#annotated assemblies
#datasets download genome taxon 7214 --annotated --exclude-seq --exclude-genomic-cds --exclude-protein --exclude-rna --dehydrated --filename annotatedAssemblies.dataset.zip


dataformat tsv genome --fields $fields --package allAssemblies.dataset.zip | tail -n +2 > all.tmp
dataformat tsv genome --fields $fields --package refonly.dataset.zip | tail -n +2  > ref_assemblies.tmp

#removing entries of nanopore genomes already available on ncbi
cat already_on_ncbi.txt |cut -f3| while read -r line
 do
 grep $line ref_assemblies.tmp >> tmp1
 grep $line all.tmp >> tmp2
 done
grep -Fxvf tmp1 ref_assemblies.tmp > NCBI_assemblies.txt
grep -Fxvf tmp2 all.tmp >  all_assemblies.tmp
rm tmp1 tmp2 all.tmp ref_assemblies.tmp


cat NCBI_assemblies.txt |cut -f1 > file1_refonly.tmp
cat all_assemblies.tmp |cut -f1 > file2_all.tmp

#extracting species without reference genome
grep -Fxvf file1_refonly.tmp file2_all.tmp| sort -u > sp_withoutref.tmp


oldIFS=$IFS
IFS==$'\n'
for i in `cat sp_withoutref.tmp`
 do
 cat all_assemblies.tmp | grep "$i" |sort -t$'\t' -k6 -nr -k11 | head -1 >> NCBI_assemblies.txt
done
IFS=$oldIFS


cat Nanopore_assemblies.txt|cut -f1|sort > f1.tmp
cat NCBI_assemblies.txt|cut -f1| sort -u > f2.tmp

#species with both nanopore and refseq assemblies
grep -Fxf f1.tmp f2.tmp > busco_sp.txt

#species with only nanopore assembly
grep -Fxvf f2.tmp f1.tmp > only_nanopore.tmp
oldIFS=$IFS
IFS==$'\n'
for s in `cat only_nanopore.tmp`
 do
  assembly=`grep "${s}" Nanopore_assemblies.txt |cut -f2`
  sp=`echo ${s}|tr ' ' '_'`
  cp ~/fly_annotation/data/Nanopore_assemblies/${assembly} ~/fly_annotation/data/CAT_genomes/${sp}.nanopore_genomic.fna.gz
done
IFS=$oldIFS

#species with only refseq assembly
grep -Fxvf f1.tmp f2.tmp > only_refseq.tmp

#download ncbi assemblies
awk -F '\t' 'FILENAME=="only_refseq.tmp"{A[$1]=$1} FILENAME=="NCBI_assemblies.txt"{if(A[$1]){print}}' only_refseq.tmp NCBI_assemblies.txt > download.tmp
 oldIFS=$IFS
 IFS==$'\n'
 for line in `cat download.tmp`
 do
  
  acc=`echo $line|awk 'BEGIN {FS=OFS="\t"} {print $3}'`
  sp=`echo $line|awk 'BEGIN {FS=OFS="\t"} {print $1}'|tr ' ' '_'`
  download_assembly_ftp $sp $acc
  mv ${sp}.${acc}_genomic.fna.gz ~/fly_annotation/data/CAT_genomes/
 done
 IFS=$oldIFS
rm download.tmp

#BUSCO comparison of assemblies

scp NCBI_assemblies.txt Nanopore_assemblies.txt s2215768@basden.bio.ed.ac.uk:~/
scp NCBI_assemblies.txt Nanopore_assemblies.txt s2215768@rafflesia.bio.ed.ac.uk:~/
scp NCBI_assemblies.txt Nanopore_assemblies.txt s2215768@sarlacc.bio.ed.ac.uk:~/

IFS==$'\n'

## >>> Submit parallel jobs on local(vera) <<<<
tail -n +66 busco_sp.txt|parallel --jobs 4 "./busco_compare_local.sh {}"

## >>> Submit parallel jobs on remote servers(rafflesia, sarlacc, basden)  <<<

#export download_assembly_ftp
#head -65 busco_sp.txt | parallel --env download_assembly_ftp -S 2/s2215768@basden.bio.ed.ac.uk,4/s2215768@rafflesia.bio.ed.ac.uk,3/s2215768@sarlacc.bio.ed.ac.uk --transferfile busco_compare_remote.sh "./busco_compare_remote.sh {}" 


IFS=$oldIFS

rm *tmp*

######
# Annotation filter




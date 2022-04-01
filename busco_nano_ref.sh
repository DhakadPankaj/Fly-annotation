#! /bin/bash

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

mkdir -p ~/fly_annotation/data/CAT_genomes
mkdir -p ~/fly_annotation/data/buscos_all

if [[ -f busco_summary_only_sp.csv ]]
 then
 echo "busco_summary_only_sp.csv already exists"
 else
 echo -e "species\tNanopore_busco\tRefseq_busco" > ./busco_summary_only_sp.csv
fi


oldIFS=$IFS
IFS==$'\n'
for sp in `cat $1`
do

 if cat NCBI_assemblies.txt |cut -f1|grep -qx "$sp"
 then

  acc=`grep "$sp" NCBI_assemblies.txt|awk 'BEGIN {FS=OFS="\t"} {print $3}'`
  spcy=`grep "$sp" NCBI_assemblies.txt|awk 'BEGIN {FS=OFS="\t"} {print $1}'|tr ' ' '_'`
  download_assembly_ftp $spcy $acc
  zcat ${spcy}.${acc}_genomic.fna.gz > ${spcy}.${acc}_genomic.fna
  busco -i ${spcy}.${acc}_genomic.fna -c 12 -o refseq.busco_${spcy} -m geno -l ~/diptera_odb10 -f --augustus_species fly --offline
  refseq_busco=`cat refseq.busco_${spcy}/short_summary.specific.diptera_odb10.refseq.busco_${spcy}.txt |grep -A2 "***** Results: *****"|tail -1|sed 's_\s__g'`
  echo -e "${spcy}\t\t${refseq_busco}" >> ./busco_summary_only_sp.csv
  tar cf - refseq.busco_${spcy} | pigz -p12 > ${spcy}.refseq.busco.tar.gz && mv ${spcy}.refseq.busco.tar.gz ~/fly_annotation/data/buscos_all/
  mv ${spcy}.${acc}_genomic.fna.gz ~/fly_annotation/data/CAT_genomes/
  rm -r refseq.busco_${spcy} ${spcy}.${acc}_genomic.fna

 elif cat Nanopore_assemblies.txt |cut -f1|grep -qx "$sp"
 then

  spcy=`grep "$sp" Nanopore_assemblies.txt| cut -f1|tr ' ' '_'`
  nano_assembly=`grep "$sp" Nanopore_assemblies.txt|cut -f2`
  zcat ~/fly_annotation/data/Nanopore_assemblies/${nano_assembly} > ${spcy}.nanopore_genomic.fna
  busco -i ${spcy}.nanopore_genomic.fna -c 12 -o nano.busco_${spcy} -m geno -l ~/diptera_odb10 -f --augustus_species fly --offline
  nano_busco=`cat nano.busco_${spcy}/short_summary.specific.diptera_odb10.nano.busco_${spcy}.txt |grep -A2 "***** Results: *****"|tail -1|sed 's_\s__g'`
  echo -e "${spcy}\t${nano_busco}" >> ./busco_summary_only_sp.csv
  tar cf - nano.busco_${spcy} | pigz -p12 > ${spcy}.nano.busco.tar.gz && mv ${spcy}.nano.busco.tar.gz ~/fly_annotation/data/buscos_all/
  mv ~/fly_annotation/data/Nanopore_assemblies/${nano_assembly} ~/fly_annotation/data/CAT_genomes/${spcy}.nanopore_genomic.fna.gz
  rm -r  nano.busco_${spcy} ${spcy}.nanopore_genomic.fna

  fi

done

IFS=$oldIFS

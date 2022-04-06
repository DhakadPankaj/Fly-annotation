

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


oldIFS=$IFS
IFS==$'\n'
for sp in `cat $1`
do
 for line in `grep "$sp" NCBI_assemblies.txt`
 do
  acc=`echo $line|awk 'BEGIN {FS=OFS="\t"} {print $3}'`
  spcy=`echo $line|awk 'BEGIN {FS=OFS="\t"} {print $1}'|tr ' ' '_'`
  download_assembly_ftp $spcy $acc
  zcat ${spcy}.${acc}_genomic.fna.gz > ${spcy}.${acc}_genomic.fna

  nano_assembly=`grep "$sp" Nanopore_assemblies.txt|cut -f2`
  zcat ~/fly_annotation/data/Nanopore_assemblies/${nano_assembly} > ${spcy}.nanopore_genomic.fna

  busco -i ${spcy}.${acc}_genomic.fna -c 12 -o refseq.busco_${spcy} -m geno -l ~/diptera_odb10 -f --augustus_species fly --offline
  busco -i ${spcy}.nanopore_genomic.fna -c 12 -o nano.busco_${spcy} -m geno -l ~/diptera_odb10 -f --augustus_species fly --offline

##busco comparison
  refseq_busco=`cat refseq.busco_${spcy}/short_summary.specific.diptera_odb10.refseq.busco_${spcy}.txt |grep -A2 "***** Results: *****"|tail -1|awk 'BEGIN { FS=OFS=":" } {print $2}'|tr '%[S' ' '|sed 's_\s__g'`
  nano_busco=`cat nano.busco_${spcy}/short_summary.specific.diptera_odb10.nano.busco_${spcy}.txt |grep -A2 "***** Results: *****"|tail -1|awk 'BEGIN { FS=OFS=":" } {print $2}'|tr '%[S' ' '|sed 's_\s__g'`


  if (( $(echo "$nano_busco < $refseq_busco" | bc -l) ))
  then
  mv ${spcy}.${acc}_genomic.fna.gz ~/fly_annotation/data/CAT_genomes/
  rm ${spcy}.${acc}_genomic.fna ${spcy}.nanopore_genomic.fna
  elif (( $(echo "$nano_busco >= $refseq_busco" | bc -l) ))
  then
  mv ~/fly_annotation/data/Nanopore_assemblies/${nano_assembly} ~/fly_annotation/data/CAT_genomes/${spcy}.nanopore_genomic.fna.gz
  rm ${spcy}.nanopore_genomic.fna ${spcy}.${acc}_genomic.fna*
  fi

  tar cf - refseq.busco_${spcy} | pigz -p12 > ${spcy}.refseq.busco.tar.gz && mv ${spcy}.refseq.busco.tar.gz ~/fly_annotation/data/buscos_all/
  tar cf - nano.busco_${spcy} | pigz -p12 > ${spcy}.nano.busco.tar.gz && mv ${spcy}.nano.busco.tar.gz ~/fly_annotation/data/buscos_all/

 rm -r  nano.busco_${spcy} refseq.busco_${spcy}
 done
done
IFS=$oldIFS

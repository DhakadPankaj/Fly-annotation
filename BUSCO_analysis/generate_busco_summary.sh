#! /bin/bash

threads=40
busco_path="~/fly_annotation/data/buscos_all/"

#parallel -j${threads} 'pigz -dc -p 1 {} | tar xf -' ::: $(ls -lah ${busco_path}* | awk '{print $NF}' | tr '\n' ' ')

#echo -e "species\tNanopore_busco\tRefseq_busco" > busco_summary.csv

oldIFS=$IFS
IFS==$'\n'

sp_list=$1
for i in `cat $sp_list`
 do
 sp=`echo $i|tr ' ' '_'`
 nano_busco=`cat nano.busco_${sp}/short_summary.specific.diptera_odb10.nano.busco_${sp}.txt |grep -A2 "***** Results: *****"|tail -1|sed 's_\s__g'`
 refseq_busco=`cat refseq.busco_${sp}/short_summary.specific.diptera_odb10.refseq.busco_${sp}.txt |grep -A2 "***** Results: *****"|tail -1|sed 's_\s__g'`
 echo -e "${sp}\t${nano_busco}\t${refseq_busco}" >> busco_summary.csv
done

IFS=$oldIFS

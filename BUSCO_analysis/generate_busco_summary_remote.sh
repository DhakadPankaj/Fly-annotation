#! /bin/bash

threads=40
busco_path="~/buscos_all/"

mkdir -p busco_summary
cd busco_summary
parallel -j${threads} 'pigz -dc -p 1 {} | tar xf -' ::: $(ls -lah ${busco_path}* | awk '{print $NF}' | tr '\n' ' ')

if [[ -f busco_summary_only_sp.csv ]]
 then
 echo "busco_summary_remote.csv already exists"
 else
 echo -e "species\tNanopore_busco\tRefseq_busco" > ./busco_summary_remote.csv
fi


oldIFS=$IFS
IFS==$'\n'

sp_list=$1
for i in `cat $sp_list`
 do
 sp=`echo $i|tr ' ' '_'`
 nano_busco=`cat nano.busco_${sp}/short_summary.specific.diptera_odb10.nano.busco_${sp}.txt |grep -A2 "***** Results: *****"|tail -1|sed 's_\s__g'`
 refseq_busco=`cat refseq.busco_${sp}/short_summary.specific.diptera_odb10.refseq.busco_${sp}.txt |grep -A2 "***** Results: *****"|tail -1|sed 's_\s__g'`
 echo -e "${sp}\t${nano_busco}\t${refseq_busco}" >> busco_summary.csv
 rm 
done

IFS=$oldIFS

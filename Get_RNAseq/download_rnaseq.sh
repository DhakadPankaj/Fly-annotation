#! /bin/bash

sp_list=$1


#ENA API to download metadata for Drosophila RNAseq. Search rules/fields are saved at ENA Rulespace(Elixir login)

#curl link without aspera links and md5 checksum
#curl -o rnaseq_fly_rawmetadata.tsv -X POST -H "Content-Type: application/x-www-form-urlencoded" -d 'result=read_run&query=tax_tree(7214)%20AND%20library_source%3D%22TRANSCRIPTOMIC%22%20AND%20library_strategy%3D%22RNA-Seq%22%20AND%20library_layout%3D%22PAIRED%22&fields=tax_id%2Cscientific_name%2Csub_species%2Caccession%2Cstudy_accession%2Crun_accession%2Csample_accession%2Csubmission_accession%2Cexperiment_accession%2Cdescription%2Cproject_name%2Ctissue_type%2Csample_description%2Ccollection_date%2Cdev_stage%2Csex%2Cisolation_source%2Cinstrument_platform%2Clibrary_layout%2Clibrary_source%2Clibrary_strategy%2Clibrary_selection%2Clibrary_name%2Csra_bytes%2Cbase_count%2Cread_count&limit=0&format=tsv' "https://www.ebi.ac.uk/ena/portal/api/search"

#curl api link with aspera links and md5 checksum
curl -o rnaseq_fly_rawmetadata.tsv -X POST -H "Content-Type: application/x-www-form-urlencoded" -d 'result=read_run&query=tax_tree(7214)%20AND%20library_source%3D%22TRANSCRIPTOMIC%22%20AND%20library_strategy%3D%22RNA-Seq%22%20AND%20library_layout%3D%22PAIRED%22&fields=tax_id%2Cscientific_name%2Csub_species%2Caccession%2Cstudy_accession%2Crun_accession%2Csample_accession%2Csubmission_accession%2Cexperiment_accession%2Cdescription%2Cproject_name%2Ctissue_type%2Csample_description%2Ccollection_date%2Cdev_stage%2Csex%2Cisolation_source%2Cinstrument_platform%2Clibrary_layout%2Clibrary_source%2Clibrary_strategy%2Clibrary_selection%2Clibrary_name%2Csra_bytes%2Cbase_count%2Cread_count%2Cfastq_aspera%2Cfastq_md5&format=tsv' "https://www.ebi.ac.uk/ena/portal/api/search"

cat rnaseq_fly_rawmetadata.tsv |grep  "_1.fastq.gz" > Only_paired_splitted_fastq_rnaseqData.tsv

rm rnaseq_fly_rawmetadata.tsv

cat $sp_list|cut -f1|taxonkit name2taxid > tmp.txt

awk -F '\t' 'FILENAME=="tmp.txt"{A[$2]=$2} FILENAME=="Only_paired_splitted_fastq_rnaseqData.tsv"{if(A[$1]){print}}' tmp.txt Only_paired_splitted_fastq_rnaseqData.tsv |\
 cut -f 1,2,5,6,12,13,15,16,18,26,27,28 |sort -n -k 1 > tmp2.txt


oldIFS=$IFS
IFS=$'\n'
for i in `awk -F '\t' '{print $2}' < tmp2.txt|uniq` 
 do
  name=`echo $i | tr ' ' '_'`
  mkdir -p $name 
  grep "$i" tmp2.txt > $name/$name.txt

 
  cd $name
  count=`cat $name.txt | wc -l`
  grep -iF -h -e "carcass" -e "Testes" -e "testis" -e "ovary" -e "ovaries" -e "body" -e "brain" -e "fat" -e "head" \
   -e "abdomen" -e "embryo" -e "larva" -e "larvae" -e "thorax" -e "adult" -e "salivary" -e "digestive" -e "gut" -e "nervous"  $name.txt > tissue_var.txt
  tissue_count=`cat tissue_var.txt| wc -l`
  cat tissue_var.txt|sort -t$'\t' -k10 -nr > size_sorted.txt
  
  for i in carcass 'Teste(s|is)' 'ovar(y|ies)' body '(brain|head)' 'fat' abdomen embryo 'larv(a|ae)' thorax  salivary digestive gut nervous
  do 
  grep -i -h -E $i $name.txt | sort -t$'\t' -k10 -nr >> $i_tmp.txt
  if [[ $(grep -i -h -E -w female $i_tmp.txt|grep -v -i -w male | wc -l) -ge "1" ]]
   then 
   grep -i -h -E -w female carcass_tmp.txt|grep -v -i -w male|head -1 >> read_sorted.tx
   else
   cat $i_tmp.txt |head -1 >> read_sorted.txt
  fi 
  if [[ $(grep -i -h -E -w male $i_tmp.txt|grep -v -i -w female | wc -l) -ge "1" ]]
   then
   grep -i -h -E -w male carcass_tmp.txt|grep -v -i -w female|head -1 >> read_sorted.tx 
   else
   cat $i_tmp.txt |head -1 >> read_sorted.txt
  fi
  rm $i_tmp.txt
  done
  cat read_sorted.txt|sort -t$'\t' -k10 -nr | uniq > read_uniq.txt

  readsorted=`cat read_uniq.txt| wc -l`

  if [[ $count -le 10 ]]
   then
   cat $name.txt | sort -t$'\t' -k10 -nr| cut -f 4 > sra.txt
    cat sra.txt |xargs -I sra ~/softwares/enaBrowserTools-0.0.3/python3/./enaDataGet -f fastq --aspera sra
    
    elif [[ $tissue_count -gt 10 ]]  && [[ $readsorted -ge 10 ]]
      then
      awk '!x[$0]++' read_sorted.txt| head -10 | cut -f 4 > sra.txt
      cat sra.txt |xargs -I sra ~/softwares/enaBrowserTools-0.0.3/python3/./enaDataGet -f fastq --aspera sra
     
     elif [[ $tissue_count -ge 10 ]]  && [[ $readsorted -lt 10 ]]
      then
      cat size_sorted.txt >> read_uniq.txt
      awk '!x[$0]++' read_uniq.txt | head -10 | cut -f 4 > sra.txt
      cat sra.txt |xargs -I sra ~/softwares/enaBrowserTools-0.0.3/python3/./enaDataGet -f fastq --aspera sra 
 
     elif [[ $tissue_count -lt 10 ]] && [[ $count -gt 10 ]]
      then
      cat  $name.txt | sort -t$'\t' -k10 -nr >> size_sorted.txt
      awk '!x[$0]++' size_sorted.txt | head -10 | cut -f 4 > sra.txt
      cat sra.txt |xargs -I sra ~/softwares/enaBrowserTools-0.0.3/python3/./enaDataGet -f fastq --aspera sra
     
   fi 


#md5 checksum and downloading retry loops for failed aspera downloads

for i in `cat sra.txt`
do

cd $i

#ls *.gz | xargs -n 1 gunzip -t 2>&1 | cut -f 2 -d: - | xargs -t -n 1 rm

loop=0

for seconds in 60 120 240
do

ena_md5_F1=`grep "$i" ../../Only_paired_splitted_fastq_rnaseqData.tsv | cut -f 28| sed 's_;_\t_'| cut -f1`
ena_md5_F2=`grep "$i" ../../Only_paired_splitted_fastq_rnaseqData.tsv | cut -f 28| sed 's_;_\t_'| cut -f2`

 if [[ -f "$i"_1.fastq.gz ]] && [[ -f "$i"_2.fastq.gz ]]
 then
 md5_F1=`md5sum "$i"_1.fastq.gz |sed 's_\s_\t_'|cut -f 1`
 md5_F2=`md5sum "$i"_2.fastq.gz |sed 's_\s_\t_'|cut -f 1`
        if [[ $md5_F1 == $ena_md5_F1 ]] && [[ $md5_F2 == $ena_md5_F2 ]]
         then
         echo "md5 checked at loop:$loop & it same as ENA file"
         else
         rm *.fastq.gz *logs
        fi
 fi

 if [[ -f "$i"_1.fastq.gz ]] && [[ -f "$i"_2.fastq.gz ]]
  then
  echo "Both $i fastq are downloaded"
  break
  else
  rm *.fastq.gz *logs
  ~/softwares/enaBrowserTools-0.0.3/python3/./enaDataGet -f fastq --aspera $i
  mv $i/* ./
  rm -r $i/
 fi

 sleep $seconds

loop=$((loop+1))

done

cd ..

done

rm tissue_var.txt  
rm size_sorted.txt
rm read_uniq.txt
rm read_sorted.txt

cd ..

done

rm tmp.txt tmp2.txt

IFS=$OLDIFS

genome="~/fly_annotation/data/CAT_genomes/"
busco_path="/data/home/s2215768/fly_annotation/data/buscos_all"
threads=40
mkdir -p ~/fly_annotation/data/rnaseq_all/assembly_annot_check/busco_summary
cd ~/fly_annotation/data/rnaseq_all/assembly_annot_check/busco_summary/

#parallel -j${threads} 'pigz -dc -p 1 {} | tar xf -' ::: $(ls -lah ${busco_path}* | awk '{print $NF}' | tr '\n' ' ')

if [[ -f busco_summary_only_sp.tsv ]]
 then
 echo "busco_summary_only_sp.csv already exists"
 else
 echo -e "species\tNanopore_busco\tRefseq_busco\tAssembly_name" >  ~/fly_annotation/data/rnaseq_all/assembly_annot_check/busco_summary_only_sp.tsv
fi


oldIFS=$IFS
IFS==$'\n'
for sp in `cat $1`
do

 if cat ../NCBI_assemblies.txt |cut -f1|grep -qx "$sp"
 then

  acc=`grep "$sp" ../NCBI_assemblies.txt|awk 'BEGIN {FS=OFS="\t"} {print $3}'`
  spcy=`grep "$sp" ../NCBI_assemblies.txt|awk 'BEGIN {FS=OFS="\t"} {print $1}'|tr ' ' '_'`
  pigz -dc ${busco_path}/${spcy}.refseq.busco.tar.gz | tar xf -
  refseq_busco=`cat refseq.busco_${spcy}/short_summary.specific.diptera_odb10.refseq.busco_${spcy}.txt |grep -A2 "***** Results: *****"|tail -1|sed 's_\s__g'`
  echo -e "${spcy}\tNA\t${refseq_busco}\t${spcy}.${acc}_genomic.fna.gz" >>  ~/fly_annotation/data/rnaseq_all/assembly_annot_check/busco_summary_only_sp.tsv
  rm refseq.busco_${spcy}
 elif cat ../Nanopore_assemblies.txt |cut -f1|grep -qx "$sp"
 then

  spcy=`grep "$sp" ../Nanopore_assemblies.txt| cut -f1|tr ' ' '_'`
  nano_assembly=${spcy}.nanopore_genomic.fna
  pigz -dc ${busco_path}/${spcy}.nano.busco.tar.gz |tar xf -
  nano_busco=`cat nano.busco_${spcy}/short_summary.specific.diptera_odb10.nano.busco_${spcy}.txt |grep -A2 "***** Results: *****"|tail -1|sed 's_\s__g'`
  echo -e "${spcy}\t${nano_busco}\tNA\t${nano_assembly}" >>  ~/fly_annotation/data/rnaseq_all/assembly_annot_check/busco_summary_only_sp.tsv
  rm nano.busco_${spcy} 
 fi

done

cd ../

IFS=$oldIFS


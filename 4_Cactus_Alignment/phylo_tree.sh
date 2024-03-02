#! /bin/bash

#This scripts builds maximum likelihood phylogenetic tree using BUSCO genes.

threads=60
busco_path="/data/home/s2215768/fly_annotation/data/busco_data/"
wd=`pwd`

cd $busco_path
#parallel -j${threads} 'pigz -dc -p 1 {} | tar xf -' ::: $(ls -lah ${busco_path}* | awk '{print $NF}' | tr '\n' ' ')
#rm *tar.gz


#Set the maximum number of taxa BUSCOs can be dup/missing/fragmented in
species_num=`echo $(ls -1d *nano.busco* | wc -l)+$(ls -1d *refseq.busco* | wc -l)|bc`
missing_no=8

#extract names of complete BUSCO genes predicted from all genomes
if [ -f complete_buscos.txt ]; then
    rm complete_buscos.txt
fi

for file in $(find ${busco_path} -name "full_table.tsv"); do
 grep -v "^#" ${file} | awk '$2=="Complete" {print $1}' >> complete_buscos.txt;
done

#filter out the complete BUSCOs that are only present less than 3 genomes
sort complete_buscos.txt | uniq -c > complete_busco_ids_with_counts.txt
awk -v spnum="${species_num}" -v missing="${missing_no}" '$1 >= (spnum - missing) {print $2}' complete_busco_ids_with_counts.txt > tmp
cat tmp | sort | uniq > final_busco_ids.txt
rm tmp

echo "Number of BUSCOs used: $(cat final_busco_ids.txt | wc -l)"

#copy exon nucleotide sequences

mkdir -p busco_nt
for dir in $(find . -type d -name "single_copy_busco_sequences"); do
  #get species names from archived busco folders
  IFS='/' read -ra a <<< "$dir"
  sppname=$(echo ${a[1]} | sed 's/\(nano\|refseq\).busco_//')

  #alternative
  comp () {
    file=${1}/${2}.fna
      if [ -f $file ]; then
          #clean up fasta headers
          (echo $(echo ">${3}" | tr '.' '_' | tr '[:lower:]' '[:upper:]' );
           cat $file | grep -v ">") > ./busco_nt/${3}_${2}.fna
      fi
  }
  export -f comp
  cat final_busco_ids.txt | parallel -j${threads} comp ${dir} {} ${sppname}
done

#generate one file per BUSCO ID
while read line; do
        cat ./busco_nt/*_${line}.fna >> ./busco_nt/${line}.fasta
done < final_busco_ids.txt
find . -wholename './busco_nt/*.fna' | xargs rm -f #regular rm doesn't work for many files


#only keep BUSCO nt sequences with length%3==0
#why are some seqs not divisible by 3?
#this needs to be addressed later
#throw out BUSCO nt sequences that length%3 != 0
mkdir -p busco_nt_3
while read line; do
    cat ./busco_nt/${line}.fasta | \
        perl -pe '$. > 1 and /^>/ ? print "\n" : chomp' - | \
        awk 'NF>0{x=$0;getline y;if ((length(y) % 3) == 0){print x"\n"y}}' \
        > ./busco_nt_3/${line}.3.fasta
done < final_busco_ids.txt



# code below can be run if excluding any BUSCO with weird CDS
rm -f align_busco_ids.txt
while read line; do
    remain=$(cat ./busco_nt_3/${line}.3.fasta | \
        perl -pe '$. > 1 and /^>/ ? print "\n" : chomp' - | \
        awk 'NF>0{x=$0;getline y;print length(y) % 3}' | \
        awk '{s+=$1} END {print s}')
    if [[ ${remain} == 0 ]]; then
        echo ${line} >> align_busco_ids.txt
    fi
done < final_busco_ids.txt


#cat final_busco_ids.txt > align_busco_ids.txt

#alignment
mkdir -p ./busco_nt_aln
alnthreads="6" #num threads per process
parallel -j$(echo "${threads}/${alnthreads}" | bc) \
    mafft --thread ${alnthreads} --genafpair --maxiterate 1000 {} '>' \
    ./busco_nt_aln/{/.}.aln.fasta ::: $(ls -lah ./busco_nt_3/*.3.fasta | awk '{print $NF}' | tr '\n' ' ')

#throw away failed jobs (big sequence, OOM)
ls -l ./busco_nt_aln/*.fasta | awk '{if($5 != 0){print $NF}}' | \
    awk 'BEGIN {FS="/"} {print $NF}' | sed "s/.3.aln.fasta//" \
    > align_busco_ids.txt

# trim alignments using fasta_site_trim.py
mkdir -p busco_nt_aln_tr
cat align_busco_ids.txt | parallel -j${threads} python ./fasta_site_trim.py --Nbase 3 --input busco_nt_aln/{}.3.aln.fasta
mv ./busco_nt_aln/*.tr ./busco_nt_aln_tr/

#Gene tree reconstruction with iqtree
mkdir -p iqtree
treethreads="6" #num threads per process
cat align_busco_ids.txt | parallel -j$(echo "${threads}/${treethreads}" | bc) \
    iqtree2 -s ./busco_nt_aln_tr/{}.3.aln.fasta.tr -bb 1000 -nt ${treethreads} -m GTR+I+G \
    -blmin 1e-300 -safe -pre iqtree_{} -mem 18G

mv iqtree_* ./iqtree

#concatenate iqtree trees: To generate a single file with ML trees for all BUSCO genes
cat iqtree/*.treefile > busco_gene_trees_ml.tree


#Species tree reconstruction: run astral on best trees
java -jar ~/fly_annotation/tools/ASTRAL/astral.5.7.8.jar \
    --input busco_gene_trees_ml.tree \
    --output busco_species_astral.tree \
    -t 3

#Compute the pairwise RF distances between 2 sets of trees:
if [ -f RF_dist.tmp ]; then
    rm RF_dist.tmp
fi

for tree in iqtree/*.treefile
do
iqtree2 -t busco_species_astral.tree --tree-dist2 ${tree}
rf=`cat busco_species_astral.tree.rfdist |tail -1|sed 's_\s__g'|sed 's/Tree0//'`
echo -e "$tree\t$rf" >> RF_dist.tmp
rm busco_species_astral.tree.rfdist
done
cat RF_dist.tmp | sort -n -k2|cut -f1| head -50|sed 's:iqtree/iqtree_::'|sed 's:.treefile:.3.aln.fasta.tr:' > selected_aln.txt

#copy selected gene alignment
if [ -d selected_aln ]; then
    rm -r selected_aln/
fi

mkdir -p selected_aln
for file in `cat selected_aln.txt`
do
cp ./busco_nt_aln_tr/${file} ./selected_aln/
done

#Concatenation of BUSCO MSAs into a Supermatrix using catfasta2phyml.pl
catfasta2phyml.pl -f -c ./selected_aln/*.3.aln.fasta.tr > out.aln.fasta

#ML phylogenetic inference using IQTREE
mkdir -p iqtree_ML
treethreads="40" #num threads per process
iqtree2 -s out.aln.fasta -nt ${treethreads} -m GTR+I+G -pre ML -safe -bb 1000 -alrt 1000 -abayes
mv ML.* ./iqtree_ML

#cat busco_species_astral.tree | sed -E 's/([A-Z\-]+[A-Z0-9]*)/\1:1/' > busco_species_astral_withtiplengths.tree

cd $wd

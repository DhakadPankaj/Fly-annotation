#! /bin/bash

#cat $1 | while read line || [[ -n $line ]];
#do
#  ~/softwares/enaBrowserTools-0.0.3/python3/./enaDataGet -f fastq --aspera $line
#done

#for i in `cat $1`
#do
# cd $i
# fasterq-dump -p -e50 $i
# cd ..
#done




for sp in `cat species_names.txt | sed 's:\s:_:'` 
do 
cd $sp
echo "${sp} RNA reads mapping starting"
mv */*fastq.gz ./
zcat *_1.fastq.gz > ${sp}_rnaseq_R1.fastq
zcat *_2.fastq.gz > ${sp}_rnaseq_R2.fastq

# BBnorm kmer based approach to remove reads in regions with depth >maxCoverage
bbnorm.sh in=${sp}_rnaseq_R1.fastq in2=${sp}_rnaseq_R2.fastq out=${sp}_rnaseq_norm_R1.fastq out2=${sp}_rnaseq_norm_R2.fastq min=2 target=100 threads=50 -Xmx150G

cat ${sp}_rnaseq_norm_R1.fastq | pigz -p50 > ${sp}_rnaseq_norm_R1.fastq.gz
cat ${sp}_rnaseq_norm_R2.fastq | pigz -p50 > ${sp}_rnaseq_norm_R2.fastq.gz

# STAR for mapping rna-seq
# must use high quality (paired end and >75bp)
cd ../../
mkdir ${sp}_genome_fasta
zcat ${sp}* > ${sp}.assembly.sm.fasta
mkdir ${sp}_genomeDir
mkdir ${sp}_aligned_reads

# generate genome index
STAR --runThreadN 50 --genomeDir ./${sp}_genomeDir/ --runMode genomeGenerate --genomeFastaFiles ${sp}.assembly.sm.fasta

# map reads
STAR --runThreadN 50 --genomeDir ./${sp}_genomeDir/ --readFilesIn ./rnaseq_data/${sp}/${sp}_rnaseq_norm_R1.fastq ./rnaseq_data/${sp}/${sp}_rnaseq_norm_R2.fastq --outFileNamePrefix ./${sp}_aligned_reads/

samtools view -Sb -@$50 ./${sp}_aligned_reads/Aligned.out.sam | sambamba sort -o ${sp}_rnaseq.bam -p -t 50 -m 150G /dev/stdin && samtools index ${sp}_rnaseq.bam 

cd ./rnaseq_data/

echo ${sp}_complete

done

echo Bam!

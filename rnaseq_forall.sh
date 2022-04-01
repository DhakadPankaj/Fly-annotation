#! /bin/bash

sp_list=$1

threads=60
mem="200G"

IFS=$'\n'

for line in `cat $sp_list` 
do 

sp=`echo $line|cut -f1| sed 's:\s:_:g'`
hal=`echo $line|cut -f2|tr a-z A-Z`
assembly_name=`echo $line|cut -f3`


if [[ -d ${sp} ]]
then
echo "BBnorm starts for ${sp}"
cd $sp
mv */*fastq.gz ./
zcat *_1.fastq.gz > ${sp}_rnaseq_R1.fastq
zcat *_2.fastq.gz > ${sp}_rnaseq_R2.fastq

# BBnorm kmer based approach to remove reads in regions with depth >maxCoverage
bbnorm.sh in=${sp}_rnaseq_R1.fastq in2=${sp}_rnaseq_R2.fastq out=${sp}_rnaseq_norm_R1.fastq out2=${sp}_rnaseq_norm_R2.fastq min=2 target=100 threads=${threads} -Xmx${mem}

cat ${sp}_rnaseq_norm_R1.fastq | pigz -p50 > ${sp}_rnaseq_norm_R1.fastq.gz
cat ${sp}_rnaseq_norm_R2.fastq | pigz -p50 > ${sp}_rnaseq_norm_R2.fastq.gz

cd ..

# STAR for mapping rna-seq
# must use high quality (paired end and >75bp)

#assembly_name=`cat ../cactus_alignment/version_3/cactus_setup/cactus_$anc.txt| grep "$hal" |sed '/^# Outgroup:/d'|sed -n 2p| cut -d "/" -f3`

mkdir -p bams

cd bams
zcat ../../cactus_alignment/version_3/genomes/${assembly_name} > ${sp}.assembly.sm.fasta

mkdir ${sp}_genome_fasta
mkdir ${sp}_genomeDir
mkdir ${sp}_aligned_reads

# generate genome index
STAR --runThreadN ${threads} --genomeDir ./${sp}_genomeDir/ --runMode genomeGenerate --genomeFastaFiles ${sp}.assembly.sm.fasta

# map reads
STAR --runThreadN ${threads} --genomeDir ./${sp}_genomeDir/ --readFilesIn ../${sp}/${sp}_rnaseq_norm_R1.fastq ../${sp}/${sp}_rnaseq_norm_R2.fastq --outFileNamePrefix ./${sp}_aligned_reads/

samtools view -Sb -@${threads} ./${sp}_aligned_reads/Aligned.out.sam | sambamba sort -o ${hal}_rnaseq.bam -p -t ${threads} -m ${mem} /dev/stdin && samtools index ${hal}_rnaseq.bam 

echo "${sp} rnaseq read mapping successfully completed"

rm -r ${sp}_genomeDir  ${sp}_aligned_reads ${sp}_genome_fasta ${sp}.assembly.sm.fasta  
rm ../${sp}/*.fastq
cd ..
fi

done


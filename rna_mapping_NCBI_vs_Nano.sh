#! /bin/bash

sp_list=$1

for sp in `cat $sp_list | sed 's:\s:_:'` 
do 

if [[ -d ${sp} ]]
then
echo "BBnorm starts for ${sp}"
cd $sp
#mv */*fastq.gz ./
zcat *_1.fastq.gz > ${sp}_rnaseq_R1.fastq
zcat *_2.fastq.gz > ${sp}_rnaseq_R2.fastq

# BBnorm kmer based approach to remove reads in regions with depth >maxCoverage
bbnorm.sh in=${sp}_rnaseq_R1.fastq in2=${sp}_rnaseq_R2.fastq out=${sp}_rnaseq_norm_R1.fastq out2=${sp}_rnaseq_norm_R2.fastq min=2 target=100 threads=50 -Xmx150G

cat ${sp}_rnaseq_norm_R1.fastq | pigz -p50 > ${sp}_rnaseq_norm_R1.fastq.gz
cat ${sp}_rnaseq_norm_R2.fastq | pigz -p50 > ${sp}_rnaseq_norm_R2.fastq.gz

# STAR for mapping rna-seq
# must use high quality (paired end and >75bp)

#Mapping on NCBI and Nanopore genomes
cd ../../

for assembly in Nanopore_genomes NCBI_genomes
do
cd $assembly

if [[ -f ${sp}.${assembly}.fasta ]]
then
echo "$assembly is available for ${sp}, running STAR for mapping reads"
mkdir ${sp}.${assembly}_genome_fasta
#zcat ${sp}* > ${sp}.assembly.sm.fasta
mkdir ${sp}.${assembly}_genomeDir
mkdir ${sp}.${assembly}_aligned_reads

# generate genome index
STAR --runThreadN 50 --genomeDir ./${sp}.${assembly}_genomeDir/ --runMode genomeGenerate --genomeFastaFiles ${sp}.${assembly}.fasta

# map reads
STAR --runThreadN 50 --genomeDir ./${sp}.${assembly}_genomeDir/ --readFilesIn ../rnaseq_all/${sp}/${sp}_rnaseq_norm_R1.fastq ../rnaseq_all/${sp}/${sp}_rnaseq_norm_R2.fastq --outFileNamePrefix ./${sp}.${assembly}_aligned_reads/

samtools view -Sb -@$50 ./${sp}.${assembly}_aligned_reads/Aligned.out.sam | sambamba sort -o ${sp}.${assembly}_rnaseq.bam -p -t 50 -m 150G /dev/stdin && samtools index ${sp}.${assembly}_rnaseq.bam 

echo "${sp} read mapping on to $assembly successfully completed"

else
echo "$assembly is not available for ${sp}"

fi
cd ..
done

fi

cd ./rnaseq_all/

done


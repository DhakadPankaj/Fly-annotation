#! /bin/bash
 
source ~/miniconda3/etc/profile.d/conda.sh
conda activate annotation

alt_rna_mapping() {
    threads=40
    mem="200G"

    sp=$(echo "$1" | tr '[:upper:]' '[:lower:]' | sed 's/.*/\u&/')
    hal=$1
    assembly_name=$3
    alt_sp=$(echo "$2" | tr '[:upper:]' '[:lower:]' | sed 's/.*/\u&/')

    if [ -s "./bams/${2}_rnaseq.bam" ]; then
        if [[ -f "./bams/${hal}_rnaseq.intron.bam" ]]; then
            echo "BAM already exists for ${sp}"
        elif [[ -d ${alt_sp} ]]; then
            echo "BBnorm starts for ${sp}"

            cd ${alt_sp}
            if [[ -f "${alt_sp}_rnaseq_norm_R1.fastq.gz" ]] && [[ -f "${alt_sp}_rnaseq_norm_R2.fastq.gz" ]]; then
                echo "bbnorm files already exist for ${alt_sp}"
                zcat "${alt_sp}_rnaseq_norm_R1.fastq.gz" > "${alt_sp}_rnaseq_norm_R1.fastq"
                zcat "${alt_sp}_rnaseq_norm_R2.fastq.gz" > "${alt_sp}_rnaseq_norm_R2.fastq"
            else
                mv */*fastq.gz ./
                zcat *_1.fastq.gz > "${alt_sp}_rnaseq_R1.fastq"
                zcat *_2.fastq.gz > "${alt_sp}_rnaseq_R2.fastq"

                # BBnorm kmer based approach to remove reads in regions with depth >maxCoverage
                bbnorm.sh in="${alt_sp}_rnaseq_R1.fastq" in2="${alt_sp}_rnaseq_R2.fastq" out="${alt_sp}_rnaseq_norm_R1.fastq" out2="${alt_sp}_rnaseq_norm_R2.fastq" min=2 target=100 threads=${threads} -Xmx${mem}
            fi
            cd ..

            # STAR for mapping rna-seq
            mkdir -p bams

            cd bams
            cp "/data/home/s2215768/fly_annotation/data/CAT_genomes/selected_genomes/${assembly_name}" "${sp}.assembly.fasta"
            mkdir "${sp}_genome_fasta"
            mkdir "${sp}_genomeDir"
            mkdir "${sp}_aligned_reads"

            # generate genome index
            STAR --runThreadN ${threads} --genomeDir "./${sp}_genomeDir/" --runMode genomeGenerate --limitGenomeGenerateRAM 350000000000 --genomeFastaFiles "${sp}.assembly.fasta"

            # map reads
            STAR --runThreadN ${threads} --genomeDir "./${sp}_genomeDir/" --readFilesIn "../${alt_sp}/${alt_sp}_rnaseq_norm_R1.fastq" "../${alt_sp}/${alt_sp}_rnaseq_norm_R2.fastq" \
                --outFileNamePrefix "./${sp}_aligned_reads/"

            samtools view -Sb -@${threads} "./${sp}_aligned_reads/Aligned.out.sam" | sambamba sort -o "${hal}_rnaseq.intron.bam" -p -t ${threads} -m ${mem} /dev/stdin && samtools index "${hal}_rnaseq.intron.bam"

            if [[ -s "./${hal}_rnaseq.intron.bam" ]]; then
                cat "../${alt_sp}/${alt_sp}_rnaseq_norm_R1.fastq" | pigz -p${threads} > "../${alt_sp}/${alt_sp}_rnaseq_norm_R1.fastq.gz"
                cat "../${alt_sp}/${alt_sp}_rnaseq_norm_R2.fastq" | pigz -p${threads} > "../${alt_sp}/${alt_sp}_rnaseq_norm_R2.fastq.gz"
                echo "${sp} rnaseq read mapping successfully completed"
            else
                echo "${sp} rnaseq read mapping failed"
            fi

            rm -r "${sp}_genomeDir" "${sp}_aligned_reads" "${sp}_genome_fasta" "${sp}.assembly.fasta"
            rm "../${alt_sp}"/*.fastq
            cd ..
        fi
    fi
}

export -f alt_rna_mapping
parallel -j 1 --colsep "\t" alt_rna_mapping :::: alt_rna_sp.tsv

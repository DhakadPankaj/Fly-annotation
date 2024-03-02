#!/bin/bash

# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate annotation

# Set the number of threads and memory
threads=30
mem="300G"

# Set variables
sp_list=$1
clade_id=7214

# Download RNA-seq metadata from ENA
URL="result=read_run&query=tax_tree($clade_id)%20AND%20library_source%3D%22TRANSCRIPTOMIC%22%20AND%20library_strategy%3D%22RNA-Seq%22%20AND%20library_layout%3D%22PAIRED%22&fields=tax_id%2Cscientific_name%2Csub_species%2Caccession%2Cstudy_accession%2Crun_accession%2Csample_accession%2Csubmission_accession%2Cexperiment_accession%2Cdescription%2Cproject_name%2Ctissue_type%2Csample_description%2Ccollection_date%2Cdev_stage%2Csex%2Cisolation_source%2Cinstrument_platform%2Clibrary_layout%2Clibrary_source%2Clibrary_strategy%2Clibrary_selection%2Clibrary_name%2Csra_bytes%2Cbase_count%2Cread_count%2Cfastq_aspera%2Cfastq_md5&format=tsv"
curl -o rnaseq_fly_rawmetadata.tsv -X POST -H "Content-Type: application/x-www-form-urlencoded" -d $URL "https://www.ebi.ac.uk/ena/portal/api/search"

# Filter metadata for paired-end reads
cat rnaseq_fly_rawmetadata.tsv | grep "_1.fastq.gz" > Only_paired_splitted_fastq_rnaseqData.tsv
rm rnaseq_fly_rawmetadata.tsv

# Process RNA-seq data for each species
RNAseq_selector() {
  name=$(echo $1 | tr ' ' '_')
  NAME=$(echo ${name} | tr [a-z] [A-Z])

  if [[ -d ${name} ]] || [[ -f ./bams/${NAME}_rnaseq.bam ]] || [[ -f ./bams/${NAME}_rnaseq.intron.bam ]]; then
    echo "RNA seq data already exists for ${name}"
    echo ${1} >> rnaseq_splist.txt
  else
    grep "^${1}$(printf '\t')" tmp2.txt > $name.txt

    RNA=$(cat $name.txt | wc -l)
    rm $name.txt

    if [[ $RNA -ge 1 ]]; then
      echo $name
      mkdir -p $name
      cd $name
      grep "^${1}$(printf '\t')" ../tmp2.txt > $name.txt
      cp $name.txt ${name}_all.txt
      cat $name.txt | tr " " "_" | awk -F "\t" '$8=="PolyA"' > polyA.tmp
      polyA=$(cat polyA.tmp | wc -l)

      x=$(echo "10" - $polyA | bc)
      if [[ $polyA -ge 10 ]]; then
        cp polyA.tmp $name.txt
      else
        grep -vFf polyA.tmp $name.txt > without_polyA.tmp
        mv without_polyA.tmp $name.txt
      fi

      count=$(cat $name.txt | wc -l)
      tissue_count=$(cat $name.txt | tr " " "_" | sort -t$'\t' -k3,3 -k9nr | awk -F"\t" 'BEGIN { OFS = "\t"} {print $1,$2,$3,$4,$5"_"$7,$6,$7,$8,$9,$10,$11,$12}' | cut -f5,6 | awk -F"\t" '!seen[$6]++' | awk -F"\t" '!seen[$5]++' | wc -l)

      if [[ $count -le 10 ]]; then
        cat polyA.tmp | cut -f 4 | sort -R | head -10 > sra.txt
        cat $name.txt | sort -t$'\t' -k9 -nr | head -$x | cut -f 4 >> sra.txt
        cat sra.txt | sort -u | xargs -I sra ~/softwares/enaBrowserTools-1.7/python3/./enaDataGet -f fastq --aspera sra
      elif [[ $tissue_count -lt 10 ]] && [[ $count -gt 10 ]]; then
        cat polyA.tmp | cut -f 4 | sort -R | head -10 > sra.txt
        cat $name.txt | sort -t$'\t' -k3,3 -k9nr | awk -F"\t" 'BEGIN { OFS = "\t"} {print $1,$2,$3,$4,$5"_"$7,$6,$7,$8,$9,$10,$11,$12}' | awk -F"\t" '!seen[$5]++' | sort -t$'\t' -R | head -$x | cut -f 4 >> sra.txt
        cat sra.txt | sort -u | xargs -I sra ~/softwares/enaBrowserTools-1.7/python3/./enaDataGet -f fastq --aspera sra
      elif [[ $tissue_count -ge 10 ]] && [[ $count -gt 10 ]]; then
        cat polyA.tmp | cut -f 4 | sort -R | head -10 > sra.txt
        cat $name.txt | sort -t$'\t' -k3,3 -k9nr | awk -F"\t" 'BEGIN { OFS = "\t"} {print $1,$2,$3,$4,$5"_"$7,$6,$7,$8,$9,$10,$11,$12}' | awk -F"\t" '!seen[$6]++' | awk -F"\t" '!seen[$5]++' | head -$x | cut -f 4 >> sra.txt
        cat sra.txt | sort -u | xargs -I sra ~/softwares/enaBrowserTools-1.7/python3/./enaDataGet -f fastq --aspera sra
      fi

      >sra_validation.tsv
      for i in $(cat sra.txt | sort -u); do
        cd $i
        loop=0
        for seconds in 60 120 240; do
          ena_md5_F1=$(awk -F "\t" -va=$i '$4==a' ../${name}_all.txt | cut -f12 | sed 's_;_\t_' | cut -f1)
          ena_md5_F2=$(awk -F "\t" -va=$i '$4==a' ../${name}_all.txt | cut -f12 | sed 's_;_\t_' | cut -f2)
          if [[ -f "$i"_1.fastq.gz ]] && [[ -f "$i"_2.fastq.gz ]]; then
            md5_F1=$(md5sum "$i"_1.fastq.gz | sed 's_\s_\t_' | cut -f 1)
            md5_F2=$(md5sum "$i"_2.fastq.gz | sed 's_\s_\t_' | cut -f 1)
            if [[ $md5_F1 == $ena_md5_F1 ]] && [[ $md5_F2 == $ena_md5_F2 ]]; then
              echo "md5 checked at loop:$loop & it is same as ENA file"
            else
              rm -r *.fastq.gz logs/
            fi
          fi
          if [[ -f "$i"_1.fastq.gz ]] && [[ -f "$i"_2.fastq.gz ]]; then
            echo "Both $i fastq are downloaded" >>../sra_validation.tsv
            break
          else
            rm -r *.fastq.gz logs/
            ~/softwares/enaBrowserTools-1.7/python3/./enaDataGet -f fastq --aspera $i
            mv $i/* ./
            rm -r $i/
          fi
          sleep $seconds
          loop=$((loop + 1))
        done
        cd ..
      done

      cd ..
    fi
  fi
}

export -f RNAseq_selector

# Process RNA-seq data in parallel
parallel -j 5 --colsep "\t" RNAseq_selector :::: $sp_list

# Update BAM_sp.txt file
PolyA_update() {
  name=$(echo $1 | tr ' ' '_')
  NAME=$(echo ${name} | tr [a-z] [A-Z])

  polyA=$(cat ${name}/polyA.tmp | wc -l)
  sra=$(cat ${name}/sra.txt | wc -l)
  if [[ $polyA -ge 2 ]] || [[ $sra -ge 5 ]]; then
    echo -e "${NAME}\t${polyA}\t${sra}" >>BAM_sp.txt
  else
    mv bams/${NAME}_rnaseq.bam bams/${NAME}_rnaseq.intron.bam
    mv bams/${NAME}_rnaseq.bam.bai bams/${NAME}_rnaseq.intron.bam.bai
  fi
}

export -f PolyA_update

# Update BAM_sp.txt file in parallel
>BAM_sp.txt
parallel -j 5 --colsep "\t" PolyA_update :::: rnaseq_splist.txt

##############################################
##           RNA-seq mapping                ##
##############################################

# Iterate through each line in the sp_list file
while IFS=$'\n' read -r line; do
  sp=$(echo "$line" | cut -f1 | tr ' ' '_')
  hal=$(echo "$line" | cut -f1 | tr 'a-z' 'A-Z')
  assembly_name=$(echo "$line" | cut -f3)

  # Check if BAM file already exists
  if [[ -f "./bams/${hal}_rnaseq.bam" ]] || [[ -f "./bams/${NAME}_rnaseq.intron.bam" ]]; then
    echo "BAM already exists for $sp"
  elif [[ -d "$sp" ]]; then
    echo "BBnorm starts for $sp"
    cd "$sp"

    # Check if normalized fastq files exist
    if [[ -f "${sp}_rnaseq_norm_R1.fastq" ]] && [[ -f "${sp}_rnaseq_norm_R2.fastq" ]]; then
      echo "bbnorm skipped"
    elif [[ -f "${sp}_rnaseq_norm_R1.fastq.gz" ]] && [[ -f "${sp}_rnaseq_norm_R2.fastq.gz" ]]; then
      zcat "${sp}_rnaseq_norm_R1.fastq.gz" > "${sp}_rnaseq_norm_R1.fastq"
      zcat "${sp}_rnaseq_norm_R2.fastq.gz" > "${sp}_rnaseq_norm_R2.fastq"
      echo "bbnorm skipped (unzipped)"
    else
      mv */*fastq.gz ./
      zcat *_1.fastq.gz > "${sp}_rnaseq_R1.fastq"
      zcat *_2.fastq.gz > "${sp}_rnaseq_R2.fastq"

      # Use BBnorm to remove reads in regions with depth > maxCoverage
      bbnorm.sh in="${sp}_rnaseq_R1.fastq" in2="${sp}_rnaseq_R2.fastq" out="${sp}_rnaseq_norm_R1.fastq" out2="${sp}_rnaseq_norm_R2.fastq" ignorebadquality \
        qin=auto qout=auto -da min=2 target=100 threads="$threads" -Xmx"$mem"
    fi
    cd ..

    if [[ -f "${sp}/${sp}_rnaseq_norm_R1.fastq" ]] && [[ -f "${sp}/${sp}_rnaseq_norm_R2.fastq" ]]; then
      # STAR for mapping RNA-seq
      mkdir -p bams
      cd bams

      # Copy the assembly file
      cp "/data/home/s2215768/fly_annotation/data/CAT_genomes/selected_genomes/${assembly_name}" "${sp}.assembly.fasta"

      # Create necessary directories
      mkdir "${sp}_genome_fasta"
      mkdir "${sp}_genomeDir"
      mkdir "${sp}_aligned_reads"

      # Generate genome index
      STAR --runThreadN "$threads" --genomeDir "./${sp}_genomeDir/" --runMode genomeGenerate --limitGenomeGenerateRAM 150000000000 --genomeFastaFiles "${sp}.assembly.fasta"

      # Map reads
      STAR --runThreadN "$threads" --genomeDir "./${sp}_genomeDir/" --readFilesIn "../${sp}/${sp}_rnaseq_norm_R1.fastq" "../${sp}/${sp}_rnaseq_norm_R2.fastq" --outFileNamePrefix "./${sp}_aligned_reads/"

      # Convert SAM to BAM and sort
      samtools view -Sb -@"$threads" "./${sp}_aligned_reads/Aligned.out.sam" | sambamba sort -o "${hal}_rnaseq.bam" -p -t "$threads" -m "$mem" /dev/stdin && samtools index "${hal}_rnaseq.bam"

      if [[ -f "./${hal}_rnaseq.bam" ]]; then
        echo "$sp rnaseq read mapping successfully completed"
        cat "../${sp}/${sp}_rnaseq_norm_R1.fastq" | pigz -p "$threads" > "../${sp}/${sp}_rnaseq_norm_R1.fastq.gz"
        cat "../${sp}/${sp}_rnaseq_norm_R2.fastq" | pigz -p "$threads" > "../${sp}/${sp}_rnaseq_norm_R2.fastq.gz"
      else
        echo "$sp rnaseq read mapping failed"
      fi

      # Clean up temporary directories and files
      rm -r "${sp}_genomeDir" "${sp}_aligned_reads" "${sp}_genome_fasta" "${sp}.assembly.fasta"
      rm "../${sp}"/*.fastq
      cd ..
    fi
  fi
done < "$sp_list"

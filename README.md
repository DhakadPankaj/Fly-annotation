# Drosophila genome annotation pipeline

The provided scripts are used to annotate ~250 Drosophila genomes using Comparative Annotation Toolkit (CAT). Here is the original publication of CAT which provides some of the background necessary to understand how CAT works: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6028123/  and the github repo: https://github.com/ComparativeGenomicsToolkit/Comparative-Annotation-Toolkit

## General steps before runnig CAT
There are steps that will need to be performed before you can run the CAT pipeline.

### CAT installation 
To setup the conda environment for annotation and install CAT and all the dependencies follow the steps and instructions given in: ```CAT/cat_installation.sh```

### Prerequisite data for better gene prediction

#### 1. Obtain genome files and mask repeats
Here is a function to download genome assembly from NCBI provided genome accession and species name. Or you could use your self assembled genomes.

Usage: ```download_assembly_ftp species_name genome_accession```
```bash
download_assembly_ftp() {
    if [[ $2 = GCF_* ]]
    then
    esearch -db assembly -query $2 </dev/null | esummary| grep 'FtpPath type="RefSeq"'|sed 's/<FtpPath type="RefSeq">//'|sed 's_</FtpPath>__'|sed 's_\s__g' \
     |while read -r url ; do
    fname=$(echo $url | grep -o 'GCF_.*' | sed 's/$/_genomic.fna.gz/') ;
    wget -O "${1}.${2}_genomic.fna.gz"  "$url/$fname" ;
    done;
    else
    esearch -db assembly -query $2 </dev/null | esummary | xtract -pattern DocumentSummary -element FtpPath_GenBank \
     | while read -r url ; do
    fname=$(echo $url | grep -o 'GCA_.*' | sed 's/$/_genomic.fna.gz/') ;
    wget -O "${1}.${2}_genomic.fna.gz" "$url/$fname" ;
    done ;
    fi
}
```

*Optional: It is best to softmask repeats in the genome as they will affect the gene predictions. 

I've used repeatModeler2 to generate individualized repeat libraries for each genome. To setup the conda environment for repeat masking follow the instructions given in: ```cactus_alignment/setup_repeatmodeler.sh``` and to run the repeatModeler2/repeatMasker you can follow the instructions provided in this workshop: https://github.com/ISUgenomics/bioinformatics-workbook/blob/master/dataAnalysis/ComparativeGenomics/RepeatModeler_RepeatMasker.md 

#### 2. Get all trancriptional data available for your genomes (ESTs, RNA-seq, Transcripts, Isoseq)
#### 3. Get protein data


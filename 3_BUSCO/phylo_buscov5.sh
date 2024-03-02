#! /bin/bash

# This script builds a maximum likelihood phylogenetic tree using BUSCO genes.

# Set default values
threads=40
busco_supermatrix="OFF"
gene_tree="OFF"
species_tree="OFF"
batchscript="OFF"

# Parse command line arguments
while true ; do
    case "$1" in
        -a|--Supermatrix)
            busco_supermatrix=$2; shift 2;;
        -t|--cores)
            threads=$2; shift 2;;
        -g|--GeneTrees)
            gene_tree=$2; shift 2;;
        -s|--SpeciesTree)
            species_tree=$2; shift 2;;
        -b|--batchScript)
            batchscript=$2; shift 2;;
        --)
            shift; break;;
        *)
            echo "Usage: $0 [-t cores] [-a ON|OFF] [-g ON|OFF] [-s ON|OFF] [-b ON|OFF]"  # Usage instructions
            exit 1;;
    esac
done

echo "Running with $threads threads/cores"

# Generate supermatrix alignment from BUSCO output
if [ $busco_supermatrix == "ON" ]; then
    python ~/fly_annotation/tools/BUSCO_phylogenomics/BUSCO_phylogenomics.py -i BUSCO/ -o BUSCO_phylogenomics/ -t $threads --percent_single_copy 99 --supermatrix_only
elif [ $busco_supermatrix == "OFF" ]; then
    echo "Supermatrix: OFF"
else
    echo "Please check -a or --Supermatrix argument... EXITING"
    exit 1
fi

# Reconstruct gene trees with iqtree
if [ $gene_tree == "ON" ]; then
    mkdir -p iqtree
    treethreads=4
    ls -1 BUSCO_phylogenomics/supermatrix/trimmed_alignments/ | sed 's/.trimmed.aln//' | \
    parallel -j$(echo "${threads}/${treethreads}" | bc) \
    iqtree2 -s ./BUSCO_phylogenomics/supermatrix/trimmed_alignments/{}.trimmed.aln -bb 1000 -nt ${treethreads} -m LG+F+G \
    -blmin 1e-300 -safe -pre iqtree_{} -mem 18G

    mv iqtree_* ./iqtree

    # Concatenate IQTREE gene trees
    cat iqtree/*.treefile > busco_gene_trees_ml.tree
elif [ $gene_tree == "OFF" ]; then
    echo "GeneTree: OFF"
else
    echo "Please check -g or --GeneTrees argument... EXITING"
    exit 1
fi

# Reconstruct species tree with ASTRAL
if [ $species_tree == "ON" ]; then
    java -jar ~/fly_annotation/tools/Astral/astral.5.15.5.jar --input busco_gene_trees_ml.tree \
    --output busco_species_astral.tree -T $threads

    # ASTRAL species Tree rescaling
    raxml-ng --evaluate --threads $threads --model LG+F+G --tree busco_species_astral.tree \
    --msa BUSCO_phylogenomics/supermatrix/SUPERMATRIX.fasta --prefix rescaled \
    --opt-model on --opt-branches on
elif [ $species_tree == "OFF" ]; then
    echo "SpeciesTree: OFF"
else
    echo "Please check -s or --SpeciesTree argument... EXITING"
    exit 1
fi

# Generate batch script for Clubashworth IQTREE2
if [ $batchscript == "ON" ]; then
    cat << EOF > run_iqtree.sh
#!/bin/bash

#$ -V
#$ -cwd
#$ -N iqtree
#$ -l h="c1|c2|c3|c4|c5|c6"
#$ -t 1-1085
#$ -tc 50
#$ -pe smp 6
#$ -e std_out/
#$ -o std_out/

source /ceph/software/miniconda/etc/profile.d/conda.sh &&
conda activate fly_annotation &&

busco_id=\`ls -1 BUSCO_phylogenomics/supermatrix/trimmed_alignments/ | sed 's/.trimmed.aln//' | awk "NR==$SGE_TASK_ID"\`

mkdir -p iqtree
mkdir -p /scratch/$USER/iqtree &&

iqtree2 -s BUSCO_phylogenomics/supermatrix/trimmed_alignments/\${busco_id}.trimmed.aln -B 1000 -nt 6 -m LG+F+G \
    -blmin 1e-300 -safe -pre iqtree_\${busco_id} -mem 18G

mv iqtree_\${busco_id}* /scratch/$USER/iqtree/

rsync -av --remove-source-files /scratch/$USER/iqtree/ ./iqtree/

EOF

#qsub run_iqtree.sh

fi
############################################
#############################################



######
# Ultrametric time tree (Alt to BEAST)
#####
# Pre-requisite: fossil date file (Can be obtained using timetree websites)

#iqtree2 -s ALN_FILE --date dro_anc_dates.config -te TREE_FILE --date-tip 0 --date-root 65 --date-ci 100 -mem 18G


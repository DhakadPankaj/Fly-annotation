#! /bin/bash

# This script runs the OMArK pipeline on protein sets of Drosophila species.
# Arguments:
#   $1: species name
# Check if species list files is provided as argument
if [ -z "$1" ]; then
    echo "Please provide a species list file"
    exit 1
fi

run_OMArk() {
    # Run OMArK pipeline
    species=$(echo $1|tr [a-z] [A-Z])
    protein_path="/data/home/s2215768/fly_annotation/complement_annotations/orthology_inference"
    if [ ! -s OMArk_output/${species}/*.sum ]; then
        echo "Running OMArK for ${species}"
        mkdir -p OMArk_output/${species}
        #omamer search --db  OMAmer_Database/LUCA.h5 --query ${protein_path}/${species}.fa --out OMAmer_Database/${species}.omamer
        omark -f OMAmer_Database/${species}.omamer -d OMAmer_Database/LUCA.h5 -r family -t 7214 -o OMArk_output/${species}/
    else
        echo "OMArK for ${species} already exists"
        return
    fi    
}

# Run OMArK pipeline for each species
export -f run_OMArk
mkdir -p OMAmer_Database
mkdir -p OMArk_output

parallel -j 40 --colsep "\t" run_OMArk {1} :::: $1
#!/bin/bash

set -euo pipefail

# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate busco

# Set the number of threads
threads=3

run_busco() {
    sp=$1
    fasta=$3

    if [ -f "BUSCO/${sp}/run_diptera_odb10/short_summary.txt" ]; then
        echo "$sp busco short_summary.txt file exists"
    else
        busco -i "${fasta}" -c 15 -o "BUSCO/${sp}" -m geno -l ~/diptera_odb10 -f --offline
    fi

    cd BUSCO/
    busco=$(cat "${sp}/run_diptera_odb10/short_summary.txt" | grep -A2 "***** Results: *****" | tail -1 | awk 'BEGIN { FS=OFS=":" } {print $2}' | tr '%[S' ' ' | sed 's_\s__g')

    if (( $(echo "$busco >= 90" | bc -l) )); then
        echo -e "${1}\t${2}\t${3}\t${busco}" > "tmp/${sp}.tmp"
    else
        echo -e "${1}\t${2}\t${3}\t${busco}" > "../bad_BUSCO/${sp}.tmp"
        mv "${sp}/" "../bad_BUSCO/${sp}"
    fi

    cd ..
}

main() {
    if [ $# -eq 0 ]; then
        echo "Usage: $0 <species_table>"
        exit 1
    fi

    # Create necessary directories
    mkdir -p BUSCO/tmp
    mkdir -p bad_BUSCO

    # Export the run_busco function to be used by parallel
    export -f run_busco

    # Read species table and run busco in parallel
    sp_table=$1
    parallel -j "$threads" --colsep "\t" run_busco :::: "$sp_table"

    # Merge temporary files into final output files
    > sp_list_busco.tsv
    cat BUSCO/tmp/*.tmp >> sp_list_busco.tsv

    > sp_list_bad_busco.tsv
    cat bad_BUSCO/*.tmp >> sp_list_bad_busco.tsv

    # Clean up temporary files
    rm -r BUSCO/tmp/
}

main "$@"

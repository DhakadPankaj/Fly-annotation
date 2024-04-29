1. `run_braker3.sh`, is used to run the BRAKER3 gene annotation pipeline for multiple genomes. 

The script requires the following inputs:

1. A `genome_info` file that contains the species name and their corresponding assemblies. This file is passed as a command-line argument when running the script.
2. A protein sequence file named `Diptera_noDupl.faa` located in the same directory as the script.
3. RNA-seq data in the form of compressed fastq files located in the directory `/data/home/s2215768/fly_annotation/data/rnaseq_all`.

The script first checks if the necessary input files are provided. If not, it exits with an error message.

The `rna_seq` function is defined to handle the RNA-seq data. It checks if the RNA-seq data for a given species is available and not already processed. If so, it extracts the data, converts the compressed fastq files to uncompressed format, and removes the original compressed files.

The `run_braker` function is the main function that runs the BRAKER3 pipeline. It takes a species name and a genome assembly as input. It creates a new directory for each species and runs the BRAKER3 pipeline with the appropriate parameters. If RNA-seq data is available for the species, it is included in the BRAKER3 run. After the BRAKER3 run, the uncompressed fastq files are compressed again and the uncompressed versions are removed.

The script uses GNU `parallel` to run multiple instances of the `run_braker` function in parallel, allowing for efficient processing of multiple genomes.

Please note that this script assumes that the necessary software (BRAKER3, AUGUSTUS, etc.) and data files are already set up and available in the specified directories. It also assumes that the Conda environment `braker3` is set up with all the necessary dependencies.
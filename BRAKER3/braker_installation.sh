#!/bin/bash

conda create -n braker3 
conda activate braker3
conda install -c bioconda perl biopython perl-app-cpanminus perl-file-spec  perl-list-util perl-module-load-conditional \
 perl-posix perl-file-homedir perl-parallel-forkmanager perl-hash-merge\
  perl-scalar-util-numeric perl-yaml perl-class-data-inheritable perl-exception-class perl-test-pod perl-file-which perl-mce perl-threaded  \
  perl-list-util perl-math-utils cdbtools perl-data-dumper 
conda install -c eumetsat perl-yaml-xs


mkdir -p ${CONDA_PREFIX}/etc/conda/activate.d
wget https://github.com/Gaius-Augustus/BRAKER/archive/refs/tags/v3.0.8.tar.gz
tar -xzf v3.0.8.tar.gz
cd BRAKER-3.0.8
nano ${CONDA_PREFIX}/etc/conda/activate.d/env_vars.sh
#add the following lines
# export PATH="/data/home/s2215768/fly_annotation/tools_braker3/BRAKER-3.0.8/scripts:$PATH"

# install Augustus
git clone https://github.com/Gaius-Augustus/Augustus.git
cd Augustus
make augustus
# Error? install augustus dependencies from:  https://github.com/Gaius-Augustus/Augustus/blob/master/docs/INSTALL.md 


# GeneMark-ETP
git clone https://github.com/gatech-genemark/GeneMark-ETP.git

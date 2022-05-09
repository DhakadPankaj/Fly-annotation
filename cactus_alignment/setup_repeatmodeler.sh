#! /bin/bash

conda create --name repeatmasker python=3.8

conda activate repeatmasker

#setup tools directoru to setup everything
export TOOLDIR="/home/s2215768/fly_annotation_eddie/tools"
mkdir -p ${TOOLDIR} && cd ${TOOLDIR}

#install dependencies(if already not present)
conda install h5py
pip install -U pytest
conda install -c conda-forge pytest-mpi
conda install -c bioconda mafft
conda install -c bioconda cd-hit
#necessary perl modules
conda install -c bioconda perl-json perl-file-which perl-uri perl-devel-size perl-lwp-protocol-https

#install TRF
git clone https://github.com/Benson-Genomics-Lab/TRF.git
cd TRF/
mkdir build
cd build && ../configure
make
cp src/trf ${CONDA_PREFIX}/bin

#install RMblast
wget http://www.repeatmasker.org/rmblast-2.11.0+-x64-linux.tar.gz
tar zxvf rmblast-2.11.0+-x64-linux.tar.gz
cp bin/* ${CONDA_PREFIX}/bin

#install repeatmasker
wget http://www.repeatmasker.org/RepeatMasker/RepeatMasker-4.1.2-p1.tar.gz
tar xvzf RepeatMasker-4.1.2-p1.tar.gz
#Run Configure Script
perl ./configure


#install RepeatScout
wget http://www.repeatmasker.org/RepeatScout-1.0.6.tar.gz
tar xvzf RepeatScout-1.0.6.tar.gz && cd RepeatScout-1.0.6/
cp RepeatScout ${CONDA_PREFIX}/bin && cp build_lmer_table ${CONDA_PREFIX}/bin

#install RECON
wget http://www.repeatmasker.org/RepeatModeler/RECON-1.08.tar.gz
tar xvzf RECON-1.08.tar.gz
cd RECON-1.08/
cd src && make && make install
cd ..
###one more step
#In "recon.pl" in the scripts directory -- on the third line, add
#the path to the binaries (the bin directory here)

#path= /exports/eddie3_homes_local/s2215768/fly_annotation_eddie/tools/RECON-1.08/bin


#install LTR_harvest(a package in genometools)
wget http://genometools.org/pub/genometools-1.6.2.tar.gz
tar xvzf genometools-1.6.2.tar.gz && cd genometools-1.6.2/
make -j4 threads=yes
make -j4 install prefix=${TOOLDIR}/gt
#path /exports/eddie3_homes_local/s2215768/fly_annotation_eddie/tools/gt/

#install LTR_retriever
wget https://github.com/oushujun/LTR_retriever/archive/refs/tags/v2.9.0.tar.gz
tar xvzf v2.9.0.tar.gz
#path= /home/s2215768/fly_annotation_eddie/tools/LTR_retriever-2.9.0

#install NINJA
wget https://github.com/TravisWheelerLab/NINJA/archive/refs/tags/0.95-cluster_only.tar.gz
tar xvzf 0.95-cluster_only.tar.gz
cd NINJA-0.95-cluster_only/NINJA && make && cd ..
cp NINJA/* ${CONDA_PREFIX}/bin
#pwd /home/s2215768/fly_annotation_eddie/tools/NINJA-0.95-cluster_only/NINJA

#download UCSC twoBit tools
for tool in faToTwoBit twoBitDup twoBitInfo twoBitMask twoBitToFa; do
  wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/${tool}
  chmod +x ${tool}
  cp ${tool} ${CONDA_PREFIX}/bin
done

#install Repeatmodeler
wget https://www.repeatmasker.org/RepeatModeler/RepeatModeler-2.0.3.tar.gz
tar xvzf RepeatModeler-2.0.3.tar.gz
cd RepeatModeler-2.0.3/

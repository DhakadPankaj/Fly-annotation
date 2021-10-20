#Create a common directory to install all dependencies/tools

export TOOLDIR="${HOME}/fly_annotation/tools"

#create a conda environment(annotation) and install required CAT dependencies available on conda 
conda create -y -n annotation -c conda-forge -c bioconda -c defaults python=3.7 pyfasta luigi seaborn pandas ete3 pysam numpy scipy bx-python bcbio-gff biopython parasail-python configobj sqlalchemy samtools bamtools wiggletools bedtools phame
conda activate annotation

#If sudo access available
#If you are not root user build these from source and update 'INCLUDE' and 'LIBRARY' paths in common.mk file of augustus directory    
sudo apt install -y libboost-iostreams-dev zlib1g-dev libgsl-dev libboost-all-dev libsuitesparse-dev liblpsolve55-dev libsqlite3-dev libmysql++-dev libbamtools-dev libboost-all-dev libcurl4-openssl-dev libssl-dev libncurses5-dev libbz2-dev liblzma-dev

#manual installation of libssl1.0.0, required for running lower version of samtools(v1.5-v1.12)
wget -q http://security.ubuntu.com/ubuntu/pool/main/o/openssl1.0/libssl1.0.0_1.0.2n-1ubuntu5.3_amd64.deb
sudo dpkg -i libssl1.0.0_1.0.2n-1ubuntu5.3_amd64.deb && rm libssl1.0.0_1.0.2n-1ubuntu5.3_amd64.deb 

#Samtools, htslib and bcftools installation
for t in "htslib" "bcftools" "samtools"; do
    git clone https://github.com/samtools/${t}.git
    cd ${t}
    autoheader
    autoconf
    ./configure
    make -j
    cd ..
done


#download executable kent toolkit from UCSC
cd ${TOOLDIR}
for tool in faToTwoBit gff3ToGenePred genePredToBed genePredToFakePsl bamToPsl transMapPslToGenePred pslPosTarget axtChain chainMergeSort pslMap pslRecalcMatch pslMapPostChain gtfToGenePred genePredToGtf pslCDnaFilter 
clusterGenes pslToBigPsl bedSort bedToBigBed wigToBigWig; do
  wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/${tool}
  chmod +x ${tool}
  mv ${tool} ${CONDA_PREFIX}/bin
done


#build bamtools (especially for augustus)
cd ${TOOLDIR}
git clone https://github.com/pezmaster31/bamtools.git
cd bamtools
mkdir -p bamtools_build && cd bamtools_build
cmake -DCMAKE_INSTALL_PREFIX=${TOOLDIR}/bamtools/bamtools_install ..
make
make install

#Add path to common.mk of augustus directory

#Install augustus (If error encountered update paths in common.mk file)
cd cd ${TOOLDIR}
git clone https://github.com/Gaius-Augustus/Augustus.git
cd Augustus && make -j && AUGUSTUS_PATH=$(pwd)/bin:$(pwd)/scripts && cd ..


#install HDF5
wget http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.1/src/hdf5-1.10.1.tar.gz
tar zxf hdf5-1.10.1.tar.gz && rm hdf5-1.10.1.tar.gz && cd hdf5-1.10.1
./configure --enable-cxx --prefix ${TOOLDIR}/hdf5 && make -j $(nproc) && make install && cd .. && rm -r ./hdf5-1.10.1
export h5prefix=-prefix=${TOOLDIR}/hdf5
cp ${TOOLDIR}/hdf5/bin/* ${CONDA_PREFIX}/bin


#install cactus
git clone https://github.com/ComparativeGenomicsToolkit/cactus.git --recursive
cd cactus
pip install --upgrade setuptools pip
pip install --upgrade -r toil-requirement.txt
pip install --upgrade .
make -j $(nproc)
cp ./bin/* ${CONDA_PREFIX}/bin
cd ..

# install CLAPACK
wget http://www.netlib.org/clapack/clapack.tgz
tar xzf clapack.tgz && rm clapack.tgz && mv CLAPACK-3.2.1 clapack
cd clapack 
cp make.inc.example make.inc && make -j f2clib && make -j blaslib && make -j lib
export CLAPACKPATH=$(pwd)
cd ..


# install Phast
git clone https://github.com/CshlSiepelLab/phast.git
cd phast/src && make -j && cd ..
cp ./bin/* ${CONDA_PREFIX}/bin && cd ..
export ENABLE_PHYLOP=1



# install sonLib
git clone https://github.com/ComparativeGenomicsToolkit/sonLib.git
pushd sonLib && make -j && popd


# install HAL
git clone https://github.com/ComparativeGenomicsToolkit/hal.git 
cd hal && make -j
chmod -R u+w ${CONDA_PREFIX}/bin
cp ./bin/* ${CONDA_PREFIX}/bin
cd ..
cp -r ./hal ${CONDA_PREFIX}/lib/python3.7/site-packages/
pip3 install newick


# install sambamba
wget -q https://github.com/biod/sambamba/releases/download/v0.7.1/\
sambamba-0.7.1-linux-static.gz \
 && gunzip sambamba-0.7.1-linux-static.gz \
 && chmod +x sambamba-0.7.1-linux-static \
 && mv sambamba-0.7.1-linux-static ${CONDA_PREFIX}/bin/sambamba

# install exonerate
wget -q http://ftp.ebi.ac.uk/pub/software/vertebrategenomics/exonerate/\
exonerate-2.2.0-x86_64.tar.gz \
 && tar zxf exonerate-2.2.0-x86_64.tar.gz \
 && chmod +x exonerate-2.2.0-x86_64/bin/* \
 && mv exonerate-2.2.0-x86_64/bin/* ${CONDA_PREFIX}/bin \
 && rm -r ./exonerate-2.2.0-x86_64*

# install CAT
git clone https://github.com/ComparativeGenomicsToolkit/\
Comparative-Annotation-Toolkit.git
pip install -e Comparative-Annotation-Toolkit



# install Ragout
pip install networkx==2.2
git clone https://github.com/fenderglass/Ragout.git
cd Ragout
python setup.py build && python setup.py install
cd ..




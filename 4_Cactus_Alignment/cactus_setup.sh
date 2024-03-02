# /bin/bash


export TOOLDIR="~/fly_annotation/tools"
mkdir -p ${TOOLDIR} && cd ${TOOLDIR}


conda create -y -n cactus python=3.8

conda activate cactus

#install dependencies
for i in wigToBigWig faToTwoBit bedToBigBed bigBedToBed bedSort hgGcPercent; do
  wget -q http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/${i} -O ${CONDA_PREFIX}/bin/${i}
  chmod ugo+x ${CONDA_PREFIX}/bin/${i}; 
done

#install hdf5
wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.12/hdf5-1.12.2/src/hdf5-1.12.2.tar.gz
tar xvzf hdf5-1.12.2.tar.gz && rm hdf5-1.12.2.tar.gz
cd hdf5-1.12.2/
./configure --enable-cxx --prefix ${TOOLDIR}/hdf5
make -j8 && make install
cd .. && rm -r ./hdf5-1.12.2
export h5prefix=-prefix=${TOOLDIR}/hdf5
cp ${TOOLDIR}/hdf5/bin/* ${CONDA_PREFIX}/bin


#install cactus
wget https://github.com/ComparativeGenomicsToolkit/cactus/releases/download/v2.0.5/cactus-bin-v2.0.5.tar.gz
tar xvzf cactus-bin-v2.0.5.tar.gz && rm cactus-bin-v2.0.5.tar.gz
cd cactus-bin-v2.0.5/
python3 -m pip install -U setuptools pip==21.3.1
python3 -m pip install -U -r ./toil-requirement.txt
python3 -m pip install -U .
cp -r ./bin/* ${CONDA_PREFIX}/bin
cp -r $(pwd)/submodules/* ${CONDA_PREFIX}/lib/python3.8/site-packages/
CACTUS_BIN="$(pwd)/bin"


# save conda dir for activation scripts
wd=${CONDA_PREFIX}
conda deactivate 

mkdir -p ${wd}/etc/conda/activate.d
mkdir -p ${wd}/etc/conda/deactivate.d

(
  echo "#! /bin/bash"
  echo "PATH=\${PATH}:${CACTUS_BIN}"
) > ${wd}/etc/conda/activate.d/env_vars.sh

(
  echo "#! /bin/bash"
  echo "export PATH=\${PATH%:*}"
) > ${wd}/etc/conda/deactivate.d/env_vars.sh

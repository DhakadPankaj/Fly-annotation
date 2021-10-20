#export TOOLDIR="${HOME}/fly_annotation/tools"

#conda create -y -n annotation -c conda-forge -c bioconda -c defaults python=3.7 pyfasta luigi seaborn pandas ete3 pysam numpy scipy bx-python bcbio-gff biopython parasail-python configobj sqlalchemy samtools bamtools wiggletools bedtools phame

# install samtools stuff
#mkdir -p ${TOOLDIR} && cd ${TOOLDIR}
#for t in "htslib" "bcftools" "samtools"; do
#    git clone https://github.com/samtools/${t}.git
#    cd ${t}
#    autoheader
#    autoconf
#    ./configure
#    make -j
#    cd ..
#done

#download executable tools from UCSC
cd tools/
for tool in faToTwoBit gff3ToGenePred genePredToBed genePredToFakePsl bamToPsl \
transMapPslToGenePred pslPosTarget axtChain chainMergeSort pslMap \
pslRecalcMatch pslMapPostChain gtfToGenePred genePredToGtf pslCDnaFilter \
clusterGenes pslToBigPsl bedSort bedToBigBed wigToBigWig; do
  wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/${tool}
  chmod +x ${tool}
  mv ${tool} ${CONDA_PREFIX}/bin
done

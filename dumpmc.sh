#!/bin/sh

VERSION=$1

SEEDS=$2

OUTDIR=${3:-10m}

mcfile=/home/graham/BeamSpectrumGenerator/${OUTDIR}/testbcg-${VERSION}-${SEEDS}.dat

echo 'Using mcfile '${mcfile}

# Make symbolic link to the MC file
if [ -L "mcfile.ini" ]
then
   rm mcfile.ini 
   ln -s ${mcfile} mcfile.ini
else
   ln -s ${mcfile} mcfile.ini
fi

# Compile code
gfortran -o dumpmc dumpmc.f

# Execute code
time ./dumpmc

# Rename shorter output file
mv fort.23 testbcg-${VERSION}-${SEEDS}.cdat
mv testbcg-${VERSION}-${SEEDS}.cdat $OUTDIR

# Tidy up
rm mcfile.ini

exit

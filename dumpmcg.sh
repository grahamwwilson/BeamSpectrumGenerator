#!/bin/sh

VERSION=$1

SEEDS=$2

OUTDIR=${3:-1m}

mcfile=/home/graham/BeamSpectrumGenerator/${OUTDIR}/testbc-${VERSION}-${SEEDS}.dat

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
mv fort.23 testbc-${VERSION}-${SEEDS}.cdat
mv testbc-${VERSION}-${SEEDS}.cdat $OUTDIR

# Tidy up
rm mcfile.ini

exit

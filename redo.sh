#!/bin/sh

VERSION=$1

SEEDS=$2

echo 'VERSION: '${VERSION}
echo 'SEEDS:   '${SEEDS}

rm testbc.hbook

# link input beamstrahlung parameters file
rm betapars.f
ln -s betapars${VERSION}.f betapars.f

# link random number seeds file
rm seeds.f
ln -s seeds${SEEDS}.f seeds.f

# recompile picking up the defined random number seeds and parameters
./cl.sh

# execute
time ./testbc

# Rename output file
mv fort.45 testbc-${VERSION}-${SEEDS}.dat
cp testbc.hbook testbc-${VERSION}-${SEEDS}.hbook

exit

#!/bin/sh

VERSION=$1

SEEDS=$2

OUTDIR=${3:-1m}

echo 'VERSION: '${VERSION}
echo 'SEEDS:   '${SEEDS}
echo 'OUTDIR:  '${OUTDIR}

rm testbcg.hbook

# link input beamstrahlung parameters file
rm betapars.f
ln -s betapars${VERSION}.f betapars.f

# link random number seeds file
rm seeds.f
ln -s seeds${SEEDS}.f seeds.f

# recompile picking up the defined random number seeds and parameters
./cl.sh testbcg

# execute
time ./testbcg

# Rename output file
mv fort.45 testbcg-${VERSION}-${SEEDS}.dat
cp testbcg.hbook testbcg-${VERSION}-${SEEDS}.hbook

mv testbcg-${VERSION}-${SEEDS}.dat $OUTDIR
mv testbcg-${VERSION}-${SEEDS}.hbook $OUTDIR

exit

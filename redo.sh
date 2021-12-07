#!/bin/sh

VERSION=$1

# link random number seeds file
rm seeds.f
ln -s seeds${VERSION}.f seeds.f

# recompile
./cl.sh

# execute
time ./testbc

# Rename output file
mv fort.41 testbc-${VERSION}.dat

exit

#!/bin/sh

VERSION=${1:-1-0}

INFILE=testbc-${VERSION}.dat

echo 'Input file '$INFILE

date
tail -n +17 testbc-${VERSION}.dat | sort -g --key=3 >testbc-${VERSION}-1M.x1sorted
date
tail -n +17 testbc-${VERSION}.dat | sort -g --key=4 >testbc-${VERSION}-1M.x2sorted
date

exit

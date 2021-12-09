#!/bin/sh

OUTDIR=$1

./redod.sh 1 0 $OUTDIR 
./redod.sh 1 1 $OUTDIR
./redod.sh 2 2 $OUTDIR
./redod.sh 2 3 $OUTDIR

exit

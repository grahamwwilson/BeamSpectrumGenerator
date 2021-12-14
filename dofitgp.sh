#!/bin/sh

date

RW=$1

./reweightfitgp.sh 0 ${RW}
./reweightfitgp.sh 1 ${RW}
./reweightfitgp.sh 2 ${RW}
./reweightfitgp.sh 3 ${RW}
./reweightfitgp.sh 4 ${RW}

./reweightfitgp.sh 5 ${RW}
./reweightfitgp.sh 6 ${RW}
./reweightfitgp.sh 7 ${RW}
./reweightfitgp.sh 8 ${RW}
./reweightfitgp.sh 9 ${RW}

date

exit

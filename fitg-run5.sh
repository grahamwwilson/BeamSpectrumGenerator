#!/bin/sh

date

RW=$1

./reweightfitg.sh 0 ${RW}
./reweightfitg.sh 1 ${RW}
./reweightfitg.sh 2 ${RW}
./reweightfitg.sh 3 ${RW}
./reweightfitg.sh 4 ${RW}

date

exit

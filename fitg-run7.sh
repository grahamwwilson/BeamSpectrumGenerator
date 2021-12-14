#!/bin/sh

date

RW=$1

./reweightfitg.sh 5 ${RW}
./reweightfitg.sh 6 ${RW}
./reweightfitg.sh 7 ${RW}
./reweightfitg.sh 8 ${RW}
./reweightfitg.sh 9 ${RW}

date

exit

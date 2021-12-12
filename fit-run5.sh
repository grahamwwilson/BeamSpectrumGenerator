#!/bin/sh

date

RW=2

./reweightfit.sh 0 ${RW}
./reweightfit.sh 1 ${RW}
./reweightfit.sh 2 ${RW}
./reweightfit.sh 3 ${RW}
./reweightfit.sh 4 ${RW}

date

exit

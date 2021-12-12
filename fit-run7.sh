#!/bin/sh

date

RW=2

./reweightfit.sh 5 ${RW}
./reweightfit.sh 6 ${RW}
./reweightfit.sh 7 ${RW}
./reweightfit.sh 8 ${RW}
./reweightfit.sh 9 ${RW}

date

exit

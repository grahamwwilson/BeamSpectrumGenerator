#!/bin/sh

I=$1
J=$2
VERSION=1

echo 'ichoice (data data-set)  = '$I
echo 'jchoice (RW MC data-set) = '$J

rm ichoice.f
ln -s inc/ichoice${I}.f ichoice.f

rm jchoice.f
ln -s inc/jchoice${J}.f jchoice.f

./cl.sh rwbinl

time ./rwbinl >reweightfit-${I}-${J}-${VERSION}.out

exit

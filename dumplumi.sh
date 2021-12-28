#!/bin/sh
#
# Read in GP lumi.ee.out file, analyze and dump contents to .cdat file
#

RUN=$1

OUTDIR=/home/graham/gpDigest

lumifile=/home/graham/gpDigest/lumi-Run${RUN}.ee.out
echo 'Using lumifile '${lumifile}
wc -l ${lumifile} >LumiFile_LineCount.dat

# Make symbolic link to the lumi file
if [ -L "lumifile.ini" ]
then
   rm lumifile.ini 
   ln -s ${lumifile} lumifile.ini
else
   ln -s ${lumifile} lumifile.ini
fi

# execute precompiled code
time ./dumplumi

# Rename shorter output file
OUTFILE=gplumi-Run${RUN}-ALL.cdat
mv fort.23 ${OUTFILE}
mv ${OUTFILE} $OUTDIR

# Tidy up
rm lumifile.ini
rm LumiFile_LineCount.dat

exit

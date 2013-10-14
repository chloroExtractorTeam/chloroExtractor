#!/bin/bash

# Author: Thomas Hackl, thomas.hackl@uni-wuerzburg.de

# reliable detect binary folder
pushd `dirname $0` > /dev/null
BD=`pwd`
popd > /dev/null

HASH=$1;
READS=$2;
THREADS=$3;
CHUNKS=$4;
FIRST_CHUNK=$5;

if [ -z $CHUNKS ] || [ $CHUNKS -lt $THREADS ] ; then
        CHUNKS=$THREADS
fi;

if [ -z $FIRST_CHUNK ] ; then
        FIRST_CHUNK=1
fi;

PRE=`basename "$READS"`;
PRE=${PRE%.fast[aq]}
PRE=${PRE%.f[aq]}
SUF="mer.cov"

(for i in $(seq $FIRST_CHUNK $CHUNKS); do
        echo '"echo chunk '"$i/$CHUNKS"'; SeqChunker -n '"$CHUNKS"' -f '"$i"' -l '"$i $READS"' | '"$BD"'/read_cov.pl --no-verbose --hash '"$HASH"' > '"$PRE.$SUF.$i"'" ';
done;) | xargs -P $THREADS -n1 bash -c;

pv $(ls -1 *."$SUF".* | sort -V) > "$PRE.$SUF"
rm *."$SUF".*;


#!/bin/bash

if [ "$#" -ne 5 ]; then
    echo "usage: $0 <cram> <in_bed> <ref> <out> <chrom>"
    exit -1
fi

OUT=$4.$5

grep -E "^$5" $2 > $OUT.i.bed

./cnv_realign.py $1 $OUT.i.bed $3 $OUT

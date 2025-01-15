#!/bin/bash

if [ "$#" -ne 4 ]; then
    echo "usage: $0 <cram> <in_bed> <ref> <out>"
    exit -1
fi

cut -f 1 <$2 | sort | uniq | parallel ./cnv_realign_chrom.sh $1 $2 $3 $4
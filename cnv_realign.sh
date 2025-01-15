#!/bin/bash

if [ "$#" -ne 4 ]; then
    echo "usage: $0 <cram> <in_bed> <ref> <out>"
    exit -1
fi


jupyter-nbconvert --to python cnv_realign.ipynb --stdout | python /dev/stdin $1 $2 $3 $4

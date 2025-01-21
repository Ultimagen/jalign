#!/bin/bash 

rm -rf out/results
mkdir out/results

echo Build aggregated bed file ...
cat out/*bed >out/results/output.bed

echo Build merged bam file ...
samtools view -H `find out -name "*.bam" | head -1` >out/results/output.sam
cat out/*sam | grep -v -E "^@" >>out/results/output.sam
samtools sort --threads 64 out/results/output.sam >out/results/output.bam
samtools index out/results/output.bam
rm out/results/output.sam


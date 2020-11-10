#!/bin/bash

cd /home/test/chip-seq/
mkdir output

## Peak Calling Using HOMER
# convert bam file into tag file using makeTagDirectory
makeTagDirectory input/ip input/ip.part.bam
makeTagDirectory input/input input/input.part.bam

# use findPeaks to call peak
findPeaks input/ip/ -style factor -o output/part.peak -i input/input/
# -style: factor (transcription factor ChIP-Seq), histone (histone modification ChIP-Seq), and so on

## Motif Analysis
# find enriched motifs in ChIP-Seq peaks
findMotifsGenome.pl output/part.peak sacCer2 output/part.motif.output -len 8

## peak calling using MACS
macs2 callpeak -t input/ip.part.bam -c input/input.part.bam --outdir output/macs_peak --name=yeast_macs_p05 --format=BAM --gsize=1.2e7 --tsize=50 --pvalue=1e-5


## homework

cd /home/test/chip-seq/
mkdir homework-output

makeTagDirectory homework/ip homework/ip.chrom_part.bam
makeTagDirectory homework/input homework/input.chrom_part.bam

findPeaks homework/ip/ -style factor -o homework-output/chrom_part.peak -i homework/input/ -F 8 -P 1e-8    # fold change>8 and p-value<1e-8 (vs control)

findMotifsGenome.pl homework-output/chrom_part.peak sacCer2 homework-output/chrom_part.motif.outout -len 8    # need p-value<1e-10 but no option found; manual select using awk is needed



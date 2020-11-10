#!/bin/bash

### Prepare bam files from ChIP-seq data

## download data
# The fastq data for yeast ChIP-seq was downloaded from GSE61210.
# Input data was downloaded from GSM1499619;
# IP data was downloaded from GSM1499607.

##　build yeast bowtie index
# Yeast sacCer2 genome data was downloaded from UCSC http://hgdownload.soe.ucsc.edu/goldenPath/sacCer2/bigZips/chromFa.tar.gz.
# Index was build with commad:
tar -xvf chromfa.tar.gz
cat *.fa >yeast.allchrom.fa
mkdir bowtie_index_yeast
bowtie-build yeast.allchrom.fa bowtie_index_yeast/sacCer2

## mapping
bowtie -p 4  -m 1  -v 3  --best --strata bowtie_index_yeast/sacCer2 \
    -q input/ip.fastq -S input/ip.sam

## sampling
# As the .sam file is too big for tutorial example, so we selected parts of them as example file.
samtools sort input/ip.sam >input/ip.sorted.bam
samtools index input/ip.sorted.bam
samtools view input/ip.sorted.bam chrI chrII chrIII -b >input/ip.part.bam

#!/bin/bash

### How to generate bam file for alternative splicing

## install hisat(ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip), bamtools, samtools

## get genome
# chromFa.tar.gz(https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/chromFa.tar.gz)
# Musmusculus.GRCm38.93.gtf.gz(ftp://ftp.ensembl.org/pub/release-93/gtf/mus_musculus/Mus_musculus.GRCm38.93.gtf.gz)

## get raw data
# ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR065/SRR065544/SRR065544_1.fastq.gz
# ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR065/SRR065544/SRR065544_2.fastq.gz
# ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR065/SRR065545/SRR065545_1.fastq.gz
# ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR065/SRR065545/SRR065545_2.fastq.gz

## make hisat index
# extract X chromosome sequence
tar -xz -f chromFa.tar.gz chrX.fa
mv chrX.fa Mus_musculus_chrX.fa

# use only X chromosome
zcat Mus_musculus.GRCm38.93.gtf.gz | grep -P '(#!)|(X\t)' > Mus_musculus_chrX.gtf

# make hisat index
hisat2_extract_splice_sites.py Mus_musculus_chrX.gtf > Mus_musculus_chrX.ss
hisat2_extract_exons.py Mus_musculus_chrX.gtf > Mus_musculus_chrX.exon

mkdir hisat2_indexes
hisat2-2.1.0/hisat2-build -p 4 --ss Mus_musculus_chrX.ss --exon Mus_musculus_chrX.exon Mus_musculus_chrX.fa hisat2_indexes/Mus_musculus_chrX

## mapping
# mapping
hisat2-2.1.0/hisat2 -p 4 --dta \
    -S SRR065544_chrX.sam -x hisat2_indexes/Mus_musculus_chrX \
    -1 SRR065544_1.fastq.gz -2 SRR065544_2.fastq.gz
hisat2-2.1.0/hisat2 -p 4 --dta \
    -S SRR065545_chrX.sam -x hisat2_indexes/Mus_musculus_chrX \
    -1 SRR065545_1.fastq.gz -2 SRR065545_2.fastq.gz

# covert to .bam
samtools sort -@ 4 -o SRR065544_chrX_raw.bam SRR065544_chrX.sam
samtools sort -@ 4 -o SRR065545_chrX_raw.bam SRR065545_chrX.sam

# filter only mapped reads
bamtools index -in SRR065544_chrX_raw.bam
bamtools index -in SRR065545_chrX_raw.bam

bamtools filter -isMapped true -in SRR065544_chrX_raw.bam -out SRR065544_chrX.bam
bamtools filter -isMapped true -in SRR065545_chrX_raw.bam -out SRR065545_chrX.bam



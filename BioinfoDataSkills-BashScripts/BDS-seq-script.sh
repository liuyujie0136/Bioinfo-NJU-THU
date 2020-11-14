#!/bin/bash

### Working with Sequence Data - scripts in book Bioinformatics Data Skills (chapter 10)


## count entries in fastq files
# @ can be quality, wc -l can be wrong for not-4-line entries
bioawk -cfastx 'END{print NR}' untreated1_chr4.fq


## quality scores in fastq files (Sanger version) in python
for b in quality:
    q=10**(-(ord(b)-33)/10)


## Inspecting and Trimming Low-Quality Bases
# Java program FastQC; R package qrqc; linux tools sickle and seqtk

sickle se -f untreated1_chr4.fq -t sanger -o untreated1_chr4_sickle.fq
seqtk trimfq untreated1_chr4.fq > untreated1_chr4_trimfq.fq 

# in R
# trim_qual.R -- explore base qualities before and after trimming
library(qrqc)
# FASTQ files
fqfiles <- c(none="untreated1_chr4.fq", sickle="untreated1_chr4_sickle.fq", trimfq="untreated1_chr4_trimfq.fq")
# Load each file in, using qrqc's readSeqFile 
# We only need qualities, so we turn off some of readSeqFile's other features. 
seq_info <- lapply(fqfiles, function(file) {readSeqFile(file, hash=FALSE, kmer=FALSE)})
# Extract the qualities as dataframe, and append a column of which trimmer (or none) was used. This is used in later plots.
quals <- mapply(function(sfq, name) {qs <- getQual(sfq); qs$trimmer <- name; qs}, seq_info, names(fqfiles), SIMPLIFY=FALSE)
# Combine separate dataframes in a list into single dataframe
d <- do.call(rbind, quals)
# Visualize qualities
p1 <- ggplot(d) + geom_line(aes(x=position, y=mean, linetype=trimmer))
p1 <- p1 + ylab("mean quality (sanger)") + theme_bw()
print(p1)
# Use qrqc's qualPlot with list produces panel plots. Only shows 10% to 90% quantiles and lowess curve
p2 <- qualPlot(seq_info, quartile.color=NULL, mean.color=NULL) + theme_bw()
p2 <- p2 + scale_y_continuous("quality (sanger)")
print(p2)


## Counting Nucleotides in python (biopython)
# nuccount.py -- tally nucleotides in a file
import sys
from collections import Counter
from readfq import readfq
IUPAC_BASES = "ACGTRYSWKMBDHVN-."
# Create a new Counter object
counts = Counter()
# for each sequence entry, add all its bases to the counter
for name, seq, qual in readfq(sys.stdin):
    counts.update(seq.upper())
# print the results
for base in IUPAC_BASES:
    print base + "\t" + str(counts[base])


## Index FASTA Files - using samtools
# create an index file named ~.fai
samtools faidx Mus_musculus.GRCm38.75.dna.chromosome.8.fa
# access a particular region
samtools faidx Mus_musculus.GRCm38.75.dna.chromosome.8.fa 8:123407082-123410744
# multiple regions
samtools faidx Mus_musculus.GRCm38.75.dna.chromosome.8.fa 8:123407082-123410744 8:123518835-123536649


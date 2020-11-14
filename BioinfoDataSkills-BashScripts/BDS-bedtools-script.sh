#!/bin/bash

## BEDTools for Ranges Data

bedtools intersect -a ranges-qry.bed -b ranges-sbj.bed
# return only overlapping regions
# like the GenomicRanges intersect() function

bedtools intersect -a ranges-qry.bed -b ranges-sbj.bed -wa
# -wa: the ranges in A(qry) that overlap B(sbj); -wb: reversed; -u: unique items in results

bedtools intersect -a query-sorted.bed -b subject-sorted.bed --sorted
# faster if sort A and B ahead

bedtools intersect -a ranges-qry.bed -b ranges-sbj.bed -wo
# return the number of overlapping bases; -v: all nonoverlapping ranges; -s: on the same strand

bedtools intersect  # for more help


bedtools slop -i ranges-qry.bed -g genome.txt -b 4
#  grow both sides of each range by 4 base pairs
#  like + or the resize() in GenomicRanges

bedtools slop -i ranges-qry.bed -g genome.txt -l 3 -r 5
# left (-l) 3 and right (-r) 5


bedtools flank -i mm_GRCm38.75_protein_coding_genes.gtf -g Mus_musculus.GRCm38_genome.txt -l 3000 -r 0 > mm_GRCm38_3kb_promoters.gtf 
# like flank(), slide windows upstream to find promoters


bedtools getfasta -fi Mus_musculus.GRCm38.75.dna_rm.toplevel_chr1.fa -bed mm_GRCm38_3kb_promoters.gtf -fo mm_GRCm38_3kb_promoters.fasta
# extract sequences for a given set of ranges


bedtools genomecov -i ranges.sorted.bed -g cov.txt
# summarizing the coverage of features along chromosome sequences (sort first)

bedtools genomecov -i ranges-cov.bed -g cov.txt -d
# per-base pair coverage, like GenomicRanges coverage()
# results: chromosome, the position on that chromosome, and the coverage at that position

bedtools genomecov -i ranges-cov.bed -g cov.txt -bg
# BedGraph format, more compact


## Other BEDTools Subcommands
bedtools annotate
# annotates how much coverage each of these files has over another input file
bedtools merge
# like GenomicRanges reduce(); merges overlapping ranges into a single range
bedtools closest
# like GenomicRanges nearest()
bedtools complement
# This is the BEDTools version of gaps() (similar to setdiff() too)
bedtools multicov
# counts the number of alignments in multiple BAM files that overlap a specified BED file. 
bedtools multiinter
# similar to intersect, but works with multiple file inputs
bedtools unionbedg
# merges multiple BedGraph files into a single file



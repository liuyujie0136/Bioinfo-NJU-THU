#!/bin/bash

### Working with Alignment Data - scripts in book Bioinformatics Data Skills (chapter 11)
### Processing SAM/BAM files (alignment/mapping) - using samtools

## SAM file sections - other ref for detail

## view header
samtools view -H celegans.bam | grep "^@RG"
samtools view -H celegans.bam # bam is also OK!


## view entire alignment section without the header (without -H)
samtools view celegans.sam


## Convert between SAM and BAM
samtools view -b celegans.sam > celegans_copy.bam
samtools view celegans.bam > celegans_copy.sam # without header!
samtools view -h celegans.bam > celegans_copy.sam # include header


## Samtools Sort 
samtools sort celegans_unsorted.bam celegans_sorted # the second argument is the output filename prefix (samtools sort will append the .bam extension for you)
-m 4G -@ 2 # -m: menmory; -@: threads


## Samtools Index
# Note: The BAM file must be sorted first, and we cannot index SAM files.
samtools index celegans_sorted.bam


## Samtools View
# usage: samtools view [options] <in.bam>|<in.sam>|<in.cram> [region ...]
# Extracting alignments from a region
samtools view NA12891_CEU_indexed.bam 1:215906469-215906652 > USH2A_sample_alns.bam
samtools view -L USH2A_exons.bed NA12891_CEU_indexed.bam # -L: many regions stored in a BED file

# Filtering alignments - according to bitwise flags, for example
samtools flags unmap # output: 0x4 4 UNMAP
samtools view -f 4 NA12891_CEU_indexed.bam
samtools view -F 4 NA12891_CEU_indexed.bam # -F: do not have any of the bits set of the supplied flag argument
samtools flags READ1,PROPER_PAIR # output: 0x42 66 PROPER_PAIR,READ1
samtools view -f 66 NA12891_CEU_indexed.bam


## Visualizing Alignments with samtools tview
# view the very beginning of a chromosome
samtools tview NA12891_CEU_indexed.bam human_g1k_v37.fasta
# -p: go to a specific region
samtools tview -p 1:215906469-215906652 NA12891_CEU_indexed.bam human_g1k_v37.fasta
# BETTER: Integrated Genomics Viewer (IGV)

## Pileups(count tracks), Variant Calling, and Base Alignment Quality(BAQ)
samtools mpileup --no-BAQ --region 1:215906528-215906567 --fasta-ref human_g1k_v37.fasta NA12891_CEU_sample.bam

# output vcf(Variant Call Format) using -v, bcf(binary~) using -g
samtools mpileup -v --no-BAQ --region 1:215906528-215906567 --fasta-ref human_g1k_v37.fasta NA12891_CEU_sample.bam > NA12891_CEU_sample.vcf.gz

# then use bfctools to call real (significant) variant sites by -v
bcftools call -v -m NA12891_CEU_sample.vcf.gz > NA12891_CEU_sample_calls.vcf.gz

# bcftools call with all sites
bcftools call -m NA12891_CEU_sample.vcf.gz | grep -v "^##" | awk 'BEGIN{OFS="\t"} {split($8, a, ";"); print $1,$2,$4,$5,$6,a[1],$9,$10}'

# enables Base Alignment Quality (BAQ) in samtools mpileup
samtools mpileup -u -v -r 1:215906528-215906567 -f human_g1k_v37.fasta NA12891_CEU_sample.bam > NA12891_CEU_sample_baq.vcf # -u: uncompressed vfc file; -r: --region; -f: --fasta-ref


## Creating Your Own SAM/BAM Processing Tools with Pysam (skipped)

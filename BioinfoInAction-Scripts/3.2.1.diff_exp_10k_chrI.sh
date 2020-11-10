#!/bin/bash

cd /home/test/diff-exp/
mkdir mapping

for i in wt1 wt2 wt1X wt2X
do
	tophat -p 4 -G yeast_annotation.gff --no-coverage-search -o mapping/"$i"_thout bowtie_index/YeastGenome Raw_reads_10k/"$i".fq  # map reads using tophat
	bamtools index -in mapping/"$i"_thout/accepted_hits.bam  # index bam files
	bamtools filter -region chrI -in mapping/"$i"_thout/accepted_hits.bam -out mapping/"$i"_thout/chrI.bam  ## Extract mapped reads on chr I only (for quick running)
	cufflinks -p 4 -o assembly/"$i"_clout  mapping/"$i"_thout/chrI.bam  # Assemble transcripts for each sample
done

# Merge all assemblies to one file containing merged transcripts
ls assembly/*/transcripts.gtf > assembly/assemblies.txt
cuffmerge -g yeast_chrI_annotation.gff -s bowtie_index/YeastGenome.fa -p 4 -o assembly/merged assembly/assemblies.txt

# Identify differentially expressed genes and transcripts
cuffdiff -o diff_expr -b bowtie_index/YeastGenome.fa -p 4 -u assembly/merged/merged.gtf mapping/wt1_thout/chrI.bam,mapping/wt2_thout/chrI.bam mapping/wt1X_thout/chrI.bam,mapping/wt2X_thout/chrI.bam

# Differentially expressed genes output
awk 'BEGIN{print "gene_id\tgene_name\tlog2(fold_change)\tp_value\tq_value"}; ($10>1||$10<-1)&&$13<0.05{print $2 "\t" $3 "\t" $10 "\t" $12 "\t" $13}' diff_expr/gene_exp.diff > diff_exp_10k_chrI.out

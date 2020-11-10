#!/bin/bash

cd /home/test/diff-exp/
mkdir mapping_all

for i in wt1 wt2 wt1X wt2X
do
	tophat -p 4 -G yeast_annotation.gff --no-coverage-search -o mapping_all/"$i"_thout bowtie_index/YeastGenome Raw_reads_10k/"$i".fq  # map reads using tophat
	cufflinks -p 4 -o assembly_all/"$i"_clout  mapping_all/"$i"_thout/accepted_hits.bam  # Assemble transcripts for each sample
done

# Merge all assemblies to one file containing merged transcripts
ls assembly_all/*/transcripts.gtf > assembly_all/assemblies.txt
cuffmerge -g yeast_annotation.gff -s bowtie_index/YeastGenome.fa -p 4 -o assembly_all/merged assembly_all/assemblies.txt

# Identify differentially expressed genes and transcripts
cuffdiff -o diff_expr_all -b bowtie_index/YeastGenome.fa -p 4 -u assembly_all/merged/merged.gtf mapping_all/wt1_thout/accepted_hits.bam,mapping_all/wt2_thout/accepted_hits.bam mapping_all/wt1X_thout/accepted_hits.bam,mapping_all/wt2X_thout/accepted_hits.bam

# Differentially expressed genes output
awk 'BEGIN{print "gene_id\tgene_name\tlog2(fold_change)\tp_value\tq_value"}; ($10>1||$10<-1)&&$13<0.05{print $2 "\t" $3 "\t" $10 "\t" $12 "\t" $13}' diff_expr_all/gene_exp.diff > diff_exp_10k_all.out



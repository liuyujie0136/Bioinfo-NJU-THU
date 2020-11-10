#!/bin/bash

# Selected Linux data tools' scripts in book Bioinformatics Data Skills (chapter7).

cat 1.gtf | grep -E -o 'gene_id "\w+"' | cut -d ' ' -f 2 | sed 's/"//g' | sort | uniq -c

bioawk -c gff '{print $feature "\t" $start "\t" $end}' 1.gtf | head

tail -n +2 text.txt
# with a + sign (e.g., +x), tail will start from the xth line. 

# write something such as alias c="clear" into ~/.bashrc, and then "source ~/.bashrc".

less 1.gtf
# in less use command "/<pattern>" to search down (forward) for string <pattern> and "?<pattern>" to search up (backward) for string <pattern>

awk -F "\t" '{print NF; exit}' 1.gtf

grep -v "^#" 1.gtf | head -n 3

grep -v "^#" 1.gtf | cut -f 1-8 | column -t  | head -n 3
# only use columnt -t to visualize data in the terminal, not to reformat data to write to a file

grep -v -w bioinfo example.txt
# constrain matches to be words

grep -B1 "AGATCGG" contam.fastq | head -n 6
# get around this context before (-B), context: after (-A), and context before and after (-C)

grep -o "Olfr.*" Mus_musculus.GRCm38.75_chr1_genes.txt | head -n 3
# extract only the matching part of the pattern with -o

grep -E -o 'gene_id "\w+"' 1.gtf | head -n 5
# -E: extended regular expressions

sort -k1V -k2n example2.bed
# sort by the number in the text with -V

join -1 1 -2 1 example_sorted.bed example_lengths.txt
# usage: join -1 <file_1_field> -2 <file_2_field> <file_1> <file_2>, where <file_1_field> and  <file_2_field> means which coloum that their content should be same after join

join -1 1 -2 1 -a 1 example_sorted.bed example_lengths_alt.txt
# use -a to specify which file is allowed to have unpairable entries

awk '$1 ~ /chr1/ && $3 - $2 > 10' example.bed
# not match: $1 !~ /chr1/ or !($1 ~ /chr1/)

bioawk -c help
bioawk -c gff '$3 ~ /gene/ && $2 ~ /protein_coding/ {print $seqname,$end-$start}' 1.gtf | head -n 4 
bioawk -c fastx '{print ">"$name"\n"revcomp($seq)}' contam.fastq | head -n 4
bioawk -c fastx '{print $name,length($seq)}' Mus_musculus.GRCm38.75.dna_rm.toplevel.fa.gz

echo "chr1:28427874-28425431" | sed -E 's/^(chr[^:]+):([0-9]+)-([0-9]+)/\1\t\2\t\3/'
echo "chr1:28427874-28425431" | sed 's/[:-]/\t/g'
echo "chr1:28427874-28425431" | sed -e 's/:/\t/' -e 's/-/\t/'
echo "chr1:28427874-28425431" | tr ':-' '\t'

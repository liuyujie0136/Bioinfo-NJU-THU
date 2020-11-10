#!/bin/bash

## Use image 'bioinfo_motif' in this section
# docker load -i bioinfo_motif.tar.gz
​# docker run -dt --name motif --restart unless-stopped -v ~/Documents/bioinfo_tsinghua_share:/data gangxu/motif:1.0
# docker exec -it -u root motif bash
cd /home/test/motif/sequence_motif/practice

## Fetch UTR or promoter sequences (in R)
# view files with cat, head, less, column, notepad, vscode, and so on

## Intersect with interested genes (in linux)

# interested 3'UTR
sort -t $'\t' -k 2 utr3p.info | join -o 1.3 2.1 1.2 1.9 1.4 1.5 1.6 1.7 1.8 1.10 -t $'\t' -1 2 -2 2 - <(cut -f 1 ../SC2_SF2.ct.dn.1_0.01.protein_coding | sort | join -t $'\t' -1 1 -2 1 - <(sort -t $'\t' -k 1 <(grep -o -P -e "gene_id.*; transcript_id.*?;" ../genome/gencode.v27.annotation.gtf | sort | uniq | sed -e 's/gene_id "//' -e 's/"; transcript_id "/\t/' -e 's/";//'))| sort -t $'\t' -k 2  ) | sort -t $'\t' -k 1 > interested_three_prime_UTR.info
# column1: chr; column2: gene; column3: transprict; column4: exon_name; column5: start; column6: end; column7: width(seq length); column8: strand; column9: exon_id; column10: exon_rank

# interested promoter
sort -t $'\t' -k 7 promoter.info | join -o 1.1 2.1 1.7 1.2  1.3 1.4 1.5 1.6 -t $'\t' -1 7 -2 2 - <(cut -f 1 ../SC2_SF2.ct.dn.1_0.01.protein_coding | sort | join -t $'\t' -1 1 -2 1 - <(sort -t $'\t' -k 1 <(grep -o -P -e "gene_id.*; transcript_id.*?;" ../genome/gencode.v27.annotation.gtf | sort | uniq | sed -e 's/gene_id "//' -e 's/"; transcript_id "/\t/' -e 's/";//' ))| sort -t $'\t' -k 2  ) | sort -t $'\t' -k 1 > interested_promoter.info
# column1: chr; column2: gene; column3: transprict; column4: start; column5: end; column6: width(seq length); column7: strand; column8: transprict_id

## Convert to bed format (in linux)

# UTR bed info
cat interested_three_prime_UTR.info | \
  awk '{print $1 "\t" $5-1 "\t" $6 "\t" $3 "\t" $2 "\t" $8}' | \
  sort -u  > interested_three_prime_UTR.bed
  # column1: chr; column2: start; column3: end; column4: gene; column5: transcript; column6: strand

# promoter bed info
cat interested_promoter.info | \
  awk '{print $1 "\t" $4-1 "\t" $5 "\t" $3 "\t" $2 "\t" $7}' | \
  sort -u  > interested_promoter.bed
  # column1: chr; column2: start; column3: end; column4: gene; column5: transcript; column6: strand

## Get genome sequence (firstly in linux then in R)

# get 3'UTR related genome sequence
bedtools getfasta -s -name -fi ../genome/GRCh38.p10.genome.fa \
  -bed interested_three_prime_UTR.bed -fo interested_three_prime_UTR.fa

# get promoter related genome sequence
bedtools getfasta -s -name -fi ../genome/GRCh38.p10.genome.fa \
  -bed interested_promoter.bed -fo interested_promoter.fa

# concatenate sequences of the same 3’ UTR; concatenate sequences of the same promoter (in R)

## Generate random sequence as background sequence (in linux except specified)

# shuffle the input sequence
fasta-shuffle-letters interested_three_prime_UTR.fa interested_three_prime_UTR.control
fasta-shuffle-letters interested_promoter.fa interested_promoter.control

# downstream 1000bp (in R)

# bedtools shuffle
bedtools shuffle -i interested_three_prime_UTR.bed \
-g ../genome/hg38.chrom.sizes >interested_three_prime_UTR_btools.bed

bedtools shuffle -i interested_promoter.bed \
-g ../genome/hg38.chrom.sizes >interested_promoter_btools.bed

## Motif enrichment (in linux)

# de novo motif discovery
cd /home/test/motif/sequence_motif/practice/
meme -dna -maxsize 2000000 \    # why 2000000? try 1000000
  -minw 4 -maxw 12 \
  -oc promoter_de_novo \
  -nmotifs 5 \
  interested_promoter.fa

# known motif enrichment
mkdir UTR_output
ame \
--control interested_three_prime_UTR.control \
--oc UTR_output/ \
interested_three_prime_UTR.fa \
../Homo_sapiens.meme \
../Ray2013_rbp_Homo_sapiens.meme

mkdir promoter_output
ame \
--control interested_promoter.control \
--oc promoter_output/ \
interested_promoter.fa \
../JASPAR2018_CORE_vertebrates_non-redundant.meme \
../HOCOMOCOv11_core_HUMAN_mono_meme_format.meme
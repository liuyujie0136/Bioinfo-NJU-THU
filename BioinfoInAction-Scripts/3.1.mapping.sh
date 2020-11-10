#!/bin/bash

cd /home/test/mapping

# mapping
bowtie -v 2 -m 10 --best --strata BowtieIndex/YeastGenome -f THA1.fa -S THA1.sam
bowtie -v 1 -m 10 --best --strata bowtie-src/indexes/e_coli -q e_coli_1000_1.fq -S e_coli_1000_1.sam

# convert file format
perl sam2bed.pl THA1.sam > THA1.bed

# clear files
grep -v chrmt THA1.bed > THA1_new.bed	# without chr MT
grep $'chrIV\t' THA1.bed > THA1_chrIV.bed	# only chr IV

cp THA1_new.bed /home/test/share
cp THA1_chrIV.bed /home/test/share


## mapping paired-end reads (Note: not in this docker!)
# docker load -i ~/Downloads/bioinfo_pairend.tar.gz
# docker run --name=bioinfo_pairend -dt -h bioinfo_docker --restart unless-stopped -v ~/Documents/bioinfo_tsinghua_share:/home/test/share gangxu/bioinfo_pairend:1.0
# docker exec -it  bioinfo_pairend bash

cd /home/test/pair_end_mapping

STAR --genomeDir ./STAR_TAIR10_index --readFilesIn Ath_1.fq Ath_2.fq --outFileNamePrefix ./output/ --outSAMtype BAM SortedByCoordinate


#!/bin/bash

## tutorial
cd /home/test/alter-spl/

# check read length
samtools view input/SRR065544_chrX.bam | cut -f 10 | perl -ne 'chomp;print length($_) . "\n"' | sort | uniq -c
samtools view input/SRR065545_chrX.bam | cut -f 10 | perl -ne 'chomp;print length($_) . "\n"' | sort | uniq -c

# run
echo "input/SRR065544_chrX.bam" > input/b1.txt
echo "input/SRR065545_chrX.bam" > input/b2.txt

python2 /usr/local/rMATS-turbo-Linux-UCS4/rmats.py --b1 input/b1.txt --b2 input/b2.txt --gtf input/Mus_musculus_chrX.gtf --od output -t paired --readLength 35


## homework
cd /home/test/share/alter-spl-homework/
samtools view SRR065546_chrX.bam | cut -f 10 | perl -ne 'chomp;print length($_) . "\n"' | sort | uniq -c
samtools view SRR065547_chrX.bam | cut -f 10 | perl -ne 'chomp;print length($_) . "\n"' | sort | uniq -c
echo "SRR065546_chrX.bam" > b1.txt
echo "SRR065547_chrX.bam" > b2.txt
python2 /usr/local/rMATS-turbo-Linux-UCS4/rmats.py --b1 b1.txt --b2 b2.txt --gtf /home/test/alter-spl/input/Mus_musculus_chrX.gtf --od output -t paired --readLength 35

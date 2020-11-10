#!/bin/bash

# sudo apt -y install ncbi-blast+
cd /home/test/blast/
mkdir output

## Pairwise sequence alignment
# Protein sequence alignment
blastp -query protein/VIM.fasta -subject protein/NMD.fasta -out output/blastp

# DNA sequence alignment
blastn -query dna/H1N1-HA.fasta -subject dna/H7N9-HA.fasta -out output/blastn

## Align a sequence to a remote database
blastp -query protein/VIM.fasta -db pdb -remote -out output/blastp_remote
blastn -query dna/H1N1-HA.fasta -db nr -remote -out output/blastn_remote

## Align a sequence to a local database
makeblastdb -dbtype nucl -in dna/YeastGenome.fa -out database/YeastGenome
blastn -query dna/Yeast.fasta -db database/YeastGenome -out output/Yeast.blastn



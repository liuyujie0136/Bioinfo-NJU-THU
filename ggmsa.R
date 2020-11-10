### ggplot2 Based Multiple Sequence Alignment Visualization

install.packages('ggmsa')

library(ggplot2)
library(ggmsa)

## Quick Example
# Protein sequences
protein_sequences <- system.file("extdata", "sample.fasta", package = "ggmsa")
ggmsa(protein_sequences, 164, 213, color = "Chemistry_AA")

# DNA sequences
nt_sequence <- system.file("extdata", "LeaderRepeat_All.fa", package = "ggmsa")
ggmsa(nt_sequence, color = "Chemistry_NT")

# RNA Sequences
miRNA_sequences <- system.file("extdata", "seedSample.fa", package = "ggmsa")
ggmsa(miRNA_sequences, color = "Chemistry_NT")

# null font 
nt_sequence <- system.file("extdata", "LeaderRepeat_All.fa", package = "ggmsa")
ggmsa(nt_sequence, font = NULL, color = "Chemistry_NT")

# Visualizing multiple sequence alignment with sequence logo
f <- system.file("extdata", "LeaderRepeat_All.fa", package = "ggmsa")
ggmsa(f, font = NULL, color = "Chemistry_NT" ) + geom_seqlogo()

# Visualizing multiple sequence alignment with GC content
f <- system.file("extdata", "LeaderRepeat_All.fa", package = "ggmsa")
ggmsa(f, font = NULL, color = "Chemistry_NT" ) + geom_GC()

# Visualizing multiple sequence alignment with ggtree
library(Biostrings)
x <- readAAStringSet(protein_sequences)
d <- as.dist(stringDist(x, method = "hamming")/width(x)[1])
library(ape)
tree <- bionj(d)
library(ggtree)
p <- ggtree(tree) + geom_tiplab()

data = tidy_msa(x, 164, 213)
p + geom_facet(geom = geom_msa, data = data,  panel = 'msa',
               font = NULL, color = "Chemistry_AA") +
  xlim_tree(1)

# Highlighting the seed in miRNA sequences
ggmsa(miRNA_sequences, font = 'DroidSansMono', none_bg = TRUE) +
  geom_seed(seed = "GAGGUAG") + theme_void()  # shaded
ggmsa(miRNA_sequences, font = 'DroidSansMono', color = "Chemistry_NT") +
  geom_seed(seed = "GAGGUAG", star = TRUE) + theme_void()  # highlighted using star

# Breaking Down MSA
ggmsa(protein_sequences, end = 300, font = NULL, color = "Chemistry_AA") + facet_msa(field = 100)  # 3 fields
ggmsa(protein_sequences, end = 400, color = "Chemistry_AA") + facet_msa(field = 100)  # 4 fields


## real data - basic

# read seq
seq=Biostrings::readDNAStringSet("phylo-ggtree/p53.fa",format="fasta")

# do multi-seq alignment
aln=muscle::muscle(seq)

# plot alignment result
library(ggplot2)
library(ggmsa)
ggmsa(aln,start=40,end=60,font='TimesNewRoman',color='Chemistry_NT')

# export plot (the file may be too large to open, be cautious!)
export::graph2ppt(file=paste0("Rplot-",gsub(":","-",gsub(" ","-",as.character(Sys.time()))),".pptx"))


## real data - with ggtree



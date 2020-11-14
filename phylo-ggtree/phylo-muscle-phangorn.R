### Phylogenetic analysis in R

library(Biostrings)
library(muscle)
library(ape)
library(phangorn)

## Read input seq using Biostrings package
# Note: p for protein, d for DNA

p.seq=readAAStringSet("NANOG.fa", format="fasta", nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
d.seq=readDNAStringSet("p53.fa", format="fasta", nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)


## Produce multiple sequence alignments using MUSCLE

p.aln=muscle::muscle(p.seq)
d.aln=muscle::muscle(d.seq)


## Convert the data format

p.phydat=as.phyDat(p.aln,type='AA')
d.phydat=as.phyDat(d.aln,type='DNA')


## UPGMA and NJ method

# compute pairwise distances
p.dm=dist.ml(p.phydat,model='JTT')
d.dm=dist.ml(d.phydat,model='JC69')

# growing a tree
p.treeUPGMA=upgma(p.dm)
p.treeNJ=NJ(p.dm)
d.treeNJ=NJ(d.dm)

# basic plot
plot(p.treeUPGMA, type="fan", main="UPGMA")
plot(p.treeNJ, type="unrooted", main="NJ")
plot(d.treeNJ)


## Maximum Parsimony analysis
# use the UPGMA and NJ tree as starting trees
  # for the maximum parsimony and maximum likelihood analyses!

# calculate parsimony score
parsimony(treeUPGMA,protein.phydat)
parsimony(treeNJ,protein.phydat)

# find lower parsimony score
treePars=optim.parsimony(treeUPGMA,protein.phydat)

# a version of the parsimony ratchet implemented, 
  # which is likely to find better trees than just doing NNI/SPR rearrangements.
treeRatchet=pratchet(protein.phydat, trace = 0)


## Maximum likelihood method

fit=pml(treeNJ, data=protein.phydat)  # compute the likelihood for a tree

fitJC=optim.pml(fit, TRUE)  # optimize branch length for Jukes-Cantor model
logLik(fitJC)

fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="?", optInv=TRUE, optGamma=TRUE, rearrangement = "NNI", control = pml.control(trace = 0))
fitGTR

# For optim.pml()
# nucleotide models: JC, F81, K80, HKY, TrNe, TrN, TPM1, K81, TPM1u, TPM2, TPM2u, TPM3, TPM3u, TIM1e, TIM1, TIM2e, TIM2, TIM3e, TIM3, TVMe, TVM, SYM and GTR. 
# amino acid models: "WAG", "JTT", "LG", "Dayhoff", "cpREV", "mtmam", "mtArt", "MtZoa", "mtREV24", "VT","RtREV", "HIVw", "HIVb", "FLU", "Blossum62", "Dayhoff_DCMut" and "JTT_DCMut"

# model test!
protein.mt=modelTest(protein.phydat,treeNJ,model='all')
dna.mt=modelTest(dna.phydat,treeNJ,c("JC","GTR"))  #?

# returns class 'pml', not standard 'phylo'
methods(class = "pml")
fitGTR.tree=fitGTR$tree  # phylo class


## Bootstrap
protein.bs=bootstrap.pml(fitGTR, bs = 100, optNni = TRUE)
plotBS(fitGTR.tree, protein.bs)





## Export plots to PPT
library(export)
t=gsub(":","-",gsub(" ","-",as.character(Sys.time())))
graph2ppt(file=paste0("Rplot-",t,".pptx"))

## save RData (important for ggtree.R!)
save.image(file="phylo-ggtree.RData") 

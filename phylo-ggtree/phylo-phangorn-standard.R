### Standard scripts for phylogenetic analysis in R using package 'phangorn'


## Standard scripts for nucleotide analysis

library(phangorn)

dat = read.phyDat("myfile")
dm = dist.ml(dat, "F81")
tree = NJ(dm)

# as alternative for a starting tree
tree <- pratchet(dat) # parsimony tree
tree <- nnls.phylo(tree, dm) # need edge weights

# 1. alternative: quick and dirty: GTR + G
fitStart = pml(tree, dat, k=4)
fit = optim.pml(fitStart, model="GTR", optGamma=TRUE, rearrangement="stochastic")

# 2. alternative: preper with modelTest
mt <- modelTest(dat, tree=tree, multicore=TRUE)
mt[order(mt$AICc),]

# choose best model from the table according to AICc
bestmodel <- mt$Model[which.min(mt$AICc)]
env = attr(mt, "env")
fitStart = eval(get("GTR+G+I", env), env)
# or let R search the table
fitStart = eval(get(bestmodel, env), env)

fit = optim.pml(fitStart, rearrangement = "stochastic", optGamma=TRUE, optInv=TRUE, model="GTR")

bs = bootstrap.pml(fit, bs=100, optNni=TRUE, multicore=TRUE)


## Standard scripts for amino acid analysis

library(phangorn)

dat = read.phyDat("myfile", type = "AA")
dm = dist.ml(dat, model="JTT")
tree = NJ(dm)

# parallel will only work safely from command line
mt <- modelTest(dat, model=c("JTT", "LG", "WAG"),multicore=TRUE)

# run all available amino acid models
mt <- modelTest(dat, model="all", multicore=TRUE)

fitStart = eval(get(mt$Model[which.min(mt$BIC)], env), env)
fitNJ = pml(tree, dat, model="JTT", k=4, inv=.2)

fit = optim.pml(fitNJ, rearrangement = "stochastic", optInv=TRUE, optGamma=TRUE)

bs = bootstrap.pml(fit, bs=100, optNni=TRUE, multicore=TRUE)



## Standard script for ML analysis

library(multicore)
library(phangorn)
dat = read.phyDat("myfile")
dm = dist.ml(dat)
tree = fastme.bal(dm)
mT = modelTest(dat, tree)
fit = pml(tree, dat, k = 4, inv = 0.2)
fit2 = optim.pml(fit, optNni = TRUE, optGamma = TRUE, optInv = TRUE, model = "GTR")
bs = bootstrap.pml(fit2, bs = 100, optNni = TRUE)
plotBS(fit2$tree, bs)

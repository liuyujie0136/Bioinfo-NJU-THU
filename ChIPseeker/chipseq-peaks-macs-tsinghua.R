# Prepare the environment of software
rm(list = ls())
options(stringsAsFactors = FALSE)

# Library packages
library("RMariaDB")
library("rtracklayer")
library("ChIPseeker")
library("org.Sc.sgd.db")
library("clusterProfiler")
library("GenomicFeatures")
library("ggplot2")
library("ggupset")

# Display available genomes list
ucscGenomes()[ , "db"]
# Display the list of table known to work
supportedUCSCtables()
# Retrieve full datasets (responds to genome ID of fasta) for Yeast from UCSC
TxDb.Ssaccer.USUC.sacCer2.ensGene <- makeTxDbFromUCSC(genome = "sacCer2", tablename = "ensGene")

# Find the position of peak in WGS
txdb <- TxDb.Ssaccer.USUC.sacCer2.ensGene
# Read the exported peak file from MACS
yeast <- readPeakFile("yeast_macs_p05_peaks.narrowPeak")
covplot(yeast, weightCol = 5)

# ChIP peaks binding TSS region
promoter <- getPromoters(TxDb = txdb, upstream = 3000, downstream = 3000)
tagMatrix <- getTagMatrix(peak = yeast, windows = promoter)
tagHeatmap(tagMatrix = tagMatrix, xlim = c(-3000, 3000), color = "red")

# Average profile of ChIP peaks binding to TSS region
plotAvgProf(tagMatrix, xlim = c(-3000, 3000), conf = 0.95, resample = 1000,
            xlab = "Genomone Region (5'->3')", ylab = "Read Count Frequency")

# Peaks annotation
peakAnno <- annotatePeak(yeast, tssRegion = c(-3000, 3000), TxDb = txdb, annoDb = "org.Sc.sgd.db")

# Visualize the genome annotation
plotAnnoPie(peakAnno)
plotAnnoBar(peakAnno)
vennpie(peakAnno)
upsetplot(peakAnno)

# Visualize distribution of TF-binding loci relative to TSS
plotDistToTSS(peakAnno)

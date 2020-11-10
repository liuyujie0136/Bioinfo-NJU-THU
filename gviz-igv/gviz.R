### Using Gviz package to mimic Genome Browser - standard example
## Install
BiocManager::install("Gviz")

## Library
library(Gviz)
library(GenomicRanges)

## acquire data
data(cpgIslands)
data(geneModels)

chr <- as.character(unique(seqnames(cpgIslands))) # get names of chromosome
gen <- genome(cpgIslands) # get genome reference name used in UCSC

## get track info (track means a small seq mapping to the gemome)
atrack <- AnnotationTrack(cpgIslands, name = "CpG")
plotTracks(atrack)

## add axis lables
gtrack <- GenomeAxisTrack()
plotTracks(list(gtrack, atrack))

## add genome picture and location
itrack <- IdeogramTrack(genome = gen,chromosome = chr)
plotTracks(list(itrack, gtrack, atrack))

## add gene model
grtrack <- GeneRegionTrack(geneModels,genome = gen,chromosome = chr, name = "Gene Model")
plotTracks(list(itrack, gtrack, atrack, grtrack))

## annotate each model
grtrack2 <- GeneRegionTrack(geneModels,genome = gen,chromosome = chr, name = "Gene Model",transcriptAnnotation = "symbol")
plotTracks(list(itrack, gtrack, atrack, grtrack2))

## show models in selected region
plotTracks(list(itrack, gtrack, atrack, grtrack),from = 26700000, to = 26750000)

## show actual genomic sequence at a given position
library(BSgenome.Hsapiens.UCSC.hg19)  # too big to download!
strack <- SequenceTrack(Hsapiens, chromosome = chr)
plotTracks(list(itrack, gtrack, atrack, grtrack, strack), from = 26591822, to = 26591852, cex = 0.8)

## add other data using (numeric) DataTrack, e.g. deepseq data
# here, generate random data for example
set.seed(255)
lim <- c(26700000, 26750000)
coords <- sort(c(lim[1], sample(seq(from= lim[1],to = lim[2]), 99), lim[2]))
dat <- runif(100, min = -10, max = 10)
dtrack <- DataTrack(data = dat, start =coords[-length(coords)],end = coords[-1], chromosome = chr, genome = gen,name ="Uniform")
plotTracks(list(itrack, gtrack, atrack, grtrack, dtrack), from = lim[1], to = lim[2])  # default: dot plot
plotTracks(list(itrack, gtrack, atrack, grtrack, dtrack), from = lim[1], to = lim[2], type = "histogram")  # change plot type

## reverse the strand (3'->5')
plotTracks(list(itrack, gtrack, atrack, grtrack), reverseStrand = TRUE)

## data grouping
data(twoGroups)
dTrack <- DataTrack(twoGroups, name ="uniform")
plotTracks(dTrack, groups = rep(c("control", "treated"), each = 3), type = c("a","p", "confint"))

## multiple samples in different lines
data(dtHoriz)
dtHoriz <- dtHoriz[1:6, ]
plotTracks(dtHoriz, type ="horiz", groups = rownames(values(dtHoriz)),showSampleNames = TRUE,cex.sampleNames = 0.6,separator = 1)

## plot read alignments data from NGS
afrom <- 2960000
ato <- 3160000
alTrack <-AlignmentsTrack(system.file(package = "Gviz","extdata","gapped.bam"), isPaired = TRUE)
bmt <- BiomartGeneRegionTrack(genome = "hg19", chromosome ="chr12",start = afrom, end = ato, filter = list(with_ox_refseq_mrna =TRUE),stacking = "dense")
plotTracks(c(bmt, alTrack), from = afrom,to = ato, chromosome ="chr12")
plotTracks(c(alTrack, bmt), from = afrom,to = ato, chromosome = "chr12", type = "coverage")  # only show peaks above

## highlighting tracks
ht <- HighlightTrack(trackList = list(atrack, grtrack, dtrack), 
                     start = c(26705000, 26720000), width = 7000,
                     chromosome = 7)
plotTracks(list(itrack, gtrack, ht), from = lim[1], to = lim[2])

## composite plots for multiple chromosomes
# generate foo data
chroms <- c("chr1", "chr2", "chr3", "chr4")
maTrack <- AnnotationTrack(range=GRanges(seqnames = chroms, 
                                         ranges = IRanges(start = 1,  width = c(100, 400, 200,1000)),
                                         strand = c("+", "+", "-", "+")), genome = "mm9", 
                           chromosome = "chr1", name = "foo")

mdTrack <- DataTrack(
  range = GRanges(seqnames = rep(chroms, c(10, 40, 20, 100)),
                  ranges = IRanges(start = c(seq(1, 100, len = 10),
                                             seq(1, 400, len = 40), 
                                             seq(1, 200, len = 20),
                                             seq(1, 1000, len = 100)), 
                                   width = 9), values = runif(170)),
  data = "values", chromosome = "chr1", genome = "mm9", name = "bar")

mgTrack <- GenomeAxisTrack(scale = 50, labelPos = "below", exponent = 3)
itrack <- IdeogramTrack(genome = "mm9",chromosome = "chr1")

# use grid for composite plot
ncols <- 2
nrows <- length(chroms) %/% ncols
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrows, ncols)))
for(i in seq_along(chroms)) {
  pushViewport(viewport(layout.pos.col = ((i - 1) %% ncols) + 1,
                        layout.pos.row = (((i) - 1) %/% ncols) + 1))
  plotTracks(list(itrack, maTrack, mdTrack, mgTrack), 
             chromosome = chroms[i], add = TRUE)
  popViewport(1)
}

# or use lattice instead of grid
library(lattice)
chroms <- data.frame(chromosome = chroms)
xyplot(1 ~ chromosome | chromosome, data = chroms, panel = function(x) {
  plotTracks(list(itrack , maTrack, mdTrack, mgTrack), 
             chromosome = x, add = TRUE, showId = FALSE) },
  scales = list(draw = FALSE), xlab = NULL, ylab = NULL)

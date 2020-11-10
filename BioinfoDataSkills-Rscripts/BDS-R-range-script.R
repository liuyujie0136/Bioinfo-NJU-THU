## Selected R scripts about range data in book Bioinformatics Data Skills (chapter9)
# Note: comments are below or on the right hand of the code with one #, titles with ##


## 0-based and 1-based
# 0-based coordinate system, with half-closed, half-open intervals: [start, end), width=end-start
# 1-based coordinate system, with closed intervals: [start, end], width=end-start+1

# BED, BAM, BCF: 0-based
# GTF, GFF, SAM, VCF, Wiggle, GenomicRanges, BLAST, GenBank/EMBL Feature Table: 1-based


## Install packages from Bioconductor (for R 3.5 or later!)
install.packages("BiocManager")
BiocManager::install("GenomicRanges")


## IRanges
library(IRanges)
x <- IRanges(start=c(40, 80), end=c(67, 114), names=letters[1:2])
x + 4L # both end extended. results: start end width [1] 36 71 36 a [2] 76 118 43 b
restrict(x, 50, 100)  # selected range
flank(x, width=7)   # Flanking ranges (on the side of). downstream when 'start=FALSE'
reduce()  # is useful when all we care about is what regions of a sequence are covered
gaps()  # returns the gaps (uncovered portions) between ranges
gaps(alns, start=1, end=60)  # all the gaps including first and last
# think about range operations as set operations like difference (setdiff()), intersection (intersect()), union (union()), and complement (which is simply the function gaps())


## Finding Overlapping Ranges
qry <- IRanges(start=c(1, 26, 19, 11, 21, 7), end=c(16, 30, 19, 15, 24, 8), names=letters[1:6])
sbj <- IRanges(start=c(1, 19, 10), end=c(5, 29, 16), names=letters[24:26])
hts <- findOverlaps(qry, sbj)
# results: each row contains the index of a query range that overlaps a subject range, and the index of the subject range it overlaps. can be accessed by queryHits() and subjectHits()
names(qry)[queryHits(hts)]
names(sbj)[subjectHits(hts)]

findOverlaps(qry, sbj, type="within")
# limit our overlap results to only include query ranges that fall entirely within subject ranges
findOverlaps(qry, sbj, select="first")  # find first overlapping subject
countQueryHits(hts)  # count query hits in subject
setNames(countQueryHits(hts), names(qry))  # name for clarity
ranges(hts, qry, sbj)  # create a set of ranges for overlapping regions
countOverlaps(qry, sbj)  # counting overlaps
subsetByOverlaps(qry, sbj)  # subset of queries that overlap subjects


## Finding Nearest Ranges and Calculating Distance 
nearest(qry, sbj)  # returns the nearest range, even if overlaps
precede(qry, sbj)  # nearest range that the query is upstream of subject
follow(qry, sbj)  # downstream
distanceToNearest(qry, sbj)  # find neighboring ranges
distance(qry, sbj)   # each pairwise distances between query and subject ranges


## Run Length Encoding and Views 
set.seed(0)
rngs <- IRanges(start=sample(seq_len(60), 10), width=7)
rngs_cov <- coverage(rngs)
# takes a set of ranges and returns their coverage
rngs_cov[rngs[9]]  # the coverage in the region 9
mean(rngs_cov[rngs[9]])  # mean coverage within this range

min_cov2 <- slice(rngs_cov, lower=2)  # regions with more than 2x coverage (return a 'view' object)
# using a subset of sequence to define new ranges (useful!)
ranges(min_cov2)  # without the right part
viewMeans(min_cov2)
# summarize the views using viewMeans(), viewMaxs(), viewApply()

length(rngs_cov)
bwidth <- 5L
end <- bwidth * floor(length(rngs_cov) / bwidth)
windows <- IRanges(start=seq(1, end, bwidth), width=bwidth) 
cov_by_wnd <- Views(rngs_cov, windows)
# view by slide windows


## GenomicRanges
library(GenomicRanges)
gr <- GRanges(seqname=c("chr1", "chr1", "chr2", "chr3"), ranges=IRanges(start=5:8, width=10), strand=c("+", "-", "-", "+")) 
seqlengths(gr) <- c(chr1=152, chr2=432, chr3=903)  # can also add 'seqlengths=...' in GRanges
# other operations on 'gr'(GRanges object):
#  seqnames(), strand(), start(), end(), width(), ranges(), length(), names(), and comparisons

grm <- GRanges(seqname=c("chr1", "chr1", "chr2", "chr3"), ranges=IRanges(start=5:8, width=10), strand=c("+", "-", "-", "+"), gc=round(runif(4), 3), num=1:4)
# can add arbitrary metadata columns 'gc','num'
mcols(grm)  # access metadata columns
grm$gc  # access specific metadata columns
mcols(grm[seqnames(grm) == "chr1"])$gc


## Grouping Data with GRangesList
gr1 <- GRanges(c("chr1", "chr2"), IRanges(start=c(32, 95), width=c(24, 123)))
gr2 <- GRanges(c("chr8", "chr2"), IRanges(start=c(27, 12), width=c(42, 34)))
grl <- GRangesList(gr1, gr2) 
doubled_grl <- c(grl, grl)  # combine
# 'grl' can also have functions for GRanges data

gr_split <- split(gr, seqnames(gr))  # split by seqnames (for large data)
unsplit(gr_split, seqnames(gr))  # unsplit
lapply(gr_split, function(x) order(width(x)))  # then can lapply and sapply
# reduce(), flank(), coverage(), and findOverlaps() can work directly with GRangesList objects


## Working with Annotation Data: GenomicFeatures and rtracklayer
library(GenomicFeatures)
BiocManager::install("TxDb.Mmusculus.UCSC.mm10.ensGene")
library(TxDb.Mmusculus.UCSC.mm10.ensGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.ensGene  # TranscriptDb for GenomicFeatures

mm_genes <- genes(txdb)  # a GRanges object
# other functions in GenomicFeatures:
#  transcripts(), exons(), cds(), and promoters()

mm_exons_by_tx <- exonsBy(txdb, by="tx")
mm_exons_by_gn <- exonsBy(txdb, by="gene") 
# all exons grouped by transcript or gene,  by=gene, tx, exon, or cds
# more functions:
#  transcriptsBy(), exonsBy(), cdsBy(), intronsBy(), fiveUTRsByTranscript(), and threeUTRsByTranscript()

qtl_region <- GRanges("chr8", IRanges(123260562, 123557264))
qtl_region_expanded <- qtl_region + 10e3  # 10e3=1e4=10kbp(e=10^)
transcriptsByOverlaps(txdb, qtl_region_expanded)
# extract feature data that only overlaps a specific region
# more functions: exonsByOverlaps(), and cdsByOverlaps()

library(rtracklayer)
mm_gtf <- import('Mus_musculus.GRCm38.75_chr1.gtf.gz')
colnames(mcols(mm_gtf))  # metadata columns, use '$' to access
# use rtracklayer package to import gtf files (can also export)
vignette("rtracklayer")  # ref in pdf

chr1_pcg <- mm_gtf[mm_gtf$type == "gene" & mm_gtf$gene_biotype == "protein_coding"]  # only one '&', '&&' is wrong!!
chr1_pcg_3kb_up <- flank(chr1_pcg, width=3000) 
chr1_pcg_3kb_up2 <- promoters(chr1_pcg, upstream=3000, downstream=0)
# extract promoters


## Connect GenomicRanges with Sequence Data
mcols(chr1_pcg_3kb_up) <- NULL
export(chr1_pcg_3kb_up, con="chr1_pcg_promoters.bed", format="BED")
# export as bed and using BEDTools, in package rtracklayer

# BiocManager::install("BSgenome"); BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")
# library(BSgenome.Mmusculus.UCSC.mm10); mm_gm <- BSgenome.Mmusculus.UCSC.mm10 
# Check: seqlevels(); seqlevelsStyle(); all(seqlevels(chr1_pcg_3kb_up) %in% seqlevels(mm_gm))
chr1_3kb_seqs <- getSeq(mm_gm, chr1_pcg_3kb_up)
writeXStringSet(chr1_3kb_seqs, file="mm10_chr1_3kb_promoters.fasta", format="fasta")
# get sequence from R packages


## GRanges in practice
gaps()  # returns all the gaps including start, two strands, etc.

## Representing the introns of transcripts
mm_introns <- intronsByTranscript(txdb)
mm_introns[names(amy1_tx)]
# using intronsByTranscript()

amy1 <- transcriptsBy(txdb, 'gene')$ENSMUSG00000074264
mm_exons <- exonsBy(txdb, "tx")
amy1_tx <- split(amy1, amy1$tx_id)
amy1_exons <- mm_exons[match(names(amy1_tx), names(mm_exons))]
amy1_introns <- setdiff(amy1_tx, amy1_exons)
# manual approach

identical(mm_introns[names(amy1_tx)], amy1_introns) 
# check two methods are identical

## Finding and Working with Overlapping Ranges
findOverlaps(dbsnp137_resized, chr1_collapsed_exons, ignore.strand=TRUE)
subsetByOverlaps(dbsnp137_resized, chr1_collapsed_exons, ignore.strand=TRUE)
countOverlaps(chr1_collapsed_exons, dbsnp137_resized, ignore.strand=TRUE)  # reversed arguments
coverage(reads)  # coverage



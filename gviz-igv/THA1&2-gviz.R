### Use Gviz to visualize THA1&2 data

## THA1_chrIV with gene model
library(Gviz)
library(GenomicRanges)
library(GenomicFeatures)

tha1=read.table("THA1_chrIV.bed")
tha1.df=as.data.frame(tha1)
tha1.gr=GRanges(seqname=tha1.df$V1,ranges=IRanges(start=tha1.df$V2,end=tha1.df$V3),strand=tha1.df$V6)

chr=as.character(unique(seqnames(tha1.gr)))
genome(tha1.gr)='sacCer3'
gen=genome(tha1.gr)

# make gene model
supportedUCSCtables(genome="sacCer3")  # yeast genome
txdb=makeTxDbFromUCSC(genome='sacCer3',tablename='ensGene')

# OTHER methods to generate a TxDb object - using GenomicFeatures
# loadDb("file.sqlite")
# makeTxDbFromBiomart
# makeTxDbFromEnsembl
# makeTxDbFromGFF

# visualization
atrack=AnnotationTrack(tha1.gr,name="THA1_chrIV")
gtrack=GenomeAxisTrack()
itrack=IdeogramTrack(genome=gen,chromosome=chr)
grtrack=GeneRegionTrack(txdb,genome=gen,chromosome=chr,name="Gene Model",transcriptAnnotation="gene")
plotTracks(list(itrack,gtrack, atrack, grtrack))
plotTracks(list(itrack,gtrack, atrack, grtrack), from=700000, to=800000)


## THA1 - Multiple chromosomes
library(Gviz)
library(GenomicRanges)

tha1=read.table("THA1_new.bed")
tha1.df=as.data.frame(tha1)
tha1.gr=GRanges(seqname=tha1.df$V1,ranges=IRanges(start=tha1.df$V2,end=tha1.df$V3),strand=tha1.df$V6) 
tha1.gr.split <- split(tha1.gr, seqnames(tha1.gr))

gen='sacCer3'

for (i in 1:16) {
  
  chr=as.character(unique(seqnames(tha1.gr.split[i])))
  
  atrack=AnnotationTrack(tha1.gr.split[i],name=fname)
  gtrack=GenomeAxisTrack()
  itrack=IdeogramTrack(genome=gen,chromosome=chr)
  plotTracks(list(itrack, gtrack, atrack))
  
  export::graph2pdf(file=paste0(paste('THA1',chr,sep='_'),"-",gsub(":","-",gsub(" ","-",as.character(Sys.time()))),".pdf"))
  
}


## THA1 - Multiple chromosomes - in one picture




## THA2_XII - homework
library(Gviz)
library(GenomicRanges)

tha2=read.table("THA2_XII.bed")
tha2.df=as.data.frame(tha2)
tha2.gr=GRanges(seqname=tha2.df$V1,ranges=IRanges(start=tha2.df$V2,end=tha2.df$V3),strand=tha2.df$V6)

chr=as.character(unique(seqnames(tha2.gr)))
genome(tha2.gr)='sacCer3'
gen=genome(tha2.gr)

atrack=AnnotationTrack(tha2.gr,name="THA2_chrXII")
gtrack=GenomeAxisTrack()
itrack=IdeogramTrack(genome=gen,chromosome=chr)
plotTracks(list(itrack,gtrack,atrack))
           
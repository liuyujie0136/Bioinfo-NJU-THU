## Install package "fgsea" from github (Rtools required)
library(devtools)
install_github("ctlab/fgsea")

## library packages
library(data.table)
library(fgsea)
library(ggplot2)

## input files introduction
# 1)test.rnk: RNK file contains a single, rank-ordered gene list (not gene set) in a simple newline-delimited text format.
# 2)test.gmx(for gsea desktop) / test.gmt(for fgsea):
  # The GMX/GMT file format is a tab-delimited file format that describes gene sets.
  # This file can contain multiple gene sets: In the GMX format, each column represents a gene set; In the GMT format, each line represents a gene set.

## Run fGSEA
# first input .rnk and .gmt 
rnk.file <- "test.rnk"
gmt.file <- "test.gmt"

# extract ranklist
ranks <- read.table(rnk.file, header=T, colClasses = c("character", "numeric"))
ranks <- setNames(ranks$Rank, ranks$GeneID)

# extract gene sets
pathways <- gmtPathways(gmt.file)

# check your ranklist and geneset
str(ranks)
str(head(pathways))

# run user-defined GSEA. Questions about warning message, refer to https://www.gitmemory.com/issue/ctlab/fgsea/79/703686364
fgseaRes <- fgsea(pathways, ranks, nperm = 10, minSize=0, maxSize=500)
head(fgseaRes[order(pval), ])

# plot for pathways
topPathwaysUp = fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown = fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways = c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(pathways[topPathways], ranks, fgseaRes, gseaParam=0.5)



## Download data with qGSEA (with example GSE19161)
library(devtools)
install_github('dongzhuoer/qGSEA')

library(qGSEA)
dir.create('raw/')
download.file(rGEO::gse_soft_ftp('GSE19161'), 'raw/GSE19161_family.soft.gz')
download.file(rGEO::gse_matrix_ftp('GSE19161'), 'raw/GSE19161_series_matrix.txt.gz')

dir.create('input/')
qGSEA::make_gsea_input(
  matrix_file = 'raw/GSE19161_series_matrix.txt.gz',
  soft_file = 'raw/GSE19161_family.soft.gz',
  output_dir = 'input/',
  gene = 'EIF4G2'
  )



## Homework
# get rnk and gmt files
rnk.file="homework.rnk"
hw.gmx=read.table("homework.gmx")
hw.gmt=as.data.frame(t(hw.gmx))
write.table(hw.gmt,file="homework.gmt",quote=F,sep="\t",row.names=F,col.names=F)
gmt.file="homework.gmt"

# extract ranklist
ranks = read.table(rnk.file, header=T, colClasses = c("character", "numeric"))
ranks = setNames(ranks$Rank, ranks$GeneID)
  # can also in this form: tmp=ranks$Rank; names(tmp)=ranks$GeneID

# extract gene sets
pathways = gmtPathways(gmt.file)

# check your ranklist and geneset
str(ranks)
str(head(pathways))

# run user-defined GSEA
fgseaRes = fgsea(pathways, ranks, nperm = 10, minSize=0, maxSize=500)
head(fgseaRes[order(pval), ])

# plot for one pathway
plotEnrichment(pathways[["Gene_set_3"]],ranks) + labs(title="Gene_set_3")


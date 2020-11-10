#Load packages
library(limma)
library(edgeR)
library(RColorBrewer)

#Load data
files=list.files(pattern='txt$')

#Read and Merge a Set of Files Containing Count Data
x=readDGE(files,columns=c(1,2))
class(x)
dim(x)
x

#Add group labels
group=as.factor(c('EmBPOE4','EmBPOE4','EmBPOE4','EmBPOE4','WT','WT','WT','WT'))
x$samples$group=group

#Data pre-processing
#Filterx$genes=rownames(x)ing
table(rowSums(x$counts==0)==8) #Genes with zero counts in all samples
keep.exprs=filterByExpr(x,group=group) #Remove genes with very low counts across all groups
x=x[keep.exprs,,keep.lib.sizes=F]
dim(x)

#Add gene ID
x$genes=rownames(x)

#Normalisation: Remove systematic bias due to technical effects
x=calcNormFactors(x,method='TMM')
x$samples$norm.factors

#Data visualisation
#Boxplots
lcpm=cpm(x,log=T) #log Counts per Million
boxplot(lcpm,las=2,col='grey',main="Distribution of samples",ylab="Log-cpm") 

#MDS plots: Multidimensional scaling plot
plotMDS(lcpm)
plotMDS(lcpm,labels=group)


#Differential expression analysis by limma
#Design matrix
design=model.matrix(~0+group)
colnames(design)=gsub('group','',colnames(design))
design #Effects are estimated for parameters specified by the design matrix

#Contrasts
contr.matrix=makeContrasts(EmBPOE4vsWT=EmBPOE4-WT,levels=colnames(design))
contr.matrix

#Estimate mean-variance relationship
v=voom(x,design,plot=T) #voom: Use mean-var trend to assign weights to each observation
v #Weights got from voom used in linear modelling removes mean-variance trend

#Fit linear model
vfit=lmFit(v,design) #Fits a linear model for each gene
vfit=contrasts.fit(vfit,contrast=contr.matrix) #Computes statistics and fold changes for comparisons of interest
efit=eBayes(vfit) #Computes more precise estimates by sharing gene information using empirical Bayes moderation
plotSA(efit,main='Final model: Mean-variance trend') #Plots residual standard deviation versus average log expression for fitted model

#Test for Differential Expression(DE)
summary(decideTests(efit)) #adjusted p-value < 0.05

#Top 10 genes from DE analysis
topTable(efit,coef='EmBPOE4vsWT')
 #LogFC: Log fold change, Fold change = ratio between two values
 #AveExpr: Average gene expression across all samples


##Test for DE relative to a threshold: when number of DE genes is large
tfit=treat(vfit,lfc=1) #threshold: logFC >= 1
dt=decideTests(tfit)
summary(dt)

vennDiagram(dt[,1],circle.col=c("turquoise","salmon"))
write.fit(tfit,dt,file="results.txt")

##Top genes from DE analysis
EmBPOE4.vs.WT=topTreat(tfit,coef='EmBPOE4vsWT',n=Inf)
head(EmBPOE4.vs.WT)

#Mean-difference plot
plotMD(tfit,column=1,status=dt[,1],main=colnames(tfit)[1],xlim=c(-8,13))

#Interactive version of MD plot using the Glimma package
library(Glimma)
glMDSPlot(lcpm,labels=group,groups=x$samples[,c(2,4)],launch=T)
glMDPlot(tfit,coef=1,status=dt,main=colnames(tfit)[1],counts=lcpm,groups=group)

#Heatmap
library(gplots)
EmBPOE4.vs.WT.topgenes=EmBPOE4.vs.WT$ID[1:100]
i=which(v$genes %in% EmBPOE4.vs.WT.topgenes)
mycol=colorpanel(1000,"blue","white","red")
heatmap.2(lcpm[i,],scale="row",labRow=v$genes[i],labCol=group,col=mycol,trace="none",density.info="none",margin=c(8,6),lhei=c(2,10),dendrogram="column")


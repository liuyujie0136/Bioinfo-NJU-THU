### Co-expression Network
### Weighted Correlation Network Analysis (WGCNA): 加权基因共表达网络分析

## Install packages and library
# BiocManager::install(c("AnnotationDbi","impute","GO.db","preprocessCore","multtest"))
install.packages("WGCNA")  # also need "stringr", "reshape2"
library(WGCNA)


## Import data
datExpr <- readRDS(file="input/input_fpkm_matrix.rds")
datTraits <- readRDS(file="input/data_traits.rds")

# check
datExpr[1:4,1:4]  # each row represents a cell line (sample), each column represents the fpkm value of a gene
dim(datExpr)  # 56 cell lines (samples), 5000 genes
datTraits[1:4,]  # each row represents a cell line (sample), each column provides the gsm number, the cell line name and the cell line subtype information
dim(datTraits)


## Pick the soft thresholding power
options(stringsAsFactors = FALSE)
enableWGCNAThreads()  #open multithreading

powers = c(c(1:10), seq(from=12, to=20, by=2))
# Call the network topology analysis function，choose a soft-threshold to fit a scale-free topology to the network
sft=pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plot the results
pdf(file="output/soft_thresholding.pdf",width=9, height=5)
par(mfrow = c(1,2))
cex1 = 0.9

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

# Red line corresponds to using an R^2 cut-off
abline(h=0.90,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

dev.off()# close and generate a pdf of two plots

# 软阈值（即权重参数，可以理解为相关系数的β次幂）取值默认为1到30，上述图形的横轴均代表软阈值，左图的纵轴数值越大，说明该网络越逼近无尺度网络，右图的纵轴表示对应的基因模块中所有基因邻接性的均值。

sft$powerEstimate


## One-step network construction and module detection
net = blockwiseModules(datExpr,
                       power = sft$powerEstimate,
                       maxBlockSize = 6000,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "FPKM-TOM",
                       verbose = 3)

table(net$colors)  # show the total modules and genes in each modules. The '0' means genes do not belong to any module


## Module visualization
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
table(mergedColors)

# Plot the dendrogram and the module colors underneath
pdf(file="output/module_visualization.pdf",width=9, height=5)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# assign all of the gene to their corresponding module hclust
dev.off()


## Quantify module similarity by eigengene correlation
# Eigengene: One of a set of right singular vectors of a gene's x samples matrix that tabulates, e.g., the mRNA or gene expression of the genes across the samples.

# Recalculate module eigengenes
moduleColors = labels2colors(net$colors)
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes

# Add the weight to existing module eigengenes
MET = orderMEs(MEs)

# Plot the relationships between the eigengenes and the trait
pdf(file="output/eigengenes_trait_relationship.pdf",width=7, height=9)
par(cex = 0.9)
plotEigengeneNetworks(MET,"", marDendro=c(0,4,1,2), 
                      marHeatmap=c(3,4,1,2), cex.lab=0.8, xLabelsAngle=90)
dev.off()


## Find the relationships between modules and traits
design = model.matrix(~0+ datTraits$subtype)
colnames(design) = levels(datTraits$subtype)
moduleColors = labels2colors(net$colors)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, design, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

pdf(file="output/module_trait_relationship.pdf",width=9, height=10)

# Display the correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))

# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(design),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.6,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()


## Select specific module (We choose the "brown" module in trait “Luminal” for further analyses)

## 1) Intramodular connectivity, module membership, and screening for intramodular hub genes
# Intramodular connectivity
connet=abs(cor(datExpr,use="p"))^6
Alldegrees1=intramodularConnectivity(connet, moduleColors)

# Relationship between gene significance and intramodular connectivity
which.module="brown"
Luminal= as.data.frame(design[,3])
names(Luminal) = "Luminal"
GS1 = as.numeric(cor(datExpr, Luminal, use = "p"))
GeneSignificance=abs(GS1)

# Generalizing intramodular connectivity for all genes on the array
datKME=signedKME(datExpr, MEs, outputColumnName="MM.")

# Display the first few rows of the data frame
head(datKME)

# Finding genes with high gene significance and high intramodular connectivity in specific modules
#abs(GS1)>.8 #adjust parameter based on actual situations
#abs(datKME$MM.black)>.8 #at least larger than 0.8
FilterGenes= abs(GS1)>0.8 & abs(datKME$MM.brown)>0.8
table(FilterGenes)

# Find 3 hub genes
hubgenes <- rownames(datKME)[FilterGenes]
hubgenes


## 2) Export the network
# Recalculate topological overlap
TOM = TOMsimilarityFromExpr(datExpr, power = 6) 
# Select module
module = "brown"
# Select module probes
probes = colnames(datExpr)
inModule = (moduleColors==module)
modProbes = probes[inModule] 
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)

# Export the network into edge and node list files Cytoscape can read
#default threshold = 0.5, we could adjust parameter based on actual situations or in Cytoscape
cyt = exportNetworkToCytoscape(modTOM,
    edgeFile = paste("CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
    nodeFile = paste("CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
    weighted = TRUE,
    threshold = 0.02,
    nodeNames = modProbes, 
    nodeAttr = moduleColors[inModule])

# Screen the top genes
nTop = 10
IMConn = softConnectivity(datExpr[, modProbes])
top = (rank(-IMConn) <= nTop)
filter <- modTOM[top, top]

cyt = exportNetworkToCytoscape(filter,
    edgeFile = paste("output/CytoscapeInput-edges-filter-", paste(module, collapse="-"), ".txt", sep=""),
    nodeFile = paste("output/CytoscapeInput-nodes-filter-", paste(module, collapse="-"), ".txt", sep=""),
    weighted = TRUE,
    threshold = 0.02,
    nodeNames = rownames(filter), 
    nodeAttr = moduleColors[inModule][1:nTop])
# we can visualize it using Cytoscape


## 3) Extract gene IDs in specific module
# Select module
module = "brown"
# Select module probes (gene ID)
probes = colnames(datExpr)
inModule = (moduleColors == module)
modProbes = probes[inModule]
write.table(modProbes,file="output/geneID_brown.txt",sep="\t",quote=F,row.names=F,col.names=F)
# We could use the gene ID list for GO/KEGG analysis


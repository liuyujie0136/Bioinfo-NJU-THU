library(DESeq2)

x <-
  read.table("GSE143259_raw_counts.txt",
             header = TRUE,
             row.names = 1)

x_lab <-
  data.frame(row.names = colnames(x), condition = c(rep("WT", 4), rep("EX", 4)))

dds <-
  DESeqDataSetFromMatrix(countData = x,
                         colData = x_lab,
                         design =  ~ condition)

dds <- DESeq(dds)
res <- results(dds)

##另一种检验，LRT，似然率检验
#ddsLRT <- DESeq(dds,test="LRT",reduced=~1)
#resLRT <- results(ddsLRT)

resdata <-
  merge(as.data.frame(res),
        as.data.frame(counts(dds, normalized = TRUE)),
        by = "row.names",
        sort = FALSE)


up_down <-
  ifelse(
    !is.na(resdata$log2FoldChange) & !is.na(resdata$padj) &
      abs(resdata$log2FoldChange) > 1 & resdata$padj < 0.05,
    ifelse(resdata$log2FoldChange > 1,
           1,-1),
    0
  )

top <- resdata[!is.na(resdata$padj) & !is.na(resdata$log2FoldChange), ]

top <- top[top$padj < 0.05 & abs(top$log2FoldChange) > 1, 1:7]

ann <- data.frame(GeneID = rownames(x))
Glimma::glMDPlot(dds, anno = ann, groups = x_lab[[1]], status = up_down)



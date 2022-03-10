f <- read.table("heatmap-test.txt", header = TRUE, row.names = 1)

pheatmap::pheatmap(f)

stats::heatmap(as.matrix(f))

gplots::heatmap.2(as.matrix(f))

library(ggplot2)
f2 <- cbind(GeneID = rownames(f), f)
f3 <-
  reshape2::melt(f2, id.vars = "GeneID", variable.name = "Sample")

ggplot(f3) +
  geom_tile(aes(x = Sample, y = GeneID, fill = value)) +
  scale_fill_gradient(low = "white", high = "steelblue")

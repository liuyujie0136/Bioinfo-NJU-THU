library(clusterProfiler)


g <- read.table("genes.txt", sep = "\t", header = FALSE)
gene = as.character(g[, 1])


# library(biomaRt)
# listMarts()
# plant<-useMart("ensembl")
# listDatasets(plant)
# mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
# listFilters(mart)


ids <- bitr(gene,fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = "org.Hs.eg.db")
head(ids)
genes = ids[,2]
head(genes)
ego <- enrichGO(gene=genes,'org.Hs.eg.db',
                ont="BP",
                pvalueCutoff=0.01,
                readable=TRUE)
enrich_gobp <- as.data.frame(ego)
write.csv(enrich_gobp,"enrich_gobp.csv",row.names = F)
barplot(ego,showCategory = 15,title = "EnrichmentGO_BP")

egCC <- enrichGO(gene=genes,'org.Hs.eg.db',
                 ont="CC",
                 pvalueCutoff=0.01,
                 readable=TRUE)
enrich_gocc <- as.data.frame(egCC)
write.csv(enrich_gocc,"enrich_gocc.csv",row.names = F)
barplot(egCC,showCategory = 15,title = "EnrichmentGO_CC")

egMF <- enrichGO(gene=genes,'org.Hs.eg.db',
                 ont="MF",
                 pvalueCutoff=0.01,
                 readable=TRUE)
enrich_gomf <- as.data.frame(egMF)
write.csv(enrich_gomf,"enrich_gomf.csv",row.names = F)
barplot(egMF,showCategory = 15,title = "EnrichmentGO_MF")

ekk <- enrichKEGG(gene = genes,
                  organism = 'hsa',
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,)
ekk <- setReadable(ekk,OrgDb = org.Hs.eg.db,keytype = "ENTREZID")
enrich_kegg <- as.data.frame(ekk)
head(enrich_kegg)
write.csv(enrich_kegg,"enrich_kegg.csv",row.names = F)
barplot(ekk,showCategory = 15,title = "Enrichment_KEGG")
dotplot(ekk,title = "Enrichment_KEGG")
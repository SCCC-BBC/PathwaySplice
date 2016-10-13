AnalysisGO <- function() {
  # library(org.Hs.eg.db)
  # data(geneList)
  # gene <- names(geneList)[abs(geneList) > 2]
  # gene.df <- bitr(gene, fromType = "ENTREZID", 
  #                 toType = c("ENSEMBL", "SYMBOL"),
  #                 OrgDb = org.Hs.eg.db)
  # 
  # ego <- enrichGO(gene          = gene,
  #                 universe      = names(geneList),
  #                 OrgDb         = org.Hs.eg.db,
  #                 ont           = "CC",
  #                 pAdjustMethod = "BH",
  #                 pvalueCutoff  = 0.01,
  #                 qvalueCutoff  = 0.05)
  # 
  # head(summary(ego))
  # 
  gene2<-row.names(Example.Go.adjusted.by.exon[[2]][which(Example.Go.adjusted.by.exon[[2]][,1]==1),])
  
  gene.df <- bitr(gene2, fromType = "ENSEMBL", 
                  toType = c("ENTREZID", "SYMBOL"),
                  OrgDb = org.Hs.eg.db)
  
  ego2 <- enrichGO(gene         = gene2,
                   OrgDb         = org.Hs.eg.db,
                   keytype       = 'ENSEMBL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.5,
                   qvalueCutoff  = 0.5)
  
  head(summary(ego2))
  dotplot(ego2, showCategory=30)
  
  ego3 <- enrichGO(gene         = gene.df$SYMBOL,
                   OrgDb         = org.Hs.eg.db,
                   keytype       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.5,
                   qvalueCutoff  = 0.5)
  
  head(summary(ego3))
  enrichMap(ego3, vertex.label.cex=1.2, layout=igraph::layout.kamada.kawai)
  cnetplot(ego3)
  plotGOgraph(ego2)
}

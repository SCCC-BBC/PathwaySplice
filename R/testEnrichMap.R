testEnrichMap <- function() {
  
  library(org.Hs.eg.db)
  library(PathwaySplice)
  library(goseq)
  library(gdata)
  library(GO.db)
  library(DOSE,quietly = TRUE)
  #library(clusterProfiler,quietly = TRUE)
  library(reshape2)
  library(igraph)
  
  #data(package="PathwaySplice")
  #data(geneList)
  data(mds)
  
  #Analysis
  Example.Go.adjusted.by.exon<-Run_pathwaysplice(mds,ad="exon_SJ",sub_feature="E",0.05,genomeID="hg19",geneID="ensGene",gene_model=hg19.gene.model,method="Wallenius")
  
  Example.Go.unadjusted<-Run_pathwaysplice(mds,ad="exon_SJ",sub_feature="E",0.05,genomeID="hg19",geneID="ensGene",gene_model=hg19.gene.model,method="Hypergeometric")
  
  #Check bias
  plotPWF2(Example.Go.adjusted.by.exon[[2]])
  #gene <- names(geneList)[abs(geneList) > 2]
  
  re.w.adjusted<-enMap2(Example.Go.adjusted.by.exon,n=4,SimilarityThreshold=0)
  
  re.w.unadjusted<-enMap2(Example.Go.unadjusted,n=4,SimilarityThreshold=0)
  
  # DE.index<-which(Example.Go.adjusted.by.exon[[2]]$DEgenes == 1)
  # gene<-unique(rownames(Example.Go.adjusted.by.exon[[2]][DE.index,]))
  # 
  # gene.df <- bitr(gene, fromType = "ENSEMBL", 
  #                 toType = c("ENTREZID", "SYMBOL"),
  #                 OrgDb = org.Hs.eg.db)
  # 
  # head(gene.df)
  # 
  #sink("test_ok_HT.txt")
  #Use modified enrichGO
  
  # ego.HT <- enrichGO(gene          = gene,
  #                    universe      = rownames(Example.Go.adjusted.by.exon[[2]]),
  #                    keytype = "ENSEMBL",
  #                    OrgDb         = org.Hs.eg.db,
  #                    ont           = "BP",
  #                    pAdjustMethod = "BH",
  #                    pvalueCutoff  = 1,
  #                    qvalueCutoff  = 0.8,GOFromGOSeq=Example.Go.adjusted.by.exon,whichway = "HT")
  # dim(ego.HT@result)
  # 
  #   re.HT<-enMap(ego.HT, vertex.label.cex=1,layout=igraph::layout.kamada.kawai,SimilarityThreshold=0.1)
  # 
  #sink()

  #sink("test_ok_WHT.txt")
  
  # ego.WHT <- enrichGO(gene          = gene,
  #                     universe      = rownames(Example.Go.adjusted.by.exon[[2]]),
  #                     keytype = "ENSEMBL",
  #                     OrgDb         = org.Hs.eg.db,
  #                     ont           = "BP",
  #                     pAdjustMethod = "BH",
  #                     pvalueCutoff  = 1,
  #                     qvalueCutoff  = 0.8,GOFromGOSeq=Example.Go.adjusted.by.exon,whichway = "WHT")
  # 
  # dim(ego.WHT@result)
  # #sink()
  # 
  # par(mfrow=c(2,1))
  # #enrichMap(ego.HT, vertex.label.cex=1.2, layout=igraph::layout.kamada.kawai)
  # #enrichMap(ego.WHT, vertex.label.cex=1.2, layout=igraph::layout.kamada.kawai)
  # 
  # re.WHT<-enMap(ego.WHT, vertex.label.cex=1, layout=igraph::layout.kamada.kawai,SimilarityThreshold=0.1)
  # 
  # re.w1<-enMap2(Example.Go.adjusted.by.exon,n=5,SimilarityThreshold=0)
  # re.w2<-enMap2(Example.Go.unadjusted,n=5,SimilarityThreshold=0)
}
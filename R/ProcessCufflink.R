#' ProcessCuffLinkResults
#'
#' @return
#' @export
#'
#' @examples
ProcessCuffLinkResults <- function() {
  getwd()

  cuff_data<-readCufflinks()

  diffGeneIDs <- getSig(cuff_data,level="genes",alpha=0.05)
  diffGenes<-getGenes(cuff_data,diffGeneIDs)

  names<-featureNames(diffGenes)
  row.names(names)=names$tracking_id
  diffGenesNames<-as.matrix(names)
  diffGenesNames<-diffGenesNames[,-1]

  names<-featureNames(diffGenes)
  row.names(names)=names$tracking_id
  diffGenesNames<-as.matrix(names)
  diffGenesNames<-diffGenesNames[,-1]


  diffGenesData<-diffData(diffGenes)
  row.names(diffGenesData)=diffGenesData$gene_id
  diffGenesData<-diffGenesData[,-1]

  # merge the two matrices by row names
  diffGenesOutput<-merge(diffGenesNames,diffGenesData,by="row.names")

  head(diffGenesOutput)

  featureNames(diffGenes)

  csDensity(genes(cuff_data))
  csScatter(genes(cuff_data),'Mut','HC')
  csVolcano(genes(cuff_data),'Mut','HC')
  mygene<-getGene(cuff_data, 'XLOC_011342')
  expressionBarplot(mygene)

  gene_diff_data<-diffData(genes(cuff_data))
  sig_gene_data<-subset(gene_diff_data,(significant=='yes'))
  nrow(sig_gene_data)

  isoform_diff_data<-diffData(isoforms(cuff_data), 'Mut', 'HC')
  sig_isoform_data<-subset(isoform_diff_data, (significant=='yes'))
  nrow(sig_isoform_data)

  tss_diff_data<-diffData(TSS(cuff_data), 'Mut', 'HC')
  sig_tss_data<-subset(tss_diff_data, (significant=='yes'))
  nrow(sig_tss_data)

  cds_diff_data<-diffData(CDS(cuff_data), 'Mut', 'HC')
  sig_cds_data<-subset(cds_diff_data, (significant=='yes'))
  nrow(sig_cds_data)

  promoter_diff_data<-distValues(promoters(cuff_data))
  sig_promoter_data<-subset(promoter_diff_data, (significant=='yes'))
  nrow(sig_promoter_data)

  splicing_diff_data<-distValues(splicing(cuff_data))
  sig_splicing_data<-subset(splicing_diff_data, (significant=='yes'))
  nrow(sig_splicing_data)

  relCDS_diff_data<-distValues(relCDS(cuff_data))
  sig_relCDS_data<-subset(relCDS_diff_data, (significant=='yes'))
  nrow(sig_relCDS_data)

  gene_diff_data<-diffData(genes(cuff_data))
  sig_gene_data<-subset(gene_diff_data, (significant=='yes'))
  write.table(sig_gene_data, "diff_genes.txt", sep = "\t", row.names = F,col.names = T, quote = F)

  dim(gene_diff_data)

  data.isoform<-read.table("isoform_exp.diff",header=T)
}

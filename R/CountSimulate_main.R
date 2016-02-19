# test count simulator
# use NBsim code from Xiaobei,
#source("http://bioconductor.org/biocLite.R")
biocLite("goseq")
library(tweeDEseqCountData)
library(edgeR)
library(goseq)
data(pickrell)
pickrell = as.matrix(exprs(pickrell.eset))

####################################
# simulation scenario 1
#  1. simulate DE and counts
#  2. unknown genes and pathways
#  3. randomly assign junctions to genes
################################

nSamp = 10
nGenes = 5000
nJunction = 30000
nPway = 1000
Pway_max_size = 100
Pway_min_size = 15
grp = as.factor(rep(0:1, each=nSamp/2))
data.sim = NBsim(foldDiff = 3, dataset = pickrell, nTags = nJunction, group = grp, verbose = TRUE, 
                       add.outlier = F, outlierMech = "S", pOutlier = 0.2, drop.extreme.dispersion = 0.1
)

data.sim.count = data.sim$counts
colnames(data.sim.count) = c(paste0("Control_", 1:(nSamp/2)), paste0("Treated_", 1:(nSamp/2)))

# assign features to genes
# for example, 1000 features assigned to 100 genes
dexp.tmp = dexp(x=seq(from=0,to=8,length.out=nGenes), rate=1)
dexp.prob = dexp.tmp/sum(dexp.tmp)
assigned.id = sample(1:nGenes, nJunction, replace=T, prob=dexp.prob)
length(table(assigned.id))

# use DESeq diff test
DESeq_2group = function(Data, Sample_names1, Sample_names2){
  cat(paste("Comparing", paste(Sample_names1,collapse = "+"), "(Treated)\n\tto", paste(Sample_names2,collapse = "+"), "(Untreated)\n"))
  N_data = nrow(Data)
  nsample1 = length(Sample_names1)
  nsample2 = length(Sample_names2)
  colID_select = match(c(Sample_names1,Sample_names2), colnames(Data))
  Data_s = Data[, colID_select]
  condition = factor(c(rep("Treated",nsample1),rep("Untreated",nsample2)))
  dds = DESeqDataSetFromMatrix(Data_s, DataFrame(condition), ~ condition)
  dds$condition = relevel(dds$condition, ref = "Untreated")
  return(results(DESeq(dds)))
}

data.sim.DE = data.frame(DESeq_2group(data.sim.count, paste0("Control_", 1:(nSamp/2)), paste0("Treated_", 1:(nSamp/2))))

# match DE juncitons with genes assigned the junction
gene.uid = unique(assigned.id)
gene.DE = matrix(0,nrow=length(gene.uid), ncol=3)
for (i in 1:length(gene.uid)){
  junc.sim.id = which(assigned.id==gene.uid[i])
  Njunc = length(junc.sim.id)
  Njunc_DE = length(which(data.sim.DE$padj[junc.sim.id] < 0.05))
  if(Njunc_DE>0){
    DE = 1
  }else{
    DE = 0
  }
  gene.DE[i,] = t(c(DE, Njunc_DE, Njunc))
}
rownames(gene.DE) = paste0("gene_", gene.uid)
colnames(gene.DE) = c("DE", "Njunc_DE", "Njunc")
gene.DE = data.frame(gene.DE)

# plot bias evidence 
gene.pwf = makespline(gene.DE[,3], gene.DE[,1])
gene.DE.pwf = data.frame(DEgenes = gene.DE[,1], bias.data = gene.DE[,3], pwf = gene.pwf)
plotPWF2(gene.DE.pwf, binsize=30)

# pathway selection and enrichment scores
goseq_wallenius = function(pwf, gene2cat){
  cat2gene = reversemapping(gene2cat)
  cats = names(cat2gene)
  DE = rownames(pwf)[pwf$DEgenes == 1]
  num_de = length(DE)
  num_genes = nrow(pwf)
  pvals = data.frame(category = cats, over_represented_pvalue = NA, 
                     under_represented_pvalue = NA, stringsAsFactors = FALSE, 
                     numDEInCat = NA, numInCat = NA)
  message("Calculating the p-values...")
  degenesnum = which(pwf$DEgenes == 1)
  cat2genenum = relist(match(unlist(cat2gene), rownames(pwf)), 
                       cat2gene)
  alpha = sum(pwf$pwf)
  pvals[, 2:3] = t(sapply(cat2genenum, function(u) {
    num_de_incat = sum(degenesnum %in% u)
    num_incat = length(u)
    avg_weight = mean(pwf$pwf[u])
    weight = (avg_weight * (num_genes - num_incat))/(alpha - num_incat * avg_weight)
    if (num_incat == num_genes) {
      weight = 1
    }
    c(dWNCHypergeo(num_de_incat, num_incat, num_genes - num_incat, num_de, weight) + pWNCHypergeo(
      num_de_incat, num_incat, num_genes - num_incat, num_de, weight, lower.tail = FALSE), 
      pWNCHypergeo(num_de_incat, num_incat, num_genes - num_incat, num_de, weight))
  }))
  degenesnum = which(pwf$DEgenes == 1)
  cat2genenum = relist(match(unlist(cat2gene), rownames(pwf)), 
                       cat2gene)
  pvals[, 4:5] = t(sapply(cat2genenum, function(u) {
    c(sum(degenesnum %in% u), length(u))
  }))
  pvals = pvals[order(pvals$over_represented_pvalue), ]
  pvals
}
# from goseq
reversemapping=function(map){
  tmp=unlist(map,use.names=FALSE)
  names(tmp)=rep(names(map),times=as.numeric(summary(map)[,1]))
  return(split(names(tmp),as.vector(tmp)))
}

# randomly create a pathway datasets
pway2gene = lapply(1:nPway,function(h)rownames(gene.DE)[sample(1:nrow(gene.DE), 
                                                               sample(Pway_min_size:Pway_max_size,1))])
names(pway2gene) = paste0("Pathway_", 1:nPway)
gene2pway = reversemapping(pway2gene)

gene.DE.nopwf = gene.DE.pwf
gene.DE.nopwf$pwf = rep(1/nrow(gene.DE.nopwf), nrow(gene.DE.nopwf))
pway.enrich = goseq_wallenius(gene.DE.nopwf, gene2pway)

# plot Avg.enrich.pval vs % genes with #junctions>10
#  combine the DE junction dataset
pway.DE.perc = rep(0, nrow(pway.enrich))
pway.DE.avgNjunc = rep(0, nrow(pway.enrich))
for(i in 1:nrow(pway.enrich)){
  pway.DE.gid = match(pway2gene[[pway.enrich$category[i]]], rownames(gene.DE))
  pway.DE.avgNjunc[i] = sum(gene.DE$Njunc[pway.DE.gid])/length(pway.DE.gid)
  pway.DE.perc[i] = length(which(gene.DE$Njunc[pway.DE.gid]>20))/length(pway.DE.gid)
}
pway.enrich$avg.Njunc = pway.DE.avgNjunc
pway.enrich$perc.Njunc = pway.DE.perc

# plots show bias
plot(y=-log(pway.enrich$over_represented_pvalue), x=pway.enrich$avg.Njunc, xlab = "Number of junctions per gene")
plot(y=-log(pway.enrich$over_represented_pvalue), x=pway.enrich$perc.Njunc, xlab="% genes with > 20 junctions")

# check
gene.DE[match(pway2gene[["Pathway_48"]], rownames(gene.DE)),]


#####################################3
# simulation scenario 2
#  1. known genes and pathways
#  2. assign junctions to genes with known database
#  3. only simulate DE and counts
###################################



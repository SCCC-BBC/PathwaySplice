# test count simulator
# use NBsim code from Xiaobei,
#source("http://bioconductor.org/biocLite.R")
biocLite("goseq")
library(tweeDEseqCountData)
library(edgeR)
library(goseq)
library("DESeq2")
source("/home/jamesban/Desktop/Projects/Research/Junctionbias/CountSimulator/AdaptedCodes_juncbias.R")
data(pickrell)

pickrell = as.matrix(exprs(pickrell.eset))
require(parallel) 
library("ggplot2")

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
rownames(gene.DE.pwf) = paste0("gene_",rownames(gene.DE.pwf))
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

goseq_sampling = function(pwf, gene2cat){
  cat2gene = reversemapping(gene2cat)
  cats = names(cat2gene)
  DE = rownames(pwf)[pwf$DEgenes == 1]
  num_de = length(DE)
  num_genes = nrow(pwf)
  pvals = data.frame(category = cats, over_represented_pvalue = NA, 
                     under_represented_pvalue = NA, stringsAsFactors = FALSE, 
                     numDEInCat = NA, numInCat = NA)
  num_DE_mask = rep(0, length(cats))
  a = table(unlist(gene2cat[DE], FALSE, FALSE))
  num_DE_mask[match(names(a), cats)] = as.numeric(a)
  num_DE_mask = as.integer(num_DE_mask)
  gene2cat = gene2cat[rownames(pwf)]
  names(gene2cat) = rownames(pwf)
  message("Running the simulation...")
  lookup = matrix(0, nrow = repcnt, ncol = length(cats))
  for (i in 1:repcnt) {
    a = table(as.character(unlist(gene2cat[order(runif(num_genes)^(1/pwf$pwf), 
                                                 decreasing = TRUE)[1:num_de]], FALSE, FALSE)))
    lookup[i, match(names(a), cats)] = a
    pp(repcnt)
  }
  message("Calculating the p-values...")
  pvals[, 2] = (colSums(lookup >= outer(rep(1, repcnt), 
                                        num_DE_mask)) + 1)/(repcnt + 1)
  pvals[, 3] = (colSums(lookup <= outer(rep(1, repcnt), 
                                        num_DE_mask)) + 1)/(repcnt + 1)
  
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
pway.enrich = goseq(gene.DE.nopwf, gene2cat=gene2pway, method="Hypergeometric")

# correcting bias using the PWF from goseq
pway.enrich = goseq(gene.DE.pwf, gene2cat=gene2pway, method="Wallenius")

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
cor.test(-log(pway.enrich$over_represented_pvalue), pway.enrich$perc.Njunc, method="s")
cor.test(-log(pway.enrich$over_represented_pvalue), pway.enrich$avg.Njunc, method="s")

# check
gene.DE[match(pway2gene[["Pathway_48"]], rownames(gene.DE)),]


#####################################3
# simulation scenario 2
#  1. known genes and pathways
#  2. assign junctions to genes with known database
#  3. only simulate DE with sampling
###################################
dir_juncs = "/home/jamesban/Desktop/Projects/Research/Junctionbias/CountSimulator/"
data.gyjuncs = read.csv(paste0(dir_juncs, "6_Samples_Count_data.csv"), row.names = 1)

# running paralell
require(parallel) 
Ncores = getOption("mc.cores", 10)

# randomly select junctions from all junctions
junc.id = grep(":J",rownames(data.gyjuncs))
data.gyjuncs.only = data.gyjuncs[junc.id, ]

# collapse to gene level
juncID2gene.list = (strsplit(gsub(":.*","",rownames(data.gyjuncs.only)), split="\\+"))
genes.wjunc = unique(unlist(juncID2gene.list))

# set parallel bins
bins = findInterval(1:length(genes.wjunc),seq(from=1, to=length(genes.wjunc)+1, length.out = Ncores+1))

# this will take a while
gene2juncID.list = mclapply(1:Ncores, function(i){
  lapply(genes.wjunc[which(bins==i)], function(h)grep(h, rownames(data.gyjuncs.only)))
}, mc.cores = Ncores)

gene2juncID.list = unlist(gene2juncID.list, recursive = F)
genes.nJunc = unlist(lapply(gene2juncID.list, function(h)length(h)))

# test with a gene match multiple combinations **result: gene maps to only single combination**

#genes.unique.test = mclapply(1:Ncores, function(i){
#  unlist(lapply(genes.wjunc[which(bins==i)], function(h){
#    tmp.name = rownames(data.gyjuncs.only)[grep(h, rownames(data.gyjuncs.only))]
#    tmp.name = gsub(":.*","",tmp.name)
#    length(unique(tmp.name))
#  }))
#}, mc.cores = Ncores)
#tmp = unlist(genes.unique.test)

# randomly select a set of junctions and assign them DE
#  set DE criteria that if a juncion is DE then the gene is DE
junc.DE = rep(0, nrow(data.gyjuncs.only))
pwf.sim.gy = data.frame(DEgenes=0, bias.data=genes.nJunc, pwf=rep(1/length(genes.wjunc),length(genes.wjunc)))
rownames(pwf.sim.gy) = genes.wjunc

# prepare the gene ontology
# run goseq **if wish to use sampling method, add repcnt=1000**
gene2go = getgo(genes = rownames(pwf.sim.gy), genome = "mm10", id="ensGene", fetch.cats = c("GO:BP"))

# simulation
GO.top10 = character(0)
nJunc.DE = 2000
nSim = 100
pwf.sim.gy.list = list()

for (j in 1:nSim){
  cat("create simulation data at i = ", j, "\n")
  junc.id.rand = sample(1:nrow(data.gyjuncs.only),nJunc.DE)
  junc.DE[junc.id.rand] = 1
  genes.DE = unlist(mclapply(1:Ncores, function(i){
    unlist(lapply(which(bins==i), function(h){
      if(sum(junc.DE[gene2juncID.list[[h]]])>0){
        s=1
      }else{
        s=0
      }
      s
    }))
  }, mc.cores = Ncores))
  pwf.sim.gy$DEgenes = genes.DE
  pwf.sim.gy.list[[j]] = pwf.sim.gy
}

# run goseq on each dataset
bins.sim = findInterval(1:nSim,seq(from=1, to=nSim+1, length.out = Ncores+1))
res.tmp = list()

# this will take a while
res.sim100 = unlist(mclapply(1:Ncores, function(i){
  res.tmp[[i]] = as.character()
  for(j in which(bins.sim==i)){
    cat("calculating for list ", j, "\n")
    GO.samp = goseq(pwf.sim.gy.list[[j]], gene2cat = gene2go, use_genes_without_cat=F, method="Hypergeometric")
    res.tmp[[i]] = c(res.tmp[[i]], GO.samp$term[1:10])
  }
  res.tmp
}, mc.cores = Ncores))

# barplot for percent occurring
GO.sim.per = table(res.sim100)/nSim
top10.GO.sim = GO.sim.per[order(GO.sim.per,decreasing=F)]

top10.GO.sim = data.frame(term = names(top10.GO.sim), perc = top10.GO.sim)
top10.GO.sim$term <- factor(top10.GO.sim$term, levels = top10.GO.sim$term)

pdf(file=paste0(dir_juncs,"/GO.BP.sim.100.pdf"), width=12, height=12)
gg = ggplot(top10.GO.sim, aes(x=term,y=perc))
gg + geom_bar(stat = "identity", width=0.7) + theme(axis.text.x = element_text(size=20), 
  axis.text.y = element_text(size=20),  panel.background =  element_rect(fill = "white", colour = NA), 
  panel.border = element_rect(fill = NA, colour="grey50"), 
  panel.grid.major =  element_line(colour = "grey90", size = 0.2),
  panel.grid.minor =  element_line(colour = "grey98", size = 0.5),
  panel.margin = unit(0.25, "lines")) + coord_flip() + scale_y_continuous(labels = scales::percent)
dev.off()


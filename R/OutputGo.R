
# test.goseq2.2.gene.symbol<-sapply(test.goseq2[,8],function(u,gene.model){
#    xx<-gene.model[match(u,as.character(gene.model[,3])),1]
#    names(xx)<-
#
#   },gene.model)
#
# head(test.goseq2.2.gene.symbol)
# dataset2<-test.goseq33
#
# dataset2[sapply(dataset2, is.list)] <-
#   sapply(dataset2[sapply(dataset2, is.list)],
#          function(x)sapply(x, function(y) paste(unlist(y),collapse=", ") ) )
#
# head(dataset2)
# write.table(dataset2,row.names = FALSE,"myfile.xls", quote=FALSE, sep="\t")

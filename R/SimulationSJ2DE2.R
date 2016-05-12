#' @title {Simulate gene sets based on the number of splicing junctions}
#
#' @decription{
#'We sample the mean of the number of splicing junction, and use this mean to simulate an array A including
#'splicing junction numbers with certain amounts. we convert this number into an array AAA between 0 and 1. Use the median of
#'AAA array as the probability for 1, we sample from [0,1] to get a list with certain length. This list is used as
#'a gene set. For example, if this gene set has 30 genes, each gene has its splicing junction number, we will
#'have 30 splicing junction number. After differential gene expression analysis, some gene is differential expressed,
#'and is labeled as 1,otherwise labeles as 0.
#'Here we simulate two scenarios:
#'Scenario1: the probability of being 1 is dependent on the median of AAA
#'Scenario2: the probability of being 1 is not dependent on the median of AAA}
#'
#'
#' @param min_num_splicing
#' @param max_num_splicing
#' @param num_gene
#'
#' @return
#' @export
#' @examples
#' re.random.DE2SJ<-SimulationSJ2DE(20,200,30)

SimulationSJ2DE<- function(min_num_splicing,max_num_splicing,num_gene) {
  crr.random<-array()
  for(j in 1:1:100){
    ans<-list()
    for(i in 1:50){
      mean_p<-sample(seq(min_num_splicing,max_num_splicing,10),1)
      temp<-data.frame()
      a=mean_p
      A=rpois(num_gene,a)
      AA=max(A)
      AAA=A/AA
      #AAA

      #p=median(AAA)
      p=runif(1)
      q=1-p
      SJNum<-A
      De.list<-sample(c(0,1),num_gene,replace = TRUE,prob=c(q,p))
      temp.data<-cbind(SJNum,De.list)

      De.prop<-length(which(De.list==1))/num_gene
      SJ.De.prop<-cbind(median(A),De.prop)

      ans[[i]]<-list(SJ.De.prop=SJ.De.prop,SJ.De.prop.data=temp.data)
    }

    ans.2<-do.call(rbind,ans)
    SJ.DE.prop<-matrix(unlist(ans.2[,1]),50,2,byrow=T)
    crr.random[j]<-cor(SJ.DE.prop)[1,2]

  }

  crr.force<-array()

  for(j in 1:100){
    ans<-list()
    for(i in 1:50){
      mean_p<-sample(seq(min_num_splicing,max_num_splicing,10),1)
      temp<-data.frame()
      a=mean_p
      A=rpois(num_gene,a)
      AA=max(A)
      AAA=A/AA
      #AAA

      p=median(AAA)
      #p=runif(1)
      q=1-p
      SJNum<-A
      De.list<-sample(c(0,1),num_gene,replace = TRUE,prob=c(q,p))
      temp.data<-cbind(SJNum,De.list)

      De.prop<-length(which(De.list==1))/30
      SJ.De.prop<-cbind(median(A),De.prop)

      ans[[i]]<-list(SJ.De.prop=SJ.De.prop,SJ.De.prop.data=temp.data)
    }

    ans.2<-do.call(rbind,ans)
    SJ.DE.prop<-matrix(unlist(ans.2[,1]),50,2,byrow=T)

    crr.force[j]<-cor(SJ.DE.prop)[1,2]

  }

  re<-list(Unrelated_SJ_DE=crr.random,related_SJ_DE=crr.force)
  re

}

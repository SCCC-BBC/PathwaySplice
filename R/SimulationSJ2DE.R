#' Title
#'
#' @param min_num_splicing
#' @param max_num_splicing
#' @param num_gene
#'
#' @return
#' @export
#'
#' @examples
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

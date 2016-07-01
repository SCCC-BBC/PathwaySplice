#' @ ProcessOutputFilesFrom_rMATS
#'
#' @ read output from rMAT
#'
#' @param input_file
#'
#' @return
#' @export
#'
#' @examples
#'
#' dir.name="/media/H_driver/PJ/rMAT/pHamard_MATS-selected-re-run-exon-centric/pHamard_MATS-selected-re-run-exon-centric/Neg_WTvKO_c0.01/"
#'
#' input.file.pattern="*Neg_WTvKO_c0.01.txt"
#'
#' sink(paste0(output.dir.name,"rMAT_output_from_PJ.txt"))
#' re.rMAT<-ProcessOutputFilesFrom_rMATS(dir.name,input.file.pattern)
#' sink()
#'
#'
ProcessOutputFilesFrom_rMATS<-function(dir.name,input.file.pattern){


  ProcessOutputFilesFrom_rMATS_read<-function(input_file){

    re=read.table(input_file,header=T)

    return(re)

  }

  file.name=paste0(dir.name,dir(dir.name,recursive = TRUE,pattern=input.file.pattern))
  file.name.2<-as.list(file.name)

  names(file.name.2)=sapply(strsplit(file.name,split="\\/"),"[[",9)

  re.out<-lapply(file.name.2,ProcessOutputFilesFrom_rMATS_read)

  re.out.2<-do.call(c,lapply(re.out, function(u){
    x<-as.character(u[which(u$FDR<0.05),]$GeneID)
    x
    #x<-as.data.frame(t(x))
    #colnames(x)<-colnames(Data4Goterm)
    #x
  }))

  print(names((re.out.2)))

  n1=length(re.out.2[grep("A3SS.MATS.ReadsOnTargetAndJunctionCounts",names(re.out.2))])
  n2=length(re.out.2[grep("A5SS.MATS.ReadsOnTargetAndJunctionCounts",names(re.out.2))])
  n3=length(re.out.2[grep("MXE.MATS.ReadsOnTargetAndJunctionCounts",names(re.out.2))])
  n4=length(re.out.2[grep("RI.MATS.ReadsOnTargetAndJunctionCounts",names(re.out.2))])
  n5=length(re.out.2[grep("SE.MATS.ReadsOnTargetAndJunctionCounts",names(re.out.2))])

  cat(n1,"\t",n2,"\t",n3,"\t",n4,"\t",n5,"\n")

  re.out.3<-list(JunctionCountOnly=unique(re.out.2[grep("JunctionCountOnly",names(re.out.2))]),
                 ReadsOnTargetAndJunctionCounts=unique(re.out.2[grep("ReadsOnTargetAndJunctionCounts",names(re.out.2))]),
                 SEMATSJunctionCountOnly=unique(re.out.2[grep("SE.MATS.JunctionCountOnly",names(re.out.2))]),
                 SEReadsOnTargetAndJunctionCounts=unique(re.out.2[grep("SE.MATS.ReadsOnTargetAndJunctionCounts",names(re.out.2))]),
                 SEs=unique(re.out.2[grep("SE.MATS",names(re.out.2))])
                 )

  return(re.out.3)

}

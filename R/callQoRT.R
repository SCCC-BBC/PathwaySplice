#' @ Use QoRTs to get the count of subfeatures in each gene
#'
#' @ cmd1
#' @ inputfile:bam file, gene annotation file,
#' @ outfile
#'
#' @return
#' @export
#'
#' @examples
#'
#' dir.name="/media/H_driver/PJ/"
#'
#' file.name=dir(dir.name,recursive = TRUE,pattern="sorted.bam")
#' file.name.whole<-paste0(dir.name,file.name)
#' file.name.selected<-file.name.whole[c(5,3,9)]
#' file.name.selected.2<-as.list(file.name.selected)
#' names(file.name.selected.2)=sapply(strsplit(file.name.selected,split="\\/"),"[[",6)
#'
#' cmd1="java -Xmx4000M -jar /home/aiminyan/QoRTs/QoRTsFullExampleData/QoRTsRelease/QoRTs.jar QC --noGzipOutput"
#' gtf1="/media/aiminyan/DATA/mus_musculus/Mus_musculus.GRCm38.83.processed.sorted.gtf"
#'
#' re.out<-lapply(file.name.selected.2,callQoRT,gtf_file=gtf1,runing_cmd=cmd1)
#'

callQoRT<-function(input_file,runing_cmd,gtf_file){

  inputfile=paste(input_file,gtf_file,sep=" ")

  outfile=paste("",sapply(strsplit(input_file,split="\\/"),"[[",2),sapply(strsplit(input_file,split="\\/"),"[[",3),
                sapply(strsplit(input_file,split="\\/"),"[[",4),sapply(strsplit(input_file,split="\\/"),"[[",6),sep="/")


  cmd2=paste(runing_cmd,inputfile,outfile,sep=" ")

  print(cmd2)


  system(cmd2, intern = TRUE, ignore.stderr = TRUE)

  return(cmd2)

}

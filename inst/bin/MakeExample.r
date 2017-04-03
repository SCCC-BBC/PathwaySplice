
#Generate a small annotation and count data set with 10 genes

cmd0 <- "awk '{print $10}' ~/PathwaySplice/inst/extdata/GTF_Files/flat.chr22.gff | sort -u"

gene.name <- try(system(cmd0,intern =TRUE))

n <- 10

cmd1 <- paste0("awk -v s1=","'",gene.name[1],"'")

cmd2 <- "'$10==s1 {print}'" 

cmd3 <- "~/PathwaySplice/inst/extdata/GTF_Files/flat.chr22.gff"

cmd4 <- paste0("> ~/PathwaySplice/inst/extdata/GTF_Files/flat.chr22.",n,".genes.gff")

cmd5 <- paste(cmd1,cmd2,cmd3,cmd4)

system(cmd5)
cmd6 <- paste0(">> ~/PathwaySplice/inst/extdata/GTF_Files/flat.chr22.",n,".genes.gff")

lapply(2:length(gene.name),function(u,gene.name,n){
  
  if(u <= n){
    
  cmd7 <- paste0("awk -v s1=","'",gene.name[u],"'")
  
  cmd8 <- paste(cmd7,cmd2,cmd3,cmd6)

  system(cmd8)
  }
  
},gene.name,n)

xx <- gene.name[1]

sample.list <- dir("~/PathwaySplice/inst/extdata/",pattern = "chr22")

cmd.list <- lapply(sample.list,function(u,xx){
  
  cat(u,"\n")
  cmd9 <- "grep"
  cmd10 <- "~/PathwaySplice/inst/extdata/"
  cmd11 <- "/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt"
  cmd12 <- ">"
  cmd13 <- paste0("/Counts.",n,".genes.txt")
  xxx <- gsub(";","",xx)
  cmd14 <- paste(cmd9,xxx,paste0(cmd10,u,cmd11),cmd12,paste0(cmd10,u,cmd13),sep = " ")
  
  system(cmd14)
  
  cmd14
  
},xx)


lapply(2:length(gene.name),function(x,gene.name,sample.list,n){
  
  
  if(x <= n){
    
  
  xx <- gene.name[x]
   
  cmd.list <- lapply(sample.list,function(y,xx){
    
    cat(y,"\n")
    cmd15 <- "grep"
    cmd16 <- "~/PathwaySplice/inst/extdata/"
    cmd17 <- "/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt"
    cmd18 <- ">>"
    cmd19 <- paste0("/Counts.",n,".genes.txt")
    xxx <- gsub(";","",xx)
    cmd20 <- paste(cmd15,xxx,paste0(cmd16,y,cmd17),cmd18,paste0(cmd16,y,cmd19),sep = " ")
    
    system(cmd20)
    
    cmd20
    
  },xx)
  
  cmd.list
  
  }
  
},gene.name,sample.list,n)

cmd21 <- paste(paste0("head -",n),"~/PathwaySplice/inst/extdata/c2.cp.v5.2.symbols.gmt.txt",sep=" ")

cmd22 <- " > ~/PathwaySplice/inst/extdata/10.cp.gmt.txt" 

cmd23 <-paste(cmd21,cmd22)

cmd23

system(cmd23)


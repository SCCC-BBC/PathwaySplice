
cmd0 <- "awk '{print $10}' ~/PathwaySplice/inst/extdata/GTF_Files/flat.chr22.gff | sort -u"

gene.name <- try(system(cmd0,intern =TRUE))

n <- 10

# "awk -v s1="ENSG00000283023;" -v s2="ENSG00000283047;" '$10==s1 || $10==s2 {print}' PathwaySplice/inst/extdata/GTF_Files/flat.chr22.gff"

cmd1 <- paste0("awk -v s1=","'",gene.name[1],"'")

cmd2 <- "'$10==s1 {print}'" 

cmd3 <- "~/PathwaySplice/inst/extdata/GTF_Files/flat.chr22.gff"

cmd4 <- paste0("> ~/PathwaySplice/inst/extdata/GTF_Files/flat.chr22.",n,".genes.gff")

cmd11 <- paste(cmd1,cmd2,cmd3,cmd4)

system(cmd11)
cmd5 <- paste0(">> ~/PathwaySplice/inst/extdata/GTF_Files/flat.chr22.",n,".genes.gff")

lapply(2:length(gene.name),function(u,gene.name,n){
  
  if(u <= n){
    
  cmd.temp <- paste0("awk -v s1=","'",gene.name[u],"'")
  
  cmd12 <- paste(cmd.temp,cmd2,cmd3,cmd5)

  system(cmd12)
  }
  
},gene.name,n)


# t2 <- try(system(cmd,intern =TRUE))
# 
# gene.name <- "ENSG00000283047"
# 
# sh_cmd1 <- "grep ENSG00000283047 PathwaySplice/inst/extdata/GTF_Files/flat.chr22.gff > PathwaySplice/inst/extdata/GTF_Files/flat.chr22.one.gene.gff"
# 
# system(sh_cmd1)

xx <- gene.name[1]

sample.list <- dir("~/PathwaySplice/inst/extdata/",pattern = "chr22")

cmd.list <- lapply(sample.list,function(u,xx){
  
  cat(u,"\n")
  cmd1 <- "grep"
  cmd2 <- "~/PathwaySplice/inst/extdata/"
  cmd3 <- "/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt"
  cmd4 <- ">"
  cmd5 <- paste0("/Counts.",n,".genes.txt")
  xxx <- gsub(";","",xx)
  cmd <- paste(cmd1,xxx,paste0(cmd2,u,cmd3),cmd4,paste0(cmd2,u,cmd5),sep = " ")
  
  system(cmd)
  
  cmd
  
},xx)


lapply(2:length(gene.name),function(x,gene.name,sample.list,n){
  
  
  if(x <= n){
    
  
  xx <- gene.name[x]
   
  cmd.list <- lapply(sample.list,function(y,xx){
    
    cat(y,"\n")
    cmd1 <- "grep"
    cmd2 <- "~/PathwaySplice/inst/extdata/"
    cmd3 <- "/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt"
    cmd4 <- ">>"
    cmd5 <- paste0("/Counts.",n,".genes.txt")
    xxx <- gsub(";","",xx)
    cmd <- paste(cmd1,xxx,paste0(cmd2,y,cmd3),cmd4,paste0(cmd2,y,cmd5),sep = " ")
    
    system(cmd)
    
    cmd
    
  },xx)
  
  cmd.list
  
  }
  
},gene.name,sample.list,n)

# sh_cmd2 <- "grep ENSG00000283047 PathwaySplice/inst/extdata/"
# 
# SRR1660308_chr22
# "/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt > PathwaySplice/inst/extdata/SRR1660308_chr22/"
# 
# "QC.spliceJunctionAndExonCounts.forJunctionSeq.one.gene.txt"
# 
# system(sh_cmd2)




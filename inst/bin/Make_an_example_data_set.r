#!/usr/bin/Rscript 

#Usage: 
 
#Rscript  Rscript ~/PathwaySplice/inst/bin/Make_an_example_data_set.r

cat("Please give the input file for sample information?\n")

input<-file('stdin', 'r')
row <- readLines(input, n=1)
#print(row)

input.dir=dirname(row)
#print(input.dir)

sample.data<-read.table(row,header=TRUE)
#print(sample.data)

sample.name<-sample.data$sample.ID

cat("Please give the gff input file ?\n")

input<-file('stdin', 'r')
gff.file <- readLines(input, n=1)

#print(gff.file)

gff.file.name=basename(gff.file)

cat("Please give the output file directory?\n")

input<-file('stdin', 'r')
output.file.dir <- readLines(input, n=1)

#print(output.file.dir)

# cat("Please define how many rows you want to sample from the input file:\n")
# 
# input<-file('stdin', 'r')
# num.row <- readLines(input, n=1)
# #print(num.row)

cmd="gshuf -n 100"

#Generate these tine files
re<-lapply(1:length(sample.name),function(u,sample.name){

   print(as.character(sample.name[u]))

  if(!dir.exists(paste0(output.file.dir,as.character(sample.name[u]),"/")))
  {
    dir.create(paste0(output.file.dir,as.character(sample.name[u]),"/"))
    #system(cmd,input1,">",output1)
  }

  input1=paste0(input.dir,"/",as.character(sample.name[u]),"/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt")
  output1=paste0(output.file.dir,"/",as.character(sample.name[u]),"/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt")

  print(input1)
  print(output1)

  system(paste(cmd,input1,">",output1,collapse = " "))

  output1

},sample.name)

if(!dir.exists(paste0(output.file.dir,"/GTF_Files/")))
{
  dir.create(paste0(output.file.dir,"/GTF_Files/"))
}

output2=paste0(output.file.dir,"/GTF_Files/",gff.file.name)

system(paste(cmd,gff.file,">",output2,collapse=" "))



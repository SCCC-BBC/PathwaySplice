
callQoRT<-function(cmd1, inputfile,outfile){
cmd2=paste0(cmd1,inputfile,outfile,collapse = " ")
print(cmd2)
system(cmd2, intern = TRUE, ignore.stderr = TRUE)
}

# callQoRT("java -Xmx4000M -jar /home/aiminyan/QoRTs/QoRTsFullExampleData/QoRTsRelease/QoRTs.jar QC --singleEnded --noGzipOutput --maxReadLength 1000 --keepMultiMapped",
#          "/media/H_driver/2015/Nimer_Cheng/164.alignments.bam /media/DATA/mus_musculus/Mus_musculus.GRCm38.83.processed.sorted.gtf",
#          "/media/H_driver/2015/Nimer_Cheng/164_junction_seq_new_gtf_7")

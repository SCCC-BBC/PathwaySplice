# Make a tiny example data set using chr4 data
```{r eval=FALSE}

dir.name="/media/H_driver/Aimin_project/GOSJ_STAR_Bam/"
file.name=dir(dir.name,recursive = TRUE,pattern="sorted.bam")
file.name.whole<-paste0(dir.name,file.name)
file.name.selected<-file.name.whole
file.name.selected.2<-as.list(file.name.selected)
names(file.name.selected.2)=sapply(strsplit(file.name.selected,split="\\/"),"[[",6)
gtf1="/media/H_driver/Aimin_project/GTF_Files/Homo_sapiens.GRCh38.84.processed.sorted.2.chr4.only.gtf"
cmd4="java -Xmx5000M -jar /home/aiminyan/QoRTs/QoRTsFullExampleData/QoRTsRelease/QoRTs.jar QC --noGzipOutput --dropChrom chr1,chr10,chr11,chr11_KI270721v1_random,chr12,chr13,chr14,chr14_GL000009v2_random,chr14_GL000194v1_random,chr14_GL000225v1_random,chr14_KI270722v1_random,chr14_KI270723v1_random,chr14_KI270726v1_random,chr15,chr15_KI270727v1_random,chr16,chr16_KI270728v1_random,chr17,chr17_GL000205v2_random,chr17_KI270729v1_random,chr18,chr19,chr1_KI270706v1_random,chr1_KI270707v1_random,chr1_KI270709v1_random,chr1_KI270710v1_random,chr1_KI270711v1_random,chr1_KI270712v1_random,chr1_KI270713v1_random,chr1_KI270714v1_random,chr2,chr20,chr21,chr22,chr22_KI270731v1_random,chr22_KI270732v1_random,chr22_KI270733v1_random,chr22_KI270734v1_random,chr22_KI270735v1_random,chr22_KI270736v1_random,chr22_KI270738v1_random,chr2_KI270716v1_random,chr3,chr3_GL000221v1_random,chr4_GL000008v2_random,chr5,chr6,chr7,chr8,chr9,chr9_KI270718v1_random,chr9_KI270719v1_random,chr9_KI270720v1_random,chrM,chrUn_GL000195v1,chrUn_GL000213v1,chrUn_GL000214v1,chrUn_GL000216v2,chrUn_GL000218v1,chrUn_GL000219v1,chrUn_GL000220v1,chrUn_GL000224v1,chrUn_KI270311v1,chrUn_KI270315v1,chrUn_KI270330v1,chrUn_KI270337v1,chrUn_KI270362v1,chrUn_KI270435v1,chrUn_KI270438v1,chrUn_KI270442v1,chrUn_KI270467v1,chrUn_KI270511v1,chrUn_KI270519v1,chrUn_KI270522v1,chrUn_KI270590v1,chrUn_KI270741v1,chrUn_KI270742v1,chrUn_KI270743v1,chrUn_KI270744v1,chrUn_KI270745v1,chrUn_KI270746v1,chrUn_KI270747v1,chrUn_KI270748v1,chrUn_KI270750v1,chrUn_KI270751v1,chrUn_KI270754v1,chrX,chrY"

re.out<-lapply(file.name.selected.2,callQoRT,output.dir="chr4_drop_other",gtf_file=gtf1,runing_cmd=cmd4)

add the count data and "Homo_sapiens.GRCh38.84.processed.sorted.4.JunctionSeq.flat.chr4.100.gff" into PathwaySplice/inst/extdata/

```
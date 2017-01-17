#!/usr/bin/Rscript 

dim(counts(Re.example))
counts.12<-counts(Re.example)
colnames(counts.12)
install_github("aiminy/Sophia")

Data=counts.12

heatmap_wPCA(Data)

heatmap_wPCA(Data,g_level=c(rep("Mu",8),rep("WT",4)))


PCGextract<-function(input,output,lengthfile){
norCounts<-read.table(input,header=T,sep="\t",row.names=1)
PCGlength<-read.table(lengthfile,sep="\t",row.names=1)
intersection<-intersect(rownames(norCounts),rownames(PCGlength))
difference<-setdiff(rownames(PCGlength),intersection)
inter.data<-norCounts[intersection,]
inter.data$length<-PCGlength[intersection,]
inter.set<-inter.data[,-ncol(inter.data)]
diff.set<-as.data.frame(matrix(0,length(difference),ncol(norCounts)))
colnames(diff.set)<-colnames(norCounts)
rownames(diff.set)<-difference
extract.data<-rbind(inter.set,diff.set)
write.table(extract.data,paste(output,"_PCGs.txt",sep=""),quote=F,sep="\t")
}

PCGextract("siPHLDA1_normalize.txt","siPHLDA1","PCGsLengthMatrix.txt")
PCGextract("siPHLDA1_diff.txt","siPHLDA1_diff","PCGsLengthMatrix.txt")
PCGextract("TCGA_CRC_diff.txt","TCGA_CRC_diff","PCGsLengthMatrix.txt")
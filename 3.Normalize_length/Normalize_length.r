norTPM<-function(input,output,lengthfile){
norCounts<-read.table(input,header=T,sep="\t",row.names=1)
PCGlength<-read.table(lengthfile,sep="\t",row.names=1)
intersection<-intersect(rownames(norCounts),rownames(PCGlength))
difference<-setdiff(rownames(PCGlength),intersection)
inter.data<-norCounts[intersection,]
inter.data$length<-PCGlength[intersection,]
inter.set<-inter.data[,-ncol(inter.data)]*1000/inter.data$length
diff.set<-as.data.frame(matrix(0,length(difference),ncol(norCounts)))
colnames(diff.set)<-colnames(norCounts)
rownames(diff.set)<-difference
norCPK.data<-rbind(inter.set,diff.set)
norTPM.data<-norCPK.data*1000000/colSums(norCPK.data)
write.table(norTPM.data,paste(output,"_norTPM.txt",sep=""),quote=F,sep="\t")
}

norTPM("siPHLDA1_normalize.txt","siPHLDA1","PCGsLengthMatrix.txt")
norTPM("TCGA_CRC_normalize.txt","TCGA_CRC","PCGsLengthMatrix.txt")

Expcor<-function(input,output,targetgene){
matrix<-read.table(input,header=T,sep="\t",row.names=1)
tarexp<-as.numeric(matrix[targetgene,])
correlation<-as.data.frame(t(apply(matrix,1,function(x){
cor<-cor.test(tarexp,as.numeric(x),method="spearman")
c(as.numeric(cor[4]),as.numeric(cor[3]))
})))
colnames(correlation)<-c("rho","pvalue")
rownames(correlation)<-rownames(matrix)
write.table(correlation,paste(output,"_correlation.txt",sep=""),quote=F,sep="\t")
}

Expcor("TCGA_CRC_norTPM.txt","TCGA_CRC","PHLDA1")
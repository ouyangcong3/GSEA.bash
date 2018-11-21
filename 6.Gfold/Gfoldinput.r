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

PCGextract("siPHLDA1_counts.txt","siPHLDA1_rawcount","PCGsLengthMatrix.txt")
PCGextract("BGI_siPHLDA1_FPKM.txt","BGI_siPHLDA1_FPKM","PCGsLengthMatrix.txt")

################################################################################################################################################
Gfoldinput<-function(expfile,pcgfile,rpkmfile){
exp<-read.table(expfile,header=T,sep="\t",row.names=1)
pcg<-read.table(pcgfile,sep="\t",row.names=1)
rpkm<-read.table(rpkmfile,header=T,sep="\t",row.names=1)
GeneSymbol<-rownames(pcg)
GeneNames<-rep("NA",nrow(pcg))
GeneExonLength<-pcg[,1]
for(i in 1:ncol(exp)){
ReadCount<-as.numeric(exp[GeneSymbol,i])
RPKM<-as.numeric(rpkm[GeneSymbol,i])
Matrix<-data.frame(GeneSymbol,GeneNames,ReadCount,GeneExonLength,RPKM)
write.table(Matrix,paste(colnames(exp)[i],".cnt",sep=""),sep="\t",quote=F,row.names=F,col.names=F)
}}

Gfoldinput("siPHLDA1_rawcount_PCGs.txt","PCGsLengthMatrix.txt","BGI_siPHLDA1_FPKM_PCGs.txt")
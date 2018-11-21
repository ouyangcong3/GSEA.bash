DESeq2diff<-function(Data,Outputname,Connum,Treatnum){
require(DESeq2)
require(limma)
rt=read.table(Data,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,-1]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
countData<-round(avereps(data),0)
condition <- factor(c(rep("con",Connum),rep("treat",Treatnum)))
dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), design= ~ condition )
dds <- DESeq(dds)
res <- results(dds)
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)
write.table(resdata[,1:7],file = paste(Outputname,"_diff.txt",sep=""),sep="\t",quote=F,row.names=F)
write.table(resdata[,-c(2:7)],file = paste(Outputname,"_normalize.txt",sep=""),sep="\t",quote=F,row.names=F)
}

DESeq2diff("_counts.txt","",2,2)
DESeq2diff("mRNA.symbol.txt","TCGA_CRC",51,647)
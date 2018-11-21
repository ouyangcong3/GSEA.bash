Averreplength<-function(input,output){
require(limma)
genelength<-read.table(input,sep="\t",header=T)
genelength=as.matrix(genelength)
rownames(genelength)=genelength[,1]
lengths=genelength[,-1]
dimnames=list(rownames(genelength),colnames(lengths))
data=matrix(as.numeric(as.matrix(lengths)),nrow=length(lengths),dimnames=dimnames)
genelengthmatrix<-round(avereps(data),0)
write.table(genelengthmatrix,paste(output,".txt",sep=""),sep="\t",quote=F,col.names=F)
}

Averreplength("Ensembl_gene_length.txt","GeneLengthMatrix")
Averreplength("Ensembl_PCGs_length.txt","PCGsLengthMatrix")
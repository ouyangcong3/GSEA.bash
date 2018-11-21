GSVAanalysis<-function(input,output,dataset,Connum,Treatnum){
require(GSVA)
require(pheatmap)
require(Biobase)
require(genefilter)
require(limma)
require(RColorBrewer)
require(GSEABase)
require("impute")
#GSVA
myC7 <- getGmt(dataset, geneIdType=EntrezIdentifier(), collectionType=BroadCollection(category="c7"), sep="\t")
exp_matrix<-as.matrix(read.table(input,header=T,row.names=1))
es <- gsva(exp_matrix, myC7, mx.diff=FALSE, verbose=FALSE, parallel.sz=1)
write.table(es[[1]],paste(output,"_GSVAresult.txt",sep=""),quote=F,sep="\t")
ESmean<-as.data.frame(rowMeans(abs(es[[1]])))

PvalueCutoff <- 0.05
logFCcutoff <- log2(2)

exp<-es[[1]]
dimnames=list(rownames(exp),colnames(exp))
exp=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
#impute missing expression data
mat=impute.knn(exp)
rt=mat$data
rt=avereps(rt)
#normalize
#boxplot(rt,col = "blue",xaxt = "n",outline = F)
#rt=normalizeBetweenArrays(as.matrix(rt))
#boxplot(rt,col = "red",xaxt = "n",outline = F)
#differential
class <- c(rep("con",Connum),rep("treat",Treatnum))
design <- model.matrix(~0+factor(class))
colnames(design) <- c("con","treat")
fit <- lmFit(rt,design)
cont.matrix<-makeContrasts(treat-con,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
allDiff=topTable(fit2,adjust='fdr',number=200000)
allDiff$AveExpr<-ESmean[rownames(allDiff),]
write.table(allDiff,paste(output,"_GSVA_limmaTab.txt",sep=""),sep="\t",quote=F)
#write table
diffSig <- allDiff[with(allDiff, (abs(logFC)>logFCcutoff & P.Value < PvalueCutoff )), ]
hmExp=rt[rownames(diffSig),]
diffExp=rbind(id=colnames(hmExp),hmExp)
rownames(diffSig)<-tolower(x=rownames(diffSig))
write.table(diffSig,paste(output,"_GSVA_diff.txt",sep=""),sep="\t",quote=F)
rownames(diffExp)<-tolower(x=rownames(diffExp))
write.table(diffExp,paste(output,"_GSVA_diffExp.txt",sep=""),sep="\t",quote=F,col.names=F)
}

#GSVAanalysis("siPHLDA1_norTPM.txt","siPHLDA1","h+c2all+c3tft+c4all+c5all+c6+c7.symbols.gmt",2,2)
GSVAanalysis("TCGA_CRC_norTPM.txt","TCGA_CRC","h+c2all+c3tft+c4all+c5all+c6+c7.symbols.gmt",51,647)
#install.packages("FactoMineR")
library(DESeq2)
library(plotly)
library("FactoMineR")
library("rmarkdown")
library("vegan")
library("manhattanly")
library("DT")
library("heatmaply")
library(gplots)
#library("FactoMineR",lib.loc="/home/wud3/R/x86_64-pc-linux-gnu-library/3.4/")
#library("rmarkdown",lib.loc="/home/wud3/R/x86_64-pc-linux-gnu-library/3.4/")
#library("vegan",lib.loc="/home/wud3/R/x86_64-pc-linux-gnu-library/3.4/")
#library("manhattanly",lib.loc="/home/wud3/R/x86_64-pc-linux-gnu-library/3.4/")
#library("DT",lib.loc="/home/wud3/R/x86_64-pc-linux-gnu-library/3.4/")
#library("heatmaply",lib.loc="/home/wud3/R/x86_64-pc-linux-gnu-library/3.4/")
rm(list=ls())
TIER2 <- Sys.getenv("TIER2")
#a1=read.table("JL-6152--01--14--2019_COUNTS_new.csv",sep=',',header=T,row.names=1,check=F,comment.char="")
args=commandArgs(TRUE)
count_file <- args[1]
sample_info <-args[2]
comparison <- args[3]
project <- args[4]
species <- args[5]
a1=read.table(count_file, sep=',',header=T,row.names=1,check=F,comment.char="")
a2=read.table(sample_info,sep=',',header=T)
a3=read.table(comparison,sep=',',header=F)
comps <- as.matrix(a3) 

# sort a2 by group and then by sample name
a2 <- a2[order( a2[,3], a2[,2] ),]

a1<-a1[,as.vector(a2[,1])]

colnames(a1) <- as.matrix(a2)[,2]
#data<-a1[-grep("ERCC-",rownames(a1)),]
countdataraw=round(a1[rowSums(a1)>0,])
countdata<-as.matrix(countdataraw)
dim(countdata)
temp <- data.frame(countdata)
temp$Gene <- row.names(temp)
#library("AnnotationDbi",lib.loc="/home/genomics/anaconda2/lib/R/library/") #add gene name to the file based on ENSEMBL
library("AnnotationDbi")
if (species == "human"){
	library("org.Hs.eg.db")
	temp$symbol = mapIds(org.Hs.eg.db,
                          keys=temp$Gene, 
                          column="SYMBOL",
                          keytype="ENSEMBL",
                          multiVals="first")
}else if (species == "mouse"){
	library("org.Mm.eg.db")
	temp$symbol = mapIds(org.Mm.eg.db,
			     keys=temp$Gene,
			     column="SYMBOL",
			     keytype="ENSEMBL",
			     multiVals="first")

}else{
	print('species name needed; input mouse or human')  
}



row.names(temp) <- paste(temp$Gene,temp$symbol, sep = "_")
n <- dim(a2)[1]
countdata <- as.matrix(temp[,1:n])
condition<-factor(as.matrix(a2)[,3])
coldata <- data.frame(row.names=colnames(countdata), condition)
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds=DESeq(dds,minReplicatesForReplace = 50)
norm=counts(dds,normalized = T)

#### get PCA output for all samples as PDF files ####
rld=rlog(dds, blind=TRUE)
rv <- rowVars(assay(rld)) # caluculate row variance
select <- order(rv, decreasing = TRUE)[seq_len(min(500,length(rv)))]
dat.norm<-t(assay(rld)[select, ])
dat.norm<-data.frame(dat.norm,condition)
res.pca=PCA(dat.norm,ncp=5,scale.unit=T,graph=F,quali.sup=ncol(dat.norm))

pdf(paste(project,"_PCA_all_samples.pdf",sep=""),16,12)
plot.PCA(res.pca,axes=c(1,2),habillage=ncol(dat.norm),cex=1)
#dev.off()
plot.PCA(res.pca,axes=c(1,2),habillage=ncol(dat.norm),cex=2,label='none')
#pdf(paste(project,"_PCA1vs3.pdf",sep=""),16,12)
plot.PCA(res.pca,axes=c(1,3),habillage=ncol(dat.norm),cex=1)
#dev.off()
plot.PCA(res.pca,axes=c(1,3),habillage=ncol(dat.norm),cex=2,label='none')
#pdf(paste(project,"_PCA2vs3.pdf",sep=""),16,12)
plot.PCA(res.pca,axes=c(2,3),habillage=ncol(dat.norm),cex=1)
#dev.off()
plot.PCA(res.pca,axes=c(2,3),habillage=ncol(dat.norm),cex=2,label='none')
dev.off()
#dev.off()

out_folder <- paste0(getwd(),"/")
##################################### get DEGs list, interative plot, and  pathway plot for each comparison ##############################################
# read how many comparisons are needed and do analysis for each comparison and save results as res1, res2, res3...

for (i in 1:dim(comps)[1]) {
  name <- paste0(comps[i,][2], "_vs_", comps[i,][3])  #"M_vs_MC"
  pattern <- c(comps[i,][2], comps[i,][3])
  countdata_1 <- countdata[, grep(paste(pattern, collapse="|"), colnames(countdata), value = TRUE)]
  condition<-factor(as.matrix(a2)[,3])
  condition <- factor(grep(paste(pattern, collapse="|"), condition, value = TRUE))
  coldata_1 <- data.frame(row.names=colnames(countdata_1), condition)
  dds_1 <- DESeqDataSetFromMatrix(countData=countdata_1, colData=coldata_1, design=~condition)
  dds_1=DESeq(dds_1,minReplicatesForReplace = 50)
  norm_1=counts(dds_1,normalized = T)
  rld_1=rlog(dds_1, blind=TRUE)
  rv_1 <- rowVars(assay(rld_1)) # caluculate row variance
  select_1 <- order(rv_1, decreasing = TRUE)[seq_len(min(500,length(rv_1)))]
  dat.norm_1<-t(assay(rld_1)[select_1, ])
  dat.norm_1<-data.frame(dat.norm_1,condition)
  res.pca_1=PCA(dat.norm_1,ncp=5,scale.unit=T,graph=F,quali.sup=ncol(dat.norm_1))
  pdf(paste(name,"_PCA_top500.pdf",sep=""),16,12)
  plot.PCA(res.pca_1,axes=c(1,2),habillage=ncol(dat.norm),cex=1)
  plot.PCA(res.pca_1,axes=c(1,2),habillage=ncol(dat.norm),cex=2,label='none')
  plot.PCA(res.pca_1,axes=c(1,3),habillage=ncol(dat.norm),cex=1)
  plot.PCA(res.pca_1,axes=c(1,3),habillage=ncol(dat.norm),cex=2,label='none')
  plot.PCA(res.pca_1,axes=c(2,3),habillage=ncol(dat.norm),cex=1)
  plot.PCA(res.pca_1,axes=c(2,3),habillage=ncol(dat.norm),cex=2,label='none')
  dev.off()

  #norm_1 <- norm[, grep(paste(pattern,collapse="|"), colnames(norm), value = TRUE)]
  assign(paste("res", i, sep = ""), results(dds_1, contrast = comps[i,]))
  assign(paste0("resdata",i), merge(as.data.frame(eval(parse(text=paste0("res", i)))), as.data.frame(norm_1), by="row.names", sort=FALSE))
  assign(paste("resdata", i, sep = ""), eval(parse(text=paste0("resdata", i)))[,-c(4,5)]) #remove lfcSE & stat columns
  # eval(parse(text=paste0("resdata", i))) get value from variable paste0("resdata", i), resdata1, resdata2, resdata3...
  write.csv(eval(parse(text=paste0("resdata", i))), file=paste(name,"_DEGs_all.csv",sep=''),row.names=F)
  n1 <- dim(subset(eval(parse(text=paste0("resdata", i))), padj<0.05))[1]
  n2 <- dim(subset(eval(parse(text=paste0("resdata", i))), padj<0.01))[1]
  if (n1 < 2000 & n1 >= 200) {
    resSig=subset(eval(parse(text=paste0("resdata", i))), padj<0.05)
    write.csv(resSig,file=paste(name,"_DEGs_padj0.05.csv",sep=''),row.names=F)
    sig_DEGs <- "p_adj < 0.05"
  }
  if (n1 >= 2000) {
    if (n2 < 2000) {
      resSig=subset(eval(parse(text=paste0("resdata", i))), padj<0.01)
      write.csv(resSig,file=paste(name,"_DEGs_padj0.01.csv",sep=''),row.names=F)
      sig_DEGs <- "p_adj < 0.01"
    } else {
      resSig=subset(eval(parse(text=paste0("resdata", i))), padj<0.01 & abs(log2FoldChange)>1)
      write.csv(resSig,file=paste(name,"_DEGs_padj0.01_FC2.csv",sep=''),row.names=F)
      sig_DEGs <- "p_adj < 0.01 and |FC| > 2"
    }
  }
  if (n1 < 200 & n1 >= 20) {
    resSig=subset(eval(parse(text=paste0("resdata", i))), padj<0.1)
    write.csv(resSig,file=paste(name,"_DEGs_padj0.1.csv",sep=''),row.names=F)
    sig_DEGs <- "p_adj < 0.1"
  }
  if (n1 < 20) {
    resSig=subset(eval(parse(text=paste0("resdata", i))), pvalue<0.1)
    write.csv(resSig,file=paste(name,"_DEGs_p0.05.csv",sep=''), row.names=F)
    sig_DEGs <- "p < 0.05"
  }
  rmarkdown::render(paste0(TIER2,"/Automatic/Interactive_report_RNAseq.Rmd"), params = list(data = a1, info = a2, comparison = a3, project = project), output_file = paste0(name, ".html"), output_dir="./")
  i=i+1 
};rm(i)


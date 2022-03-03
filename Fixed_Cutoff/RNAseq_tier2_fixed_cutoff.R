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
#a1=read.table("AP-5782--11--08--2018_COUNTS.csv",sep=',',header=T,row.names=1,check=F,comment.char="")
args=commandArgs(TRUE)
count_file <- args[1]
sample_info <-args[2]
comparison <- args[3]
project <- args[4]
sig_DEGs_padj <- args[5]
sig_DEGs_FC <- args[6]

a1=read.table(count_file, sep=',',header=T,row.names=1,check=F,comment.char="")
a2=read.table(sample_info,sep=',',header=T)
a3=read.table(comparison,sep=',',header=F)
#a2=read.table("AP-5782--11--08--2018_sample_info.csv",sep=',',header=T)
#a3=read.table("AP-5782--11--08--2018_comparisons.csv",sep=',',header=F)
#project = "AP-5782--11--08--2018"
comps <- as.matrix(a3)
a2 <- a2[order( a2[,3], a2[,2] ),] # sort a2 by group and then by sample name
a1<-a1[,as.vector(a2[,1])] #reorder a1

colnames(a1) <- as.matrix(a2)[,2]
data<-a1[-grep("ERCC-",rownames(a1)),]
countdataraw=round(data[rowSums(data)>0,])
countdata<-as.matrix(countdataraw)
dim(countdata)


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

pdf(paste("./",project,"_PCA_all_samples.pdf",sep=""),16,12)
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

#out_folder <- paste0(getwd(),"/")
##################################### get DEGs list, interative plot, and  pathway plot for each comparison ##############################################
# read how many comparisons are needed and do analysis for each comparison and save results as res1, res2, res3...

for (i in 1:dim(comps)[1]) {
  name <- paste0(comps[i,][2], "_vs_", comps[i,][3])  #"M_vs_MC"
  cutoff1_temp <- gsub(" ", "", sig_DEGs_padj, fixed = TRUE) #delete space from string
  temp <- unlist(strsplit(cutoff1_temp, "<")) #temp[1] is either pvalue or padj; temp[2] is the cutoff
  folder_path <- paste0(name, "_", temp[1], temp[2], "_FC_", sig_DEGs_FC)
  dir.create(folder_path)
  out_folder <- paste0(getwd(),"/",folder_path,"/")

  #dir.create(name)
  #out_folder <- paste0(getwd(),"/",name,"/")
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
  pdf(paste(out_folder,name,"_PCA_top500.pdf",sep=""),16,12)
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
  temp_resdata<- eval(parse(text=paste0("resdata", i)))
  colnames(temp_resdata)[5]<-"p_adj"
  write.csv(temp_resdata, file=paste(out_folder,name,"_DEGs_all.csv",sep=''),row.names=F)

#sig_DEGs_padj <- "pvalue<0.05"
#sig_DEGs_padj <- "padj<0.05"

  cutoff1_temp <- gsub(" ", "", sig_DEGs_padj, fixed = TRUE) #delete space from string
  temp <- unlist(strsplit(cutoff1_temp, "<")) #temp[1] is either pvalue or padj; temp[2] is the cutoff value

  resSig=subset(temp_resdata, temp_resdata[,temp[1]]<as.numeric(temp[2]) & abs(log2FoldChange)>=log2(as.numeric(sig_DEGs_FC)))
  if (sig_DEGs_FC==1) {
    sig_DEGs <- sig_DEGs_padj
    write.csv(resSig,file=paste(out_folder, name,"_DEGs_",sig_DEGs_padj,".csv",sep=''),row.names=F)

  }
  else if (sig_DEGs_FC>1) {
    sig_DEGs <- paste0(sig_DEGs_padj, " and ", "|FC| > ", sig_DEGs_FC)
    write.csv(resSig,file=paste(out_folder, name,"_DEGs_",sig_DEGs_padj,"_FC",sig_DEGs_FC,".csv",sep=''),row.names=F)
  }

  rmarkdown::render(paste0(TIER2,"/Fixed_Cutoff/Interactive_report_RNAseq_fixed_cutoff.Rmd"), params = list(data = a1, info = a2, comparison = a3, project = project, cutoff1 = sig_DEGs_padj, cutoff2 = sig_DEGs_FC), output_file = paste0(name, ".html"),output_dir=out_folder)
  i=i+1
};rm(i)

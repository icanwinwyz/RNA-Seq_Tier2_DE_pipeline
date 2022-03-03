library(ggplot2)

args=commandArgs(TRUE)
kegg<-args[1]
kegg<-read.delim(kegg,check.names=F,header=T,row.names=NULL,quote="",stringsAsFactors=FALSE)

data<-data.frame()
for (i in 2:4){
	data<-rbind(data,read.delim(args[i],check.names=F,header=T,row.names=NULL,quote="",stringsAsFactors=FALSE))
}
#test<-kegg[,c(2,3,10,12)]
test<-kegg[,c(2,3,10,5)]
#test<-subset(test,PValue<0.1)
#colnames(test)<-c("Term","Gene_Number","Fold_Enrichment_Score","adjusted_pvalue")
colnames(test)<-c("Term","Gene_Number","Fold_Enrichment_Score","pvalue")
name<-args[5]

pdf(paste(name,"DE_gene_KEGG.pdf",sep = "_"),8,10)
#ggplot(test,aes(x=Fold_Enrichment_Score,y=Term,size=Gene_Number,colour=adjusted_pvalue))+geom_point()+scale_colour_gradient(low="red",high="blue")+xlab("Fold Enrichment Score")+ylab("Term")+ggtitle("Enriched KEGG pathway")
ggplot(test,aes(x=Fold_Enrichment_Score,y=Term,size=Gene_Number,colour=pvalue))+geom_point()+scale_colour_gradient(low="red",high="blue")+xlab("Fold Enrichment Score")+ylab("Term")+ggtitle("Enriched KEGG pathway")
dev.off()

#bmp(paste(name,"DE_gene_KEGG.bmp",sep = "_"))
#ggplot(test,aes(x=Fold_Enrichment_Score,y=Term,size=Gene_Number,colour=adjusted_pvalue))+geom_point()+scale_colour_gradient(low="red",high="blue")+xlab("Fold Enrichment Score")+ylab("Term")+ggtitle("Enriched KEGG pathway")
#dev.off()


#jpeg(paste(name,"DE_gene_KEGG.jpeg",sep = "_"),quality = 100)
#ggplot(test,aes(x=Fold_Enrichment_Score,y=Term,size=Gene_Number,colour=adjusted_pvalue))+geom_point()+scale_colour_gradient(low="red",high="blue")+xlab("Fold Enrichment Score")+ylab("Term")+ggtitle("Enriched KEGG pathway")
#dev.off()


#png(paste(name,"DE_gene_KEGG.png",sep = "_"))
#ggplot(test,aes(x=Fold_Enrichment_Score,y=Term,size=Gene_Number,colour=adjusted_pvalue))+geom_point()+scale_colour_gradient(low="red",high="blue")+xlab("Fold Enrichment Score")+ylab("Term")+ggtitle("Enriched KEGG pathway")
#dev.off()


#tiff(paste(name,"DE_gene_KEGG.tiff",sep = "_"),width=9,height=8,units='in',res=300)
#ggplot(test,aes(x=Fold_Enrichment_Score,y=Term,size=Gene_Number,colour=adjusted_pvalue))+geom_point()+scale_colour_gradient(low="red",high="blue")+xlab("Fold Enrichment Score")+ylab("KEGG Term")+theme(axis.text.y=element_text(size=14),axis.title.y=element_text(size=14),axis.title.x=element_text(size=14))+ggtitle("Enriched KEGG pathway")
#ggplot(test,aes(x=Fold_Enrichment_Score,y=Term,size=Gene_Number,colour=adjusted_pvalue))+geom_point()+scale_colour_gradient(low="red",high="blue")+xlab("Fold Enrichment Score")+ylab("KEGG Term")+theme(axis.text.y=element_text(size=14),axis.title.y=element_text(size=14),axis.title.x=element_text(size=14))
#dev.off()


#svg(paste(name,"DE_gene_KEGG.svg",sep = "_"),width=9,height=8)
#ggplot(test,aes(x=Fold_Enrichment_Score,y=Term,size=Gene_Number,colour=adjusted_pvalue))+geom_point()+scale_colour_gradient(low="red",high="blue")+xlab("Fold Enrichment Score")+ylab("Term")+ggtitle("Enriched KEGG pathway")
#dev.off()

#data<-subset(data,Benjamini<0.1)
data<-subset(data,PValue<0.1)
data[,1]<-sub("GOTERM_","",data[,1])
data[,1]<-sub("_DIRECT","",data[,1])
#data<-data[,c(2,3,10,12,1)]
data<-data[,c(2,3,10,5,1)]
#colnames(data)<-c("Term","Gene_Number","Fold_Enrichment_Score","adjusted_pvalue","Category")
colnames(data)<-c("Term","Gene_Number","Fold_Enrichment_Score","pvalue","Category")
data$Category<-gsub("GOTERM_BP_FAT","Biological Process",data$Category)
data$Category<-gsub("GOTERM_CC_FAT","Cellular Component",data$Category)
data$Category<-gsub("GOTERM_MF_FAT","Molecular Function",data$Category)
nameorder<-data$Term[order(data$Category)]
data$Term<-factor(data$Term,levels=nameorder)
pdf(paste(name,"DE_gene_GO_term.pdf",sep="_"),15,35)
#ggplot(data,aes(x= Fold_Enrichment_Score,y=Term,size=Gene_Number,colour= adjusted_pvalue))+geom_point()+scale_colour_gradient(low="red",high="blue")+xlab("Fold Enrichment Score")+ylab("GO Term")+ggtitle("Enriched GO Term")+facet_grid(Category~.,scales="free_y",space="free_y")
ggplot(data,aes(x= Fold_Enrichment_Score,y=Term,size=Gene_Number,colour=pvalue))+geom_point()+scale_colour_gradient(low="red",high="blue")+xlab("Fold Enrichment Score")+ylab("GO Term")+ggtitle("Enriched GO Term")+facet_grid(Category~.,scales="free_y",space="free_y")
dev.off()

#tiff(paste(name,"DE_gene_GO_term.tiff",sep = "_"),width=30,height=50,units='in',res=250)
#ggplot(test,aes(x=Fold_Enrichment_Score,y=Term,size=Gene_Number,colour=adjusted_pvalue))+geom_point()+scale_colour_gradient(low="red",high="blue")+xlab("Fold Enrichment Score")+ylab("KEGG Term")+theme(axis.text.y=element_text(size=14),axis.title.y=element_text(size=14),axis.title.x=element_text(size=14))+ggtitle("Enriched KEGG pathway")
#ggplot(test,aes(x=Fold_Enrichment_Score,y=Term,size=Gene_Number,colour=adjusted_pvalue))+geom_point()+scale_colour_gradient(low="red",high="blue")+xlab("Fold Enrichment Score")+ylab("KEGG Term")+theme(axis.text.y=element_text(size=14),axis.title.y=element_text(size=14),axis.title.x=element_text(size=14))
#dev.off()
#ggplot(data,aes(x= Fold_Enrichment_Score,y=Term,size=Gene_Number,colour= adjusted_pvalue))+geom_point()+scale_colour_gradient(low="red",high="blue")+xlab("Fold Enrichment Score")+ylab("GO Term")+facet_grid(Category~.,scales="free_y",space="free_y")+scale_size_area(max_size=20)+theme(axis.text.y=element_text(size=35,face="bold",colour="black"),axis.title.y=element_text(size=40,face="bold"),axis.title.x=element_text(size=40,face="bold"),legend.title=element_text(size=25),legend.text=element_text(size=25),strip.text=element_text(face="bold",size=30))
#dev.off()

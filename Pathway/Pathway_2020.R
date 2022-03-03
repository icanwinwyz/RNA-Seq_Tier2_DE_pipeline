library(ggplot2)

args=commandArgs(TRUE)
kegg<-args[1]
kegg<-read.delim(kegg,check.names=F,header=T,row.names=NULL,quote="",stringsAsFactors=FALSE)

data<-data.frame()
for (i in 2:4){
  data<-rbind(data,read.delim(args[i],check.names=F,header=T,row.names=NULL,quote="",stringsAsFactors=FALSE))
}
#test<-kegg[,c(2,3,10,12)]
test<-kegg[,c(2,3,6,10,5,12)]
#colnames(test)[4] <- "adj-p"
#test<-subset(test,PValue<0.1)
#colnames(test)<-c("Term","Gene_Number","Fold_Enrichment_Score","adjusted_pvalue")
colnames(test)<-c("Term","Gene_Number","Gene","Fold_Enrichment_Score","pvalue","adj-p")
name<-args[5]


pdf(paste(name,"DE_gene_KEGG.pdf",sep = "_"),8,10) #5
#ggplot(test,aes(x=Fold_Enrichment_Score,y=Term,size=Gene_Number,colour=adjusted_pvalue))+geom_point()+scale_colour_gradient(low="red",high="blue")+xlab("Fold Enrichment Score")+ylab("Term")+ggtitle("Enriched KEGG pathway")
ggplot(test,aes(x=Fold_Enrichment_Score,y=Term,size=Gene_Number,colour=pvalue))+geom_point()+scale_colour_gradient(low="red",high="blue")+xlab("Fold Enrichment Score")+ylab("Term")+ggtitle("Enriched KEGG pathway")
dev.off()


#data<-subset(data,Benjamini<0.1)
data<-subset(data,PValue<0.05)
data[,1]<-sub("GOTERM_","",data[,1])
data[,1]<-sub("_DIRECT","",data[,1])
#data<-data[,c(2,3,10,12,1)]
data<-data[,c(2,3,6,10,5,12,1)]
#colnames(data)[5] <- "adj-p"
#colnames(data)<-c("Term","Gene_Number","Fold_Enrichment_Score","adjusted_pvalue","Category")
colnames(data)<-c("Term","Gene_Number","Genes","Fold_Enrichment_Score","pvalue","adj-p","Category")
data$Category<-gsub("GOTERM_BP_FAT","Biological Process",data$Category)
data$Category<-gsub("GOTERM_CC_FAT","Cellular Component",data$Category)
data$Category<-gsub("GOTERM_MF_FAT","Molecular Function",data$Category)
nameorder<-data$Term[order(data$Category)]
data$Term<-factor(data$Term,levels=nameorder)
pdf(paste(name,"DE_gene_GO_term.pdf",sep="_"),15,25) #25
#ggplot(data,aes(x= Fold_Enrichment_Score,y=Term,size=Gene_Number,colour= adjusted_pvalue))+geom_point()+scale_colour_gradient(low="red",high="blue")+xlab("Fold Enrichment Score")+ylab("GO Term")+ggtitle("Enriched GO Term")+facet_grid(Category~.,scales="free_y",space="free_y")
ggplot(data,aes(x= Fold_Enrichment_Score,y=Term,size=Gene_Number,colour=pvalue))+geom_point()+scale_colour_gradient(low="red",high="blue")+xlab("Fold Enrichment Score")+ylab("GO Term")+ggtitle("Enriched GO Term")+facet_grid(Category~.,scales="free_y",space="free_y")
dev.off()

#library(xlsx)
#write.xlsx2(test,paste0(name,"_DEGs_pathway_enrichment_.xlsx"), sheetName="KEGG", row.names=FALSE, append=TRUE)
#write.xlsx2(data,paste0(name,"_DEGs_pathway_enrichment_.xlsx"), sheetName="GO_term", row.names=FALSE, append=TRUE)

write.csv(test,paste0(name,"_DEGs_KEGG_enrichment.csv"))
write.csv(data,paste0(name,"_DEGs_GO_term_enrichment.csv"))

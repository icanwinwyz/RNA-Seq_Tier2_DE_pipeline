library(ggplot2)
library(cowplot)
args=commandArgs(TRUE)

# Do KEGG pathway ---------------------------------------------------------
# Read in data
kegg<-read.delim(args[1],check.names=F,header=T,row.names=NULL,quote="",stringsAsFactors=FALSE)
# Subset and sort by p-value
kegg<-kegg[,c(2,3,6,10,5,12)]
colnames(kegg)<-c("Term","Gene_Number","Gene","Fold_Enrichment_Score","pvalue","adj-p")
kegg$Term <- factor(kegg$Term, 
                    levels=kegg$Term[order(-kegg$pvalue)])
# Write to CSV
write.csv(kegg,paste0(args[5],"_DEG_KEGG_enrichment.csv"))

# Make the plot
ggplot(kegg[kegg$pvalue<0.05,],
       aes(x=Fold_Enrichment_Score,
           y=Term,
           size=Gene_Number,
           colour=pvalue))+
  geom_point()+
  scale_colour_gradient(low="red",
                        high="blue",
                        name="p-value")+
  scale_size_continuous(name="Num. Genes") +
  xlab("Fold Enrichment Score")+
  ylab("Term")+
  theme_minimal_grid() +
  ggtitle("Enriched KEGG pathway")

# Save with automatic sizing for height
ggsave2(filename = paste0(args[5], "_DEG_KEGG.pdf"),
        height=dim(kegg)[1]*0.25+2.5,
        width=10)

# Do GO Terms -------------------------------------------------------------
# Read in data
data<-data.frame()
for (i in 2:4){
  data<-rbind(data,read.delim(args[i],
                              check.names=F,
                              header=T,
                              row.names=NULL,
                              quote="",
                              stringsAsFactors=FALSE))
}

# Clean category names
data[grepl("BP", data[,1]),1] <- "BP"
data[grepl("CC", data[,1]),1] <- "CC"
data[grepl("MF", data[,1]),1] <- "MF"
data<-data[,c(2,3,6,10,5,12,1)]
colnames(data)<-c("Term","Gene_Number","Genes","Fold_Enrichment_Score","pvalue","adj-p","Category")
data$Term<-factor(data$Term,
                  levels=data$Term[order(data$Category, -data$pvalue)])

# Stop if not GO terms
if(dim(data)[1]==0){
  stop("No significant GO terms. No plots or files will be produced")
}

# Write to CSV
write.csv(data,paste0(args[5],"_DEG_GO_term_enrichment.csv"))

# Make the plot
ggplot(data[data$pvalue<0.05,],
       aes(x= Fold_Enrichment_Score,
           y=Term,
           size=Gene_Number,
           colour=pvalue))+
  facet_grid(Category~.,
             scales="free_y",
             space="free_y")+
  geom_point()+
  scale_colour_gradient(low="red",
                        high="blue",
                        name="p-value")+
  scale_size(name="Num. Genes")+
  xlab("Fold Enrichment Score")+
  ylab("GO Term")+
  ggtitle("Enriched GO Terms")+
  theme_minimal_grid() +
  theme(strip.background = element_rect(
    color="black", 
    fill="#CCCCCC", 
    size=1, 
    linetype="solid"))

# Save
ggsave2(filename = paste0(args[5], "_DEG_GO_term.pdf"),
        width=2+max(nchar(as.character(data[data$pvalue<0.05,1])))*0.1,
        height=1+dim(data)[1]*.25,
        limitsize = FALSE)

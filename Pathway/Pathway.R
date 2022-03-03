# This script is a stand-alone pathway analysis script that takes several command arguments
# you should first activate the RNAtier2 conda environment to make sure all of the necessary packages are loaded

### Example
# ./Pathway.R mouse DEG_list.csv Path/To/Results



# Error Checking ----------------------------------------------------------
# Conda environment
if(system("which R", intern = T)!="/home/genomics/anaconda3/envs/RNAtier2/bin/R"){
  message("WARNING:\nYou should activate the RNAtier2 conda environment prior to running this script.")
}

## Command Args
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("Speficy the input DEG list, species, and (optionally) a results path.", call.=FALSE)
}

# Set up ------------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse, quietly = T)
  library(clusterProfiler, quietly = T)
  library(magrittr, quietly = T)
  library(cowplot, quietly = T)
})

### Check species options
if(args[1]=="mouse"){
  suppressPackageStartupMessages(require(org.Mm.eg.db, quietly=T))
  message("Species set as ", args[1])
  assign("orgDB", org.Mm.eg.db)
  meta <- c(org.Mm.eg_dbInfo()[org.Mm.eg_dbInfo()[,1]=="ORGANISM",2], 
            org.Mm.eg_dbInfo()[org.Mm.eg_dbInfo()[,1]=="GOEGSOURCEDATE",2],            
            org.Mm.eg_dbInfo()[org.Mm.eg_dbInfo()[,1]=="KEGGSOURCEDATE",2],
            "mmu")
} else if(args[1]=="human"){
  suppressPackageStartupMessages(require(org.Hs.eg.db, quietly=T))
  message("Species set as ", args[1])
  assign("orgDB", org.Hs.eg.db)
  meta <- c(org.Hs.eg_dbInfo()[org.Hs.eg_dbInfo()[,1]=="ORGANISM",2], 
            org.Hs.eg_dbInfo()[org.Hs.eg_dbInfo()[,1]=="GOEGSOURCEDATE",2],            
            org.Hs.eg_dbInfo()[org.Hs.eg_dbInfo()[,1]=="KEGGSOURCEDATE",2],
            "hsa")
} else if(args[1]=="rat"){
  suppressPackageStartupMessages(require(org.Rn.eg.db, quietly = T))
  message("Species set as ", args[1])
  assign("orgDB", org.Rn.eg.db)
  meta <- c(org.Rn.eg_dbInfo()[org.Rn.eg_dbInfo()[,1]=="ORGANISM",2], 
            org.Rn.eg_dbInfo()[org.Rn.eg_dbInfo()[,1]=="GOEGSOURCEDATE",2],            
            org.Rn.eg_dbInfo()[org.Rn.eg_dbInfo()[,1]=="KEGGSOURCEDATE",2],
            "rno")
} else if(args[1]=="dog"){
  suppressPackageStartupMessages(require(org.Cf.eg.db, quietly = T))
  message("Species set as ", args[1])
  assign("orgDB", org.Cf.eg.db)
  meta <- c(org.Cf.eg_dbInfo()[org.Cf.eg_dbInfo()[,1]=="ORGANISM",2], 
            org.Cf.eg_dbInfo()[org.Cf.eg_dbInfo()[,1]=="GOEGSOURCEDATE",2],            
            org.Cf.eg_dbInfo()[org.Cf.eg_dbInfo()[,1]=="KEGGSOURCEDATE",2],
            "cfa")
} else {
  stop("Species ", args[1], " not recognized. Specify human, mouse, rat, or dog.")
}

### Check DEG format
tryCatch(DEG <- read_csv(args[2],
                         col_types = cols()),
         error=function(e){
           stop("Provide a valid CSV file of DEGs.")
         })
if(any(grepl(pattern = "_", x = DEG[1:5,]))){
  message("Parsing concatenated gene names into gene symbols")
  DEG %<>% 
    mutate(Row.names=str_split(Row.names, pattern = "_", simplify = T)[,2])
}
message("Note: DEG list will be used as-is. No filtering of DEGs (e.g. by p-value) will be applied.")

### Check/set results path
if(is.na(args[3])){
  message("Results path not set explicitly.\nOutputting results to ./Pathway/")
  dir.create("./Pathway", showWarnings = F)
  outdir <- "./Pathway/"
} else {
  if(str_sub(args[3],str_length(args[3]))!="/"){
    outdir <- paste0(args[3],"/")
  } else {
    outdir <- args[3]
  }
  message("Outputting results to ",outdir)
  dir.create(outdir,showWarnings = F)
}

# Run GO analysis ------------------------------------------------------------

g1 <- DEG %>% 
  dplyr::select(1) %>% 
  deframe() %>% 
  enrichGO(ont = "ALL",
           OrgDb = orgDB,
           keyType = "SYMBOL",
           pvalueCutoff = Inf,
           qvalueCutoff = Inf) %>% 
  data.frame() 

# Run KEGG analysis -------------------------------------------------------
k1 <- DEG %>% 
  dplyr::select(1) %>% 
  rename_at(1, ~"Row.names") %>% 
  mutate(Row.names=AnnotationDbi::mapIds(orgDB, keys = Row.names, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")) %>% 
  drop_na() %>% 
  deframe() %>% 
  enrichKEGG(organism = meta[4]) %>% 
  setReadable(OrgDb = orgDB, keyType = "ENTREZID") %>% 
  data.frame()

# Write CSV ---------------------------------------------------------------
if(dim(g1)[1]>0){
  g2 <- g1 %>% 
    rowwise() %>% 
    mutate(Fold_Enrichment_Score=(eval(parse(text=GeneRatio)) / eval(parse(text=BgRatio)))) %>% 
    dplyr::arrange(pvalue) %>% 
    dplyr::select(ONTOLOGY, ID, Description, Fold_Enrichment_Score, Count, pvalue, p.adjust, geneID) %T>% 
    write_csv(paste0(outdir,"GO.csv"))
} else {
  message("No enriched GO terms found.\nNo GO plots or lists will be generated.")
}

if(dim(k1)[1]>0){
  k2 <- k1 %>% 
    rowwise() %>% 
    mutate(Fold_Enrichment_Score=(eval(parse(text=GeneRatio)) / eval(parse(text=BgRatio)))) %>% 
    dplyr::arrange(pvalue) %>% 
    dplyr::select(ID, Description, Fold_Enrichment_Score, Count, pvalue, p.adjust, geneID) %T>% 
    write_csv(paste0(outdir,"KEGG.csv"))
} else {
  message("No enriched KEGG terms found.\nNo KEGG plots or lists will be generated.")
}

# Plot --------------------------------------------------------------------
## GO
if(exists("g2")){
  # raw pvalue
  g3 <- g2 %>% 
    dplyr::filter(pvalue<0.05) %>% 
    mutate(Concatenated = paste0(ID,"~",Description)) %>% 
    mutate(Concatenated = str_trunc(Concatenated, 60)) %>% 
    mutate(Concatenated = factor(Concatenated,
                                 levels = .$Concatenated[order(-.$pvalue)])) 
  if(dim(g3)[1]>0){
    g3 %>% ggplot(aes_string(x="Fold_Enrichment_Score",
                             y="Concatenated",
                             size="Count",
                             colour="pvalue"))+
      facet_grid(ONTOLOGY~.,
                 scales="free_y",space="free_y")+
      geom_point()+
      scale_colour_gradient(low="red",
                            high="blue",
                            name="Raw p-value")+
      scale_size(name="Num. Genes")+
      labs(x="Fold Enrichment",
           y="Term",
           title="Enriched GO Terms",
           caption=bquote(italic(.(meta[1]))*', Sourced:'~.(meta[2]))) +
      theme_minimal_grid() +
      theme(strip.background = element_rect(
        color="black",
        fill="#CCCCCC",
        size=1,
        linetype="solid"))
    ggsave2(filename = paste0(outdir,"GO_p0.05.pdf"),
            height=dim(g3)[1]*0.25+2.5,
            width=10,
            limitsize = F)
  } else {
    message("No significant GO terms found at p<0.05.\nNo plots will be generated.")
    dev.off()
  }
  # adjusted pvalue
  g4 <- g2 %>% 
    dplyr::filter(p.adjust<0.05) %>% 
    mutate(Concatenated = paste0(ID,"~",Description)) %>% 
    mutate(Concatenated = str_trunc(Concatenated, 60)) %>% 
    mutate(Concatenated = factor(Concatenated,
                                 levels = .$Concatenated[order(-.$p.adjust)])) 
  if(dim(g4)[1]>0){
    g4 %>% 
      ggplot(.,aes_string(x="Fold_Enrichment_Score",
                           y="Concatenated",
                           size="Count",
                           colour="p.adjust"))+
          facet_grid(ONTOLOGY~.,
                     scales="free_y",space="free_y")+
          geom_point()+
          scale_colour_gradient(low="red",
                                high="blue",
                                name="Adj. p-value")+
          scale_size(name="Num. Genes")+
          labs(x="Fold Enrichment",
               y="Term",
               title="Enriched GO Terms",
               caption=bquote(italic(.(meta[1]))*', Sourced:'~.(meta[2]))) +
      theme_minimal_grid() +
          theme(strip.background = element_rect(
            color="black",
            fill="#CCCCCC",
            size=1,
            linetype="solid"))
    ggsave2(filename = paste0(outdir,"GO_padj0.05.pdf"),
            height=dim(g4)[1]*0.25+2.5,
            width=10,
            limitsize = F)
  } else
    message("No significant GO terms found at adjusted p<0.05.\nNo plot will be generated.")
}

## KEGG
## GO
if(exists("k2")){
  # raw pvalue
  k3 <- k2 %>% 
    dplyr::filter(pvalue<0.05) %>% 
    mutate(Concatenated = paste0(ID,"~",Description)) %>% 
    mutate(Concatenated = str_trunc(Concatenated, 60)) %>% 
    mutate(Concatenated = factor(Concatenated,
                                 levels = .$Concatenated[order(-.$pvalue)])) 
  if(dim(k3)[1]>0){
    k3 %>% ggplot(aes_string(x="Fold_Enrichment_Score",
                             y="Concatenated",
                             size="Count",
                             colour="pvalue"))+
      geom_point()+
      scale_colour_gradient(low="red",
                            high="blue",
                            name="Raw p-value")+
      scale_size(name="Num. Genes")+
      labs(x="Fold Enrichment",
           y="Term",
           title="Enriched KEGG Terms",
           caption=bquote(italic(.(meta[1]))*', Sourced:'~.(meta[3]))) +
      theme_minimal_grid() +
      theme(strip.background = element_rect(
        color="black",
        fill="#CCCCCC",
        size=1,
        linetype="solid"))
    ggsave2(filename = paste0(outdir,"KEGG_p0.05.pdf"),
            height=dim(k3)[1]*0.25+2.5,
            width=10,
            limitsize = F)
  } else {
    message("No significant KEGG terms found at p<0.05.\nNo plots will be generated.")
    dev.off()
  }
  # adjusted pvalue
  k4 <- k2 %>% 
    dplyr::filter(p.adjust<0.05) %>% 
    mutate(Concatenated = paste0(ID,"~",Description)) %>% 
    mutate(Concatenated = str_trunc(Concatenated, 60)) %>% 
    mutate(Concatenated = factor(Concatenated,
                                 levels = .$Concatenated[order(-.$p.adjust)])) 
  if(dim(k4)[1]>0){
    k4 %>% 
      ggplot(.,aes_string(x="Fold_Enrichment_Score",
                          y="Concatenated",
                          size="Count",
                          colour="p.adjust"))+
      geom_point()+
      scale_colour_gradient(low="red",
                            high="blue",
                            name="Adj. p-value")+
      scale_size(name="Num. Genes")+
      labs(x="Fold Enrichment",
           y="Term",
           title="Enriched KEGG Terms",
           caption=bquote(italic(.(meta[1]))*', Sourced:'~.(meta[3]))) +
      theme_minimal_grid() +
      theme(strip.background = element_rect(
        color="black",
        fill="#CCCCCC",
        size=1,
        linetype="solid"))
    ggsave2(filename = paste0(outdir,"KEGG_padj0.05.pdf"),
            height=dim(k4)[1]*0.25+2.5,
            width=10,
            limitsize = F)
  } else
    message("No significant KEGG terms found at adjusted p<0.05.\nNo plot will be generated.")
}

# No idea why this is create, but this will remove it
if(file.exists("./Rplots.pdf")){
  file.remove("./Rplots.pdf")
}


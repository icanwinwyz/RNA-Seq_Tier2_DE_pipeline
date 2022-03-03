# This should be run on Titan in the RNAtier2 Conda environment. 
# The above shebang *should* work for running this R snippet, but it will not work for the later GSEA script.

### Example
# ./GSEA_Prep.R DEG_list.csv Path/To/Results

# Error Checking ----------------------------------------------------------
# Conda environment
if(system("which R", intern = T)!="/home/genomics/anaconda3/envs/RNAtier2/bin/R"){
  message("WARNING:\nYou should activate the RNAtier2 conda environment prior to running this script.")
}

# Command Args
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("Speficy the input DEG list, species, and (optionally) a results path.", call.=FALSE)
}

# Check/set results path
if(is.na(args[2])){
  message("Results path not set explicitly.\nOutputting results to ./Pathway/")
  dir.create("./Pathway", showWarnings = F)
  outdir <- "./Pathway/"
} else {
  if(str_sub(args[2],str_length(args[2]))!="/"){
    outdir <- paste0(args[2],"/")
  } else {
    outdir <- args[2]
  }
  message("DEGs will be exported to ",outdir)
  dir.create(outdir,showWarnings = F)
}

message("Note: DEG list will be used as-is. No filtering of DEGs (e.g. by p-value) will be applied.")

# Set up ------------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse, quietly = T)
  library(magrittr, quietly = T)
})

# Reformat ------------------------------------------------------------------
### Check DEG format and read in
tryCatch(DEG <- read_csv(args[1],
                         col_types = cols()),
         error=function(e){
           stop("Provide a valid CSV file of DEGs.")
         })

if(any(grepl(pattern = "_", x = DEG[1:5,]))){
  message("Parsing concatenated gene names into gene symbols")
  DEG %<>% 
	tidyr::separate(Row.names, sep = "_",
		into = c("ENSG", "Symbol"),
                extra="merge") %>%
	dplyr::select(Symbol,log2FoldChange) %>%
	dplyr::arrange(-log2FoldChange) %>%
	readr::write_tsv(paste0(outdir,"GSEA_Export.rnk"),
		col_names = F,
		quote = "none")
}

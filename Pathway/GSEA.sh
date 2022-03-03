#!/bin/bash

# Example ./GSEA.sh human ./GSEA_Export.rnk ./Path/to/results/directory/

# set env variable for program
source activate RNAtier2
GSEA_PATH=/home/genomics/apps/GSEA/
export GSEA=$GSEA_PATH/v4.2.2/gsea-cli.sh
if [ $0 != "human" ]; then
	stop("At present only GSEA of human samples is supported")
else
	export msigdb="$GSEA_PATH/Datasets/c3.all.v7.5.1.symbols.gmt.txt,$GSEA_PATH/Datasets/c5.all.v7.5.1.symbols.gmt.txt,$GSEA_PATH/Datasets/h.all.v7.5.1.symbols.gmt.txt"
fi

# set analysis ins and outs
export RNK_IN=$2
export PATH_OUT=$3

# Run GSEA
$GSEA GSEAPreranked \
  -gmx $msigdb \
  -collapse No_Collapse \
  -mode Abs_max_of_probes \
  -norm meandiv \
  -nperm 1000 \
  -rnd_seed timestamp \
  -rnk $RNK_IN \
  -scoring_scheme weighted \
  -rpt_label GSEA_Results \
  -create_svgs false \
  -include_only_symbols true \
  -make_sets true \
  -plot_top_x 20 \
  -set_max 500 \
  -set_min 15 \
  -zip_report true \
  -out $PATH_OUT



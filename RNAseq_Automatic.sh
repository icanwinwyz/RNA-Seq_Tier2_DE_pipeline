#########################################################################
# File Name: RNAseq_Automatic.sh
# Author: Di Wu
# mail: di.wu@cshs.org
# Created Time: Feb., 2019
#########################################################################
#!/bin/bash
#!/bin/bash
display_usage() {
	echo -e "NAME:\n  RNAseq_analysis."
	echo -e "\nDESCRIPTION:\n   This pipeline will use COUNT file to do downstream differentially expressed genes (DEGs) analysis."
	echo -e "\nUsage:\n   bash $TIER2/RNAseq_Automatic.sh Count_file.csv Sample_info.csv comparisons.csv "
    
    echo "Input options:"
    echo "   -h|--help    show this help"
    
    echo "Input files:"
    echo "   The first input is Count_file.csv, check https://github.com/dxw5099/RNAseq_downstream/blob/master/demo_COUNTS.csv  as an example format"
    echo "   The second input is Sample_info.csv, check https://github.com/dxw5099/RNAseq_downstream/blob/master/demo_sample_info.csv as an example format"
    echo "   The third input is Conparisons.csv, check https://github.com/dxw5099/RNAseq_downstream/blob/master/demo_comparisons.csv as an example format"
    echo "   The fourth input is project_ID, for example: AP-5782–11–08–2018"
	}
# check whether user had supplied -h or --help . If yes display usage
if [[ ($1 == "--help") ||  ($1 == "-h") ]]
then
	display_usage
	exit 0
fi

# get PCA plots for all samples and DEGs table, interactive report for each comparison
Rscript $TIER2/Automatic/RNAseq_tier2.R $1 $2 $3 $4

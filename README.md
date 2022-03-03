# Set Up
All of the analyses in this repository should be *run on Titan* and they all require you to first activate the RNATier2 conda environment to ensure that all of the required R packages are loaded.
```bash
source activate RNAtier2
```

## DEG Analysis

You can get help text for any of the DEG functions with the `-h` flag.
```bash
$ source activate RNAtier2
$ bash $TIER2/RNAseq_Automatic.sh -h 
$ bash $TIER2/RNAseq_Fixed_Cutoff.sh -h
$ bash $TIER2/RNAseq_FFPE.sh -h
```
### Automatic Cutoff

Most typically, we run the `RNAseq_Automatic.sh` script. This script requires 4 inputs:

- A count matrix produced by the Tier1 RNAseq pipeline, where rows are genes and columns are samples.
- A csv file of sample information containing (in order)  sample_ID, Sample_Name, and Group_info for each sample to be processed. 
  - These sample IDs *and their order* must match those in the columns of the counts matrix. 
  - The Group_info must be a part of the Sample_Name, since it will be used to grep which samples to include in the analysis. For example, if Group_info is G1 then Sample_Name should be G1_S1 or S1_G1. If Sample_Name is S1_Group1 or Group1_S1, no results will be generated. 
- A csv file of comparisons to be made using the same Group_info names as the previous sample information sheet. The order of groups in the comparison will inform the calculation of log<sub>2</sub> fold change of expression.
- A string specifying the project ID to name a results folder.

```bash
$ source activate RNAtier2
$ bash $TIER2/RNAseq_Automatic.sh \
  counts.csv \
  sample_info.csv \
  comparisons.csv \
  project_ID
```

This will produce several outputs:
```
./
├── Group1_vs_Group2
│   ├── Group1_vs_Group2.html
│   ├── Group1_vs_Group2_DEGs_all.csv
│   ├── Group1_vs_Group2_DEGs_padj0.05.csv
│   ├── Group1_vs_Group2_DEGs_padj0.01.csv
│   ├── Group1_vs_Group2_DEGs_padj0.01_FC2.csv
│   ├── Group1_vs_Group2_DEGs_padj0.1.csv
│   ├── Group1_vs_Group2_DEGs_p0.05.csv
│   ├── Group1_vs_Group2_heatmap.pdf
│   ├── Group1_vs_Group2_MAplot.pdf
│   ├── Group1_vs_Group2_PCA_top500.pdf
│   └── Group1_vs_Group2_Volcano.pdf
├── project_ID_PCA_all_samples.pdf
```
The number of folders will match the number of comparisons specified in the `comparison.csv` file. The specific DEG lists output will depend on the number of DEGs in each category. The script will always output the `_all.csv` files but will create different filtered lists to present a manageable number of DEGs. The folder will also contain a heatmap, MA plot, PCA, and volcano plot. Finally, the folder will contain a summary html file with all of the information.

In the main project folder there is also a PCA plot with all samples include. In contrast the the PCA plots in the comparison-specific folders, the PCA plot here will show all samples from the `sample_info.csv` sheet,even if they are not included in any of the comparisons listed in the `comparisons.csv` file

### Fixed Cutoff

Running the pipeline with a fixed cutoff for filtering is very similar to the automatic cutoff outlined above. The key difference is that it requires you to specify hard cutoffs for *both* adjusted p-value and log<sub>2</sub> fold-change.

For example, to filter by an adjusted p-value less than 0.05 and a log<sub>2</sub> fold-change greater than 1:
```bash
$ source activate RNAtier2
$ bash $TIER2/RNAseq_Fixed_Cutoff.sh \
  count.csv \
  sample_info.csv \
  comparison.csv \
  project_ID \
  "p_adj<0.05" \
  1
```

Note, the the adjusted p-value cutoff must be quoted and does not accept other criteria to filter on. The outputs for this version of the pipeline are identical to the automatic cutoff, except that the filter DEG lists will only include the pre-set criteria.


### Human FFPE samples

Running the pipeline for formalin-fixed paraffin-embedded (FFPE) samples is very similar to running it with the automatic cutoffs. Currently this version of the pipeline only works on human samples since it depends on conversion of ENSEMBL IDs for gene symbols.

```bash
$ source activate RNAtier2
$ bash $TIER2/RNAseq_FFPE.sh \
  count.csv \
  sample_info.csv \
  comparison.csv \
  project_ID
```

## Pathway Analysis

We have upgraded our pathway analysis  to move away from [DAVID](https://david.ncifcrf.gov/) and [GSEA](https://www.gsea-msigdb.org/gsea/index.jsp) in favor of [ClusterProfileR](https://yulab-smu.top/clusterProfiler-book/), an R package that draws on other more frequently updated databases and is also more high throughput.

For projects that were started before May 2021, please continue to use DAVID to ensure continuity with previous results. For all newer projects, use the ClusterProfileR analysis pipelines.

### ClusterProfileR GO and KEGG Enrichment

Currently, the ClusterProfileR pipeline can perform GO and KEGG enrichments from four species (human, mouse, rat, and dog). It was designed to take the any of the DEG lists from our standard RNAtier2 pipeline as inputs, but it can also operate on on csv file contained gene symbols in the first column. Note that this pathway analysis pipeline does not filter the input DEG list, so passing in the `_DEGs_all.csv` output is likely to produce noise results if there are many genes.

The pipeline is run as follows:

```bash
$ source activate RNAtier2
$ Rscript $TIER2/Pathway/Pathway.R \
  species \
  DEG_List.csv \
  Path/To/Results
```

- The species must be specified in lower-case and can be human, mouse, rat, or dog currently.
- The DEG_list.csv must be in comma-separated format and contain the gene symbols in the first column.
- The path to the results will be created if it does not yet exist, and will overwrite existing results if present.

The pipeline produces up to five output files:

```
Path/To/Results/
├── GO.csv
├── GO_p0.05.pdf
├── GO_padj0.05.pdf
├── KEGG.csv
├── KEGG_p0.05.pdf
└── KEGG_padj0.05.pdf
```

The `GO.csv` and `KEGG.csv` files will contain the unfiltered list of term names, ID, ontologies (for GO terms), fold enrichment scores, gene names, and gene counts. For both GO and KEGG, the pipeline will produce a dot plot filtered by either raw p-value<0.05 `p0.05.pdf` or adjusted p-value<0.05 `padj0.05.pdf`. If any of these files would be blank or contain no terms, the pipeline will print a warning and not create the file.

### ClusterProfileR GSEA

Although GSEA can be run with any set of genes, we currently only support GSEA of GO terms. In contrast to the GO and KEGG enrichment analysis, GSEA also requires information on log fold change of expression in order to rank the genes. This script replaces the [GSEA software](https://www.gsea-msigdb.org/gsea/index.jsp) we had used previously. The pipeline will attempt to locate the column, searching for the first column whose name contains the word `log`. This can be overridden by including a column index number instead.

```bash
$ source activate RNAtier2
$ Rscript $TIER2/Pathway/GSEA.R \
  species \
  DEG_List.csv \
  Pathway/To/Results \
  2 # optional column index for log fold change
```

The pipeline produces up to three output files:

```
Path/To/Results/
├── GSEA.csv
├── GSEA_p0.05.pdf
└── GSEA_padj0.05.pdf
```


### DAVID Figure Generation (Prior to May 2021)

After the GO and KEGG results have been downloaded from the DAVID webserver, run the figure generation script as follows:

```bash
$ source activate RNAtier2
$ Rscript $TIER2/Pathway/Pathway_DAVID.R \
  KEGG.txt \
  GO_BP.txt \
  GO_CC.txt \
  GO_MF.txt \
  project_ID
```

This will automatically generate the GO and KEGG enrichment dot plots and write lists of the enriched terms to csv files. The csv files will contain all terms regardless of significance, but the plots will be filtered to only contain those terms with a raw p-value less than 0.05.

```
project_ID/
├── project_ID_DEG_GO_term.pdf
├── project_ID_DEG_GO_term_enrichment.csv
├── project_ID_DEG_KEGG.pdf
└── project_ID_DEG_KEGG_enrichment.csv
```

# Example

```bash
# activate the conda environment
$ source activate RNAtier2
# Navigate to output folder
$ cd Demo
$ cd mkdir Results
$ cd Results
# Run DEG analysis
$ bash $TIER2/RNAseq_Automatic.sh \
  ../demo_COUNTS.csv \
  ../demo_sample_info.csv \
  ../demo_comparisons.csv \
  demo_results
# Run pathway analysis
$ cd Demo/Results/M_vs_MC
$ Rscript $TIER2/Pathway/Pathway.R \
  human \
  M_vs_MC_DEGs_padj0.05.csv \
  Pathway_Results
# Run GSEA analysis
$ Rscript $TIER2/Pathway/GSEA.R \
  human \
  M_vs_MC_DEGs_padj0.05.csv \
  Pathway_Results
```

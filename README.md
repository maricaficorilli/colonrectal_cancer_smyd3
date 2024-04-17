## Tailoring a novel colorectal cancer stem cell-targeted therapy by inhibiting the methyltransferase SMYD3

The following script utilizes samples of cell lines (HCT-116) as input to conduct a differential expression analysis (DEA) using edgeR Bioconductor/R package (version 3.42.4), 
followed by an enrichment analysis on DE genes using clusterProfiler Bioconductor/R package (version 4.8.3). 

In DEA-Results folder, there are additional files on the DEA results of the various comparisons. In this case we used a log fold change (logFC) >= 1 or <= -1 and False Discovery Rate (FDR) < 0.01.

NOTE: It is necessary to download the raw genes matrix from the following code GSE264157 on GEO -NCBI (https://www.ncbi.nlm.nih.gov/geo/) before using the script.

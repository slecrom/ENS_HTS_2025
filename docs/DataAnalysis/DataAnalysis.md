# Data analysis

1. [Introduction](#introduction)
2. [Search for differentially expressed genes](#search-for-differentially-expressed-genes)
3. [Functional analyses of differentially expressed genes](#functional-analyses-of-differentially-expressed-genes)


***

## Introduction

!!! caution "Objective of this practical session"

	During this practical session, you will learn:	

	* To perform statistical analysis of the gene expression matrix in order to identify differentialy expressed genes between two conditions.
	
	* To find biological relevant target genes

***

## Search for differentially expressed genes

In their article (Guida et al., 2011), the authors repeated the experiment 6 times for normoxic condition (with O2) and 4 times for hypoxic conditions (without O2). 

!!! tip "What you have to do"
	- Search for differentially expressed genes using limma and DESeq2 R packages for microarrays and RNAseq experiments.
	- How many genes are selected with different adjusted p-value thresholds (5%, 1%, etc.)?

### Set up your working environment

Connect to Rstudio server of the IFB. Look at the [tutorial on how to connect to IFB-core Rstudio server](../IFB_OpenOnDemand.md) to see how to proceed.

### Save the working notebook in your personal environment

1. In "*File > Open File...*" enter the path `/shared/projects/2528_ens_master2lf_fgat/data/tutorials/data_analysis.Rmd` to open the notebook containing all the code needed for the practical.

&nbsp;
2. Save it into your personal folder on your IFB account using "*File > Save As*"
   
### Follow the instruction of the notebook to conduct the analysis

You can also visualize the final [report version](data_analysis_report.html).

You can find help on how to use R markdown on the [R markdown project webpage](https://rmarkdown.rstudio.com/lesson-2.html).

***

## Functional analyses of differentially expressed genes

Several tools are available online to evaluate the biological relevance of the gene sets you select after the differential analysis. For example you can you use the [GoTermFinder tool](<http://www.candidagenome.org/cgi-bin/GO/goTermFinder>) dedicated to *Candida* yeast species to retrieve functional annotation. Be careful to choose the genome of *Candida parapsilosis* in the species selection dialog box. You can also obtain information on a specific gene using the [Candida Genome Database](http://www.candidagenome.org/).

!!! tip "What you have to do"
	- What are the functions of the genes in your lists (identified at the previous step).
	- Are they relevant with the studied biological system described by Guida et al. in their publication?
	- Compare your results with those presented in the original publication

To go further with functional analyses, you can read the [Nine quick tips for pathway enrichment analysis article](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1010348)
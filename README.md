---
title: "Differential methylation analysis of MBDcap-seq dataset using MBDDiff package"
author: "Yuanhang Liu"
date: "`r Sys.Date()`"
vignette: >
  %\VignetteIndexEntry{Differential methylation analysis of MBDcap-seq dataset using MBDDiff package}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8} 
output: 
  BiocStyle::html_document:
    toc: true
---

```{r style, echo = FALSE, results = 'asis'}
BiocStyle::markdown()
```

## Introduction 

Methyl-CpG binding domain-based capture
followed by high throughput sequencing (MBDCap-seq) is widely
used to examine DNA methylation pattern genome-wide. Current
MBDCap-seq data analysis approaches focus on measurement of
methylated CpG sequence reads, without considering genomic
characteristics and tissue-specific context and their impact to the
amount of methylated DNA measurement (signal) and background
fluctuation (noise). Therefore, specific software needs to be
developed to process MBDCap-seq datasets. Here we presented a
novel algorithm, termed MBDDiff, implemented as an R package that
is designed specifically for processing MBDCap-seq datasets.
MBDDiff contains three modules: quality assessment of datasets and
quantification of DNA methylation; determination of differential
methylation of promoter regions; and visualization functionalities.

## Installation 

Currently, MBDDiff can be installed from Github by: 
```{r,eval=FALSE}
install.packages('devtools')
library(devtools)
install_github('Liuy12/MBDDiff')
install_github('Liuy12/XBSeq')
```
```{r,message = FALSE, warning=FALSE, eval=FALSE}
library("MBDDiff")
```
We will submit the package to Bioconductor at earliest time possible. 

## Processing MBDcap-seq datasets using MBDdiff

#### Retriving annotation files for promoter regions

Genomic annotation files can be downloaded from UCSC database and processed within MBDdiff to locate promoter regions. Isoforms with same promoter regions will only be counted once. For instance, to get promoter annotation with human genome build 'hg19', we will do,

```{r,message = FALSE, warning=FALSE, eval=FALSE}
Promoter_anno <- GetPromoterAnno('hg19')
```

We have already generated promoter annotation file for several organisms including, human 'hg19' and 'hg38' build, mouse 'mm10' and 'mm9' build, rat 'rn6' and 'rn5' build. These files can be downloaded from Github: MBDDiff_files

#### Counting reads mapped to promoter regions

Reads mapped to each promoter regions can be counted by using bedtools coveragebed function. We implemented a wrapper function in MBDDiff which basically calls bedtools internally. You will need to provide the path to '.bam' and '.bed' files:

```{r,message = FALSE, warning=FALSE, eval=FALSE}
Promoter_counts <- bedtoolsr(bamdir = '/path/to/bam/files', bed = '/path/to/bed/file')
```

The output is a list contains two elements, a matrix containing promoter mapped reads for each sample and a vector containing region length information. 

#### Identify background regions for MBDcap-seq

We have already carried out following steps to identify background regions for several organisms including, human 'hg19' and 'hg38' build, mouse 'mm10' and 'mm9' build, rat 'rn6' and 'rn5' build. These files can be downloaded from Github: MBDDiff_files. 

If you need to or want to learn how to construct background regions for your organism of interest, please carry out following steps.

To identify regions that potentially contribute to background noise, we built a 100 bp tiling window across the whole genome. In order to identify the regions for measuring background noise, we applied following procedures:

1. Filtering step: In all the 100 bp windows, exclude any window that resides in promoter regions, predicted CpG islands regions and windows that contains ambiguous bases (gaps)

Firstly, construct potentially functional regions to exclude from 100bp windows

Then, extract sequence files for the newly constructed 100bp windows. The output fasta file will be used to calcualte GC content statistics in whole genome sliding windows:

Finally, filter 100 bp windows by excluding those that intersect with functional elements

2. Construct preliminary background regions based on GC
content: for each promoter region, identify 80 100bp
windows nearby that has low in GC content (< 40%) and
also relatively proximal to the corresponding TSS:

All the above procedures can be done via one function:

```{r,message = FALSE, warning=FALSE, eval=FALSE}
bed_100bp_bg <- IdentifyBackground(organism = 'hg19', bed_path = 'path/to/bed', binsize = 100, promo_bed = Promoter_anno, cores = 4)
```

Details regarding each this function can be found in the manual page. 

#### Counting reads mapped to background regions

We will follow a similar procedure as we did for promoter regions:

```{r,message = FALSE, warning=FALSE, eval=FALSE}
Background_counts <- bedtoolsr(bamdir = '/path/to/bam/files', bed = '/path/to/bg_bed/file')
```

#### Summerize background mapped reads to each promoter

For each promoter region, choose 40 out of
80 100 bp windows that are relatively low in TPM (Transcript Per Million)

```{r,message = FALSE, warning=FALSE, eval=FALSE}
Background_counts_summe <- Summerizebg(Background_counts[[1]], Background_counts[[2]])
```

#### Test for differential methylation 

```{r,message = FALSE, warning=FALSE, eval=FALSE}
Condition <- c(rep('C1', 3), rep('C2', 3))
TestStat <- MBDDiff(Promoter_counts, Background_counts_summe, Condition)
```

## Visualization functionalities of MBDDiff

We provide several visualization functionalities to help you visulize sample relationships based on methylation profile; distribution of promoter vesus background mapped reads, etc. 

#### Enrichment of GC content grouped by methylation levels

Do this if you haven't already done so:

```{r,message = FALSE, warning=FALSE, eval=FALSE}
whole_genome_bin_count <- bedtoolsr(bed, bam)
whole_genome_bin_count_TPM <- CalculateTPM(whole_genome_bin_count[[1]], whole_genome_bin_count[[2]])
```

Then, we can generate enrichment of GC content grouped by methylation levels: 

```{r,message = FALSE, warning=FALSE, eval=FALSE}
MethyEnrich(whole_genome_bin_count_TPM, fa)
```

#### Distribution of promoter counts and background noise

Extract XBSeqDataSet object from MBDDiff function and then plot the distribution using XBplot:

```{r,message = FALSE, warning=FALSE, eval=FALSE}
MBD <- TestStat[[1]]
XBplot(MBD, Samplenum = 1, unit = 'LogTPM', Genelength = Promoter_counts[[2]])
```

#### 3d principal component analysis

Firstly extract normalized reads and then generate 2d/3d pca plot:

```{r,message = FALSE, warning=FALSE, eval=FALSE}
Norm_count <- counts(MBD, normalized = TRUE)
pcaplot(Norm_count, cv.Th = 0.1, method = 'pca', dimension = c(1,2,3) )
```

#### MAplot of promoters

```{r,message = FALSE, warning=FALSE, eval=FALSE}
MAplot(TestStat[[2]])
```

#### Heatmap of selected promoters

```{r,message = FALSE, warning=FALSE, eval=FALSE}
DE_index <- with(TestStat[[2]], which(baseMean > quantile(baseMean)[2] & 
                  abs(log2FoldChange) > 1 &
                  padj < 0.1))
heatmap.3(Norm_count[DE_index,])
```

#### Dynamic visualization

For all figures generated above, there is one more option that the user can specify to generate one equivalent dynamic figure. For instance to generate the equivalent dynamic figure of methylation enrichment, you can do:

```{r,message = FALSE, warning=FALSE, eval=FALSE}
MethyEnrich(whole_genome_bin_count_TPM, fa, interactive = TRUE)
```

## Bug reports
Report bugs as issues on our [GitHub repository](https://github.com/Liuy12/MBDDiff/issues) or you can report directly to my email: liuy12@uthscsa.edu.
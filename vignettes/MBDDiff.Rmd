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
Promoter_counts <- bedtoolsr('/path/to/bam/files', '/path/to/bed/file')
```

The output is a list contains two elements, a matrix containing promoter mapped reads for each sample and a vector containing region length information. 

#### Identify background regions for MBDcap-seq

We have already carried out following steps to identify background regions for several organisms including, human 'hg19' and 'hg38' build, mouse 'mm10' and 'mm9' build, rat 'rn6' and 'rn5' build. These files can be downloaded from Github: MBDDiff_files. 

If you need to or want to learn how to construct background regions for your experiment, please carry out following steps.

To identify regions that potentially contribute to
background noise, we built a 100 bp tiling window across the
whole genome:

```{r,message = FALSE, warning=FALSE, eval=FALSE}
bed_100bp <- CreateWindows(chromlen, binsize)
```

In order to identify the regions
for measuring background noise, we applied following
procedures:

1. Filtering step: In all the 100 bp windows, exclude any
window that resides in promoter regions, predicted CpG
islands regions and windows that contains ambiguous bases
(gaps)

Firsly, extract sequence files for the newly constructed 100bp windows:

```{r,message = FALSE, warning=FALSE, eval=FALSE}
Getfasta(fa, bed, path = 'path/to/save/file')
```

Then construct regions to exclude from 100bp windows:

```{r,message = FALSE, warning=FALSE, eval=FALSE}
Exclude_regions <- ConstructExRegions(organism, Promoter_anno) 
```

Finally, filter 100 bp windows by excluding those that intersect with functional elements

```{r,message = FALSE, warning=FALSE, eval=FALSE}
bed_100bp_filtered <- FilterRegions(bed_100bp, 'path/to/fasta', Exclude_regions) 
```

2. Construct preliminary background regions based on GC
content: for each promoter region, identify 80 100bp
windows nearby that has low in GC content (< 40%) and
also relatively proximal to the corresponding TSS:

```{r,message = FALSE, warning=FALSE, eval=FALSE}
bed_100bp_bg <- IdentifyBackground(bed_100bp_filtered, Promoter_anno)
```

#### Counting reads mapped to background regions

We will follow a similar procedure as we did for promoter regions:

```{r,message = FALSE, warning=FALSE, eval=FALSE}
Background_counts <- bedtoolsr('/path/to/bam/files', '/path/to/bg_bed/file')
```

#### Summerize background mapped reads to each promoter

For each promoter region, choose 40 out of
80 100 bp windows that are relatively low in TPM (Transcript Per Million)

```{r,message = FALSE, warning=FALSE, eval=FALSE}
Background_counts_summe <- Summerizebg(Background_counts, gaplength)
```

#### Test for differential methylation 

```{r,message = FALSE, warning=FALSE, eval=FALSE}
Condition <- c(rep('C1', 3), rep('C2', 3))
TestStat <- MBDDiff(Promoter_counts, Background_counts_summe, Condition)
```

#### Some visualization functionalities

We provide several visualization functionalities to help you visulize sample relationships based on methylation profile; distribution of promoter vesus background mapped reads, etc. 

## Bug reports
Report bugs as issues on our [GitHub repository](https://github.com/Liuy12/MBDDiff/issues) or you can report directly to my email: liuy12@uthscsa.edu.
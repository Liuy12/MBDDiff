\name{PromoterCount}
\alias{Background}
\alias{Promoter}
\alias{Promoterlength}
\title{
   Example Datasets used in the manual pages
}
\description{
   Example datasets used in manual pages by carrying out bedtools procedure for promoter regions (Promoter) and background regions (Background) and promoter length information (Promoterlength).
}
\usage{
   data(PromoterCount)
}
\format{
   ExampleData contains two data.frames . Two of them are expression matrix. One is called 'Observed'. One is called 'Background'. For the two data.frames, rows represent methylation levels for promoter regions. Columns represent each sample. Both the two data.frames have total of 33178 number of rows and 6 number
   of columns. There is also another numeric vector containing the promoter length information named as 'Promoterlength'
}
\details{
   Datasets were obtained from GEO database with accession number GSE58020. The authors carried out MBDCap-seq to investigate DNA methylation profile for tripe negative breast cancers (TNBCs). Details regarding how to generate background annotation file can be found in vignette
}
\value{
   Two data.frames and one numeric vector as described in format section. 
}
\author{
   Yuanhang Liu
}
\references{
https://github.com/Liuy12/MBDDiff
}
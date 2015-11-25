\name{GetPromoterAnno}
\alias{GetPromoterAnno}
\title{
Retrive promoter annotation from UCSC database
}
\description{
Retrive promoter annotation of user-specified organism from UCSC database
}
\usage{
GetPromoterAnno(organism, save = F, Dir = NULL)
}
\arguments{
  \item{organism}{
  the organism you want to extract promoter annotation
  }
  \item{save}{
  Logical value indicating whether you want to save the output locally
  }
  \item{Dir}{
  If save is set to TRUE, you also need to provide directory that you would like to save the output to. 
  }
}
\details{
  This function firstly retrive 'refflat' table information from UCSC database through their MySQL server. Then for each gene, records with the same TSS are collapsed. Genes with multiple records are marked by the same gene name followed by '#' and a number.
}
\value{
A data.frame containing promoter annotation information
\item{geneName}{Official gene symbols}
\item{chrom}{chromsome number}
\item{strand}{The strand of the gene}
\item{txStart}{Start location of transcripts}
\item{txEnd}{End location of transcripts}
}
\references{
  https://github.com/Liuy12/MBDDiff
}
\author{
  Yuanhang Liu
}
\examples{
  Promoter_anno <- GetPromoterAnno('hg19')
}
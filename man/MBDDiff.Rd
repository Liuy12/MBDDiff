\name{MBDDiff}
\alias{MBDDiff}
\title{
  Performing differential methylation analysis based on Poisson-Negative binomial convolution model
}
\description{
  This function is basically a wrapper function for \code{\link{XBSeq}}
}
\usage{
MBDDiff(promoter, background, conditions, method = "pooled", sharingMode = "maximum", fitType = "local", pvals_only = FALSE, paraMethod = "NP")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{promoter}{
  A data.frame or matrix that contains promoter methylation levels
  }
  \item{background}{
  A data.frame or matrix of the same dimension as promoter that contains background noise for each promoter region
  }
  \item{conditions}{
  A character or factor vector that contains the experimental design information
  }
  \item{method}{
  Method used to estimate SCV
  }
  \item{sharingMode}{
  Mode of sharing of information
  }
  \item{fitType}{
  Option to fit mean-SCV relation
  }
  \item{pvals_only}{
  Logical; Specify whether to extract pvalues only
  }
  \item{paraMethod}{
  Method to use for estimation of distribution parameters, 'NP' or 'MLE'
  }
}
\details{
  For further details please refer to \code{\link{XBSeq}}
}
\value{
A list of two elements:
\item{MBD}{
  A \code{\link{XBSeqDataSet}} object
}
A data.frame of the following elements:
\item{id}{
   rownames of XBSeqDataSet
}
\item{baseMean}{
   The basemean for all promoters
}
\item{baseMeanA}{
   The basemean for condition 'A'
}
\item{baseMeanB}{
   The basemean for condition 'B'
}
\item{foldChange}{
   The fold change compare condition 'B' to 'A'
}
\item{log2FoldChange}{
   The log2 fold change
}
\item{pval}{
   The p value for all promoters
}
\item{padj}{
   The adjusted p value for all promoters
}
}
\references{
  https://github.com/Liuy12/MBDDiff
}
\author{
Yuanhang Liu
}
\seealso{
  \code{\link{estimateRealCount}},  \code{\link{XBSeqDataSet}}, \code{\link{estimateSCV}}, \code{\link{XBSeqTest}}
}
\examples{
  \dontrun{
    data(PromoterCount)
    conditions <- c(rep('C1', 3), rep('C2', 3))
    MBDDiff(Promoter, Background, conditions)
    }
}
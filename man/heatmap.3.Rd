\name{heatmap.3}
\alias{heatmap.3}
\title{More customized heatmap that allow control over color key}
\description{
    More customized heatmap that is designed based on heatmap.2 that gives user control over color key
}
\usage{
heatmap.3(Exprs, sel = F, thres_mean, thres_var, numbreaks = 100, col = c("blue", "white", "red"), breakratio = c(2, 1, 2), colsidebar, Colv = F, Rowv = T, scale = "row", labRow = F, labCol = F, dendrogram = "row")
}
\arguments{
  \item{Exprs}{
  A numeric matrix or data.frame
  }
  \item{sel}{
  Whether you want to pre-filter the genes first based on average expression and variance
  }
  \item{thres_mean}{
  If 'sel' is set to TRUE, what threshold you wanna use for average expression
  }
  \item{thres_var}{
  If 'sel' is set to TRUE, what threshold you wanna use for variance
  }
  \item{numbreaks}{
  The number of breaks for colorkey. Default to 100
  }
  \item{col}{
  A vector of length 3 of colors you want to use for colorkey. Default to c("blue", "white", "red")
  }
  \item{breakratio}{
  How you wanna the colorkey to be partitioned according to the colors.
  }
  \item{colsidebar}{
  Whether you want to add a column side bar above the heatmap
  }
  \item{Colv}{
  Logical value indicating whether you want to perform clustering on columns
  }
  \item{Rowv}{
  Logical value indicating whether you want to perform clustering on rows
  }
  \item{scale}{
  Whether you want to scale the input matrix row-wise or column wise so that each row or column has mean equals to 0 and variance equals to 1. Possible values are 'row' and 'column'
  }
  \item{labRow}{
  Logical value indicating whether you want to label rows based on row names of input matrix
  }
  \item{labCol}{
   Logical value indicating whether you want to label columns based on column names of input matrix 
  }
  \item{dendrogram}{
  One of values in "both", "row", "column", "none". Indicating how you wanna the dendrogram to be drawn
  }
}
\details{
   heatmap.3 internally calls heatmap.2 from gplots package, which uses layout to arragent the plot elements. Consequentially, it can not be used in a multi column/row layout using layout(...), par(mfrow=...) or (mfcol=...).
}
\references{
    https://github.com/Liuy12/MBDDiff
}
\author{
    Yuanhang Liu
}
\seealso{
    \link{\code{heatmap.2}}
}
\examples{
    data(PromoterCount)
    heatmap.3(PromoterCount[1:100,])
}
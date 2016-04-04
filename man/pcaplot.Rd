\name{pcaplot}
\alias{pcaplot}
\title{2d and 3d sample relationship plot}
\description{
    Generate 2d or 3d sample relationship plot based on principal component analysis, multidimensional scaling, etc. 
}
\usage{
pcaplot(x, subset = NULL, cv.Th = 0.1, var.Th = 0, mean.Th = 0, standardize = TRUE, method = c("cluster", "mds", "pca"), dimension = c(1, 2, 3), color = "black", princurve = FALSE, lwd = 1, starts = NULL, col.curve = "red", text = TRUE, main = NULL, psi = 4, type = "p", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  A matrix of numeric values
  }
  \item{subset}{
  A numeric value indicating the number of genes that is randomly selected for pca analysis. Default to NULL, which means no subsetting procedure will be done
  }
  \item{cv.Th}{
  If subset = NULL, a numeric value indicating threshold of coeffcient of variation in selecting genes
  }
  \item{var.Th}{
  A numeric value indicating threshold of variation in selecting genes. This is only used when subset and cv.Th are both set to NULL
  }
  \item{mean.Th}{
  Similar to var.Th, a numeric value indicating threshold of mean value in selecting genes. This is only used when subset and cv.Th are both set to NULL
  }
  \item{standardize}{
  Whether to standardize samples so that each sample has mean 0 and variance 1. Default to TRUE
  }
  \item{method}{
  One of "cluster", "mds", "pca". Please refer to details section
  }
  \item{dimension}{
  Numeric vector indicating the number of dimensions you would like to generate the figure
  }
  \item{color}{
  Color for points when 'mds' or 'pca' is chosen as method
  }
  \item{princurve}{
  Logical value indicating whether to generate a principal curve. Please refer to details
  }
  \item{lwd}{
  The line width for principal curve 
  }
  \item{starts}{
  Providing the starting point for principal curve. Please refer to details
  }
  \item{col.curve}{
  The color of principal curve
  }
  \item{text}{
  Logical value indicating whether text is added as label to the figure
  }
  \item{main}{
  Main title for the figure
  }
  \item{psi}{
  Integer value indicating point size
  }
  \item{type}{
  For the default method, a single character indicating the type of item to plot. Supported types are: 'p' for points, 's' for spheres, 'l' for lines, 'h' for line segments from z = 0, and 'n' for nothing.
  }
  \item{\dots}{
  Further arguments will be ignored
  }
}
\details{
    If method = 'cluster', \code{hclust} is used; if method = 'mds', \code{cmdscale} is used; if method = 'pca', \code{prcomp} is used. 

    If princurve is set to TURE, then fits a principal curve which describes a smooth curve that passes through the middle of the data x in an orthogonal sense. This curve is a nonparametric generalization of a linear principal component. For details of principal curve, please refer to \code{\link{principal.curve}}. When princurve is set to TRUE, you need to provide a starting point for principal curve as \code{starts} argument. \code{starts} is basically a logical vector of the same length as number of samples, that tells you which sample will be used as starting point.
}
\value{
If method = 'cluster': a 'hclust' object
If method = 'mds' or 'pca', a data.frame containing user specified number of principal components. 
}
\references{
    https://github.com/Liuy12/MBDDiff
}
\author{
    Yuanhang Liu
}
\seealso{\code{\link{cmdscale}}, \code{\link{hclust}}, \code{\link{prcomp}}, \code{\link{principal.curve}}}
\examples{
    \dontrun{
        data(PromoterCount)
        Condition <- c(rep('C1', 3), rep('C2', 3))
        TestStat <- MBDDiff(Promoter, Background, Condition)
        MBD <- TestStat[[1]]
        Norm_count <- counts(MBD, normalized = TRUE)
        pcaplot(Norm_count, cv.Th = 0.1, method = 'pca', dimension = c(1,2,3), princurve = TRUE, starts = c(1,1,1,0,0,0)) 
        }
}
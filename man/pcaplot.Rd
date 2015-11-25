\name{pcaplot}
\alias{pcaplot}
\title{2d and 3d sample relationship plot}
\description{
    Generate 2d or 3d sample relationship plot based on principal component analysis, multidimensional scaling, etc. 
}
\usage{
pcaplot(x, subset = NULL, cv.Th = 0.1, var.Th = 0, mean.Th = 0, standardize = TRUE, method = c("cluster", "mds", "pca"), dimension = c(1, 2, 3), color = "black", princurve = F, lwd = 1, normals = NULL, col.curve = "red", text, main = NULL, psi = 4, type = "p", interactive = F, ...)
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
  A numeric value indicating threshold of coeffcient of variation 
  }
  \item{var.Th}{
%%     ~~Describe \code{var.Th} here~~
}
  \item{mean.Th}{
%%     ~~Describe \code{mean.Th} here~~
}
  \item{standardize}{
%%     ~~Describe \code{standardize} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{dimension}{
%%     ~~Describe \code{dimension} here~~
}
  \item{color}{
%%     ~~Describe \code{color} here~~
}
  \item{princurve}{
%%     ~~Describe \code{princurve} here~~
}
  \item{lwd}{
%%     ~~Describe \code{lwd} here~~
}
  \item{normals}{
%%     ~~Describe \code{normals} here~~
}
  \item{col.curve}{
%%     ~~Describe \code{col.curve} here~~
}
  \item{text}{
%%     ~~Describe \code{text} here~~
}
  \item{main}{
%%     ~~Describe \code{main} here~~
}
  \item{psi}{
%%     ~~Describe \code{psi} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
}
  \item{interactive}{

  }
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
    
}
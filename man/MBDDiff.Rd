\name{MBDDiff}
\alias{MBDDiff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MBDDiff(promoter, background, conditions, method = "pooled", sharingMode = "maximum", fitType = "local", pvals_only = FALSE, paraMethod = "NP")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{promoter}{
%%     ~~Describe \code{promoter} here~~
}
  \item{background}{
%%     ~~Describe \code{background} here~~
}
  \item{conditions}{
%%     ~~Describe \code{conditions} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{sharingMode}{
%%     ~~Describe \code{sharingMode} here~~
}
  \item{fitType}{
%%     ~~Describe \code{fitType} here~~
}
  \item{pvals_only}{
%%     ~~Describe \code{pvals_only} here~~
}
  \item{paraMethod}{
%%     ~~Describe \code{paraMethod} here~~
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (promoter, background, conditions, method = "pooled", 
    sharingMode = "maximum", fitType = "local", pvals_only = FALSE, 
    paraMethod = "NP") 
{
    XBSeq(promoter, background, conditions, method, sharingMode, 
        fitType, pvals_only, paraMethod)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\name{ConstructExRegions}
\alias{ConstructExRegions}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ConstructExRegions(organism, promoter_anno, CpG = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{organism}{
%%     ~~Describe \code{organism} here~~
}
  \item{promoter_anno}{
%%     ~~Describe \code{promoter_anno} here~~
}
  \item{CpG}{
%%     ~~Describe \code{CpG} here~~
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
function (organism, promoter_anno, CpG = T) 
{
    if (CpG) {
        ucsc_db <- src_mysql(organism, "genome-mysql.cse.ucsc.edu", 
            user = "genome")
        CpGisland <- as.data.frame(tbl(ucsc_db, sql("SELECT chrom, chromStart, chromEnd FROM cpgIslandExtUnmasked")))
        colnames(promoter_anno) <- colnames(CpGisland) <- c("chrom", 
            "start", "end")
        bind_rows(promoter_anno, as.data.frame(CpGisland))
    }
    else promoter_anno
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

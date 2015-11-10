\name{FilterRegions}
\alias{FilterRegions}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
FilterRegions(bed, fasta, Exbed)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bed}{
%%     ~~Describe \code{bed} here~~
}
  \item{fasta}{
%%     ~~Describe \code{fasta} here~~
}
  \item{Exbed}{
%%     ~~Describe \code{Exbed} here~~
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
function (bed, fasta, Exbed) 
{
    registerDoMC(40)
    t1 <- Sys.time()
    Exclude_index <- ExcludeIntersection(bed, Exbed)
    t2 <- Sys.time()
    t2 - t1
    letter_freq <- CountFreqency(fasta)
    Exclude_index1 <- which(letter_freq$ATCG != 1)
    Exclude_index_all <- Reduce(dplyr::union, list(Exclude_index, 
        Exclude_index2, Exclude_index3))
    return(bind_cols(bed[-Exclude_index_all, ], letter_freq$CG[-Exclude_index_all]))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

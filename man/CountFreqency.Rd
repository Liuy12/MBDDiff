\name{CountFreqency}
\alias{CountFreqency}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CountFreqency(fasta, CG = T, ATCG = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fasta}{
%%     ~~Describe \code{fasta} here~~
}
  \item{CG}{
%%     ~~Describe \code{CG} here~~
}
  \item{ATCG}{
%%     ~~Describe \code{ATCG} here~~
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
function (fasta, CG = T, ATCG = T) 
{
    fasta_file <- readDNAStringSet(fasta, use.names = T)
    Freqmat <- data.frame()
    if (CG) 
        Freqmat$CG <- letterFrequency(fasta_file, letters = "CG", 
            as.prob = T)
    if (ATCG) 
        Freqmat$ATCG <- letterFrequency(fasta_file, letters = "ATCG", 
            as.prob = T)
    return(Freqmat)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

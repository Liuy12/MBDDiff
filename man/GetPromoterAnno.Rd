\name{GetPromoterAnno}
\alias{GetPromoterAnno}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GetPromoterAnno(organism, save = F, Dir = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{organism}{
%%     ~~Describe \code{organism} here~~
}
  \item{save}{
%%     ~~Describe \code{save} here~~
}
  \item{Dir}{
%%     ~~Describe \code{Dir} here~~
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
function (organism, save = F, Dir = NULL) 
{
    ucsc_db <- src_mysql(organism, "genome-mysql.cse.ucsc.edu", 
        user = "genome")
    Ref_Flat <- tbl(ucsc_db, sql("SELECT geneName, chrom, strand, txStart, txEnd FROM refFlat"))
    Promoter_Anno <- GetPromoters(as.data.frame(Ref_Flat))
    if (save) {
        if (is.null(Dir)) 
            stop("please provide a path where annotation file will be saved")
        write.table(Promoter_Anno, paste(Dir, "/Promoter_Anno.bed", 
            sep = ""), quote = F, sep = "\t", row.names = F)
    }
    return(Promoter_Anno)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

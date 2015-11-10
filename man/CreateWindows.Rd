\name{CreateWindows}
\alias{CreateWindows}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CreateWindows(chromlen, binsize)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chromlen}{
%%     ~~Describe \code{chromlen} here~~
}
  \item{binsize}{
%%     ~~Describe \code{binsize} here~~
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
function (chromlen, binsize) 
{
    Wholegenome_bin <- foreach(i = 1:length(chromlen), .options.multicore = list(preschedule = FALSE), 
        .combine = bind_rows, .inorder = TRUE, .verbose = TRUE, 
        .errorhandling = "stop", .multicombine = TRUE) \%dopar\% 
        {
            cat(i, "\t")
            bin <- seq(0, chromlen[i], by = binsize)
            temp <- data.frame(chrom = rep(names(chromlen)[i], 
                length(bin)), chromStart = bin, chromEnd = c(bin[2:length(bin)], 
                chromlen[i]), name = paste(names(chromlen)[i], 
                "_", 1:length(bin), sep = ""))
            temp
        }
    return(Wholegenome_bin)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

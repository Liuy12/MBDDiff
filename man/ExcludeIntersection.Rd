\name{ExcludeIntersection}
\alias{ExcludeIntersection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ExcludeIntersection(data1, data2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data1}{
%%     ~~Describe \code{data1} here~~
}
  \item{data2}{
%%     ~~Describe \code{data2} here~~
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
function (data1, data2) 
{
    foreach(i = 1:nrow(data2), .options.multicore = list(preschedule = FALSE), 
        .combine = c, .inorder = TRUE, .verbose = TRUE, .errorhandling = "stop", 
        .multicombine = TRUE) \%dopar\% {
        cat(i, "\t")
        index <- which(data1[, 1] == data2[i, 1])
        index1 <- which(data1[index, 2] > data2[i, 2] & data1[index, 
            3] < data2[i, 3])
        index2 <- which(data1[index, 2] < data2[i, 2] & data1[index, 
            3] > data2[i, 2])
        index3 <- which(data1[index, 2] < data2[i, 3] & data1[index, 
            3] > data2[i, 3])
        index4 <- index[c(index1, index2, index3)]
        if (length(index4)) 
            index4
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

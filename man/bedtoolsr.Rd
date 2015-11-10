\name{bedtoolsr}
\alias{bedtoolsr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bedtoolsr(bamdir, bed)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bamdir}{
%%     ~~Describe \code{bamdir} here~~
}
  \item{bed}{
%%     ~~Describe \code{bed} here~~
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
function (bamdir, bed) 
{
    bamfiles <- grep(".bam", dir(bamdir), value = T)
    if (length(bamfiles)) 
        stop("There are no bam files in the directory")
    else {
        bamfiles <- paste(bamdir, "/", bamfiles, sep = "")
        temp1 <- c()
        for (i in length(bamfiles)) {
            command <- paste("coverageBed -abam", bamfiles[i], 
                "-b", bed, sep = " ")
            cat("processing: ", bamfiles[i], sep = "\t")
            temp <- fread(input = command, data.table = F)
            temp <- arrange(temp, V1, V2, V3)
            temp1 <- bind_cols(temp1, temp$V6)
        }
        rownames(temp1) <- temp$V4
        return(list(count = temp1, gaplength = temp[, (ncol(temp) - 
            1)]))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

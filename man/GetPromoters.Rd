\name{GetPromoters}
\alias{GetPromoters}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GetPromoters(uctable, upstream = 2000, downstream = 2000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{uctable}{
%%     ~~Describe \code{uctable} here~~
}
  \item{upstream}{
%%     ~~Describe \code{upstream} here~~
}
  \item{downstream}{
%%     ~~Describe \code{downstream} here~~
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
function (uctable, upstream = 2000, downstream = 2000) 
{
    uctable <- uctable \%>\% filter(sapply(strsplit(chrom, "_"), 
        length) == 1)
    uniqid <- uctable \%>\% select(geneName) \%>\% distinct(geneName)
    uctable1 <- data.frame()
    for (i in 1:length(uniqid$geneName)) {
        cat(i, "\n")
        index <- grep(paste("^", uniqid$geneName[i], "$", sep = ""), 
            uctable$geneName)
        temp <- unique(uctable$chrom[index])
        temp1 <- unique(uctable$strand[index])
        for (j in 1:length(temp)) {
            for (k in 1:length(temp1)) {
                index1 <- which(uctable[index, 2] == temp[j] & 
                  uctable[index, 3] == temp1[k])
                if (length(index1)) {
                  if (temp1[k] == "+") 
                    temp2 <- unique(uctable[index[index1], 4])
                  else temp2 <- unique(uctable[index[index1], 
                    5])
                  temp3 <- data.frame(Symbol = paste(uniqid$geneName[i], 
                    "#", 1:length(temp2), sep = ""), Chrom = rep(temp[j], 
                    length(temp2)), Strand = rep(temp1[k], length(temp2)), 
                    TSS = temp2, stringsAsFactors = F)
                  uctable1 <- rbind(uctable1, temp3)
                }
            }
        }
    }
    Promoter4k <- data.frame(chrom = uctable1$Chrom, chromStart = ifelse(uctable1$Strand == 
        "+", uctable1$TSS - upstream, uctable1$TSS - downstream), 
        chromEnd = ifelse(uctable1$Strand == "+", uctable1$TSS + 
            downstream, uctable1$TSS + upstream), name = uctable1$Symbol, 
        Strand = uctable1$Strand)
    return(Promoter4k)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

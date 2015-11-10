\name{MethyEnrich}
\alias{MethyEnrich}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MethyEnrich(bed, bam, fa)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bed}{
%%     ~~Describe \code{bed} here~~
}
  \item{bam}{
%%     ~~Describe \code{bam} here~~
}
  \item{fa}{
%%     ~~Describe \code{fa} here~~
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
function (bed, bam, fa) 
{
    whole_genome_bin_count <- bedtoolsr(bed, bam)
    whole_genome_bin_count_TPM <- CalculateTPM(whole_genome_bin_count[[1]], 
        whole_genome_bin_count[[2]])
    GCcon <- CountFreqency(fa, CG = T, ATCG = F)$CG
    temp <- quantile(whole_genome_bin_count_TPM)
    cuts <- cut(whole_genome_bin_count_TPM, breaks = unique(temp), 
        labels = c("0-25\%", "25-50\%", "50-75\%", "75-100\%"))
    temp1 <- bind_cols(temp, as.data.frame(cuts))
    colnames(temp1) <- c("RPKM", "CG", "cuts")
    ggplot() + geom_density(aes(x = CG, color = cuts), adjust = 2, 
        data = temp1) + scale_colour_manual(values = c("black", 
        "blue", "purple", "gray"))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

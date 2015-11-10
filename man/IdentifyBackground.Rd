\name{IdentifyBackground}
\alias{IdentifyBackground}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
IdentifyBackground(bed, promo_bed)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bed}{
%%     ~~Describe \code{bed} here~~
}
  \item{promo_bed}{
%%     ~~Describe \code{promo_bed} here~~
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
function (bed, promo_bed) 
{
    registerDoMC(80)
    t1 <- Sys.time()
    Background_region <- foreach(i = 1:nrow(promo_bed), .options.multicore = list(preschedule = FALSE), 
        .combine = bind_rows, .inorder = TRUE, .verbose = TRUE, 
        .errorhandling = "stop", .multicombine = TRUE) \%dopar\% 
        {
            cat("promoter", i, "\t")
            preset <- filter(bed, Chrom == promo_bed[i, 1])
            preset <- mutate(preset, proximity = abs(Start - 
                promo_bed[i, 7]))
            preset <- arrange(preset, proximity)
            k <- 0
            bgregion <- c()
            for (j in 1:nrow(preset)) {
                if (preset[j, 5] < 0.4) {
                  k <- k + 1
                  bgregion <- bind_rows(bgregion, preset[j, ])
                }
                if (k == 80) 
                  break
            }
            bgregion
        }
    t2 <- Sys.time()
    t2 - t1
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\name{Methydist}
\alias{Methydist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Methydist(promoter, background, gaplength_p, gaplength_b)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{promoter}{
%%     ~~Describe \code{promoter} here~~
}
  \item{background}{
%%     ~~Describe \code{background} here~~
}
  \item{gaplength_p}{
%%     ~~Describe \code{gaplength_p} here~~
}
  \item{gaplength_b}{
%%     ~~Describe \code{gaplength_b} here~~
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
function (promoter, background, gaplength_p, gaplength_b) 
{
    Promoter_TPM <- CalculateTPM(promoter, gaplength)
    Background_TPM <- CalculateTPM(background, gaplength_b)
    temp <- data.frame(TPM <- Promoter_TPM, group <- rep("Promoter", 
        length(Promoter_TPM)))
    colnames(temp) <- c("TPM", "group")
    temp1 <- data.frame(TPM <- Background_TPM, group <- rep("Background", 
        length(Background_TPM)))
    colnames(temp1) <- c("TPM", "group")
    temp2 <- bind_rows(temp, temp1)
    ggplot() + geom_histogram(aes(x = count, y = ..density.., 
        fill = group, alpha = group), data = temp2, binwidth = 1, 
        position = "identity") + xlim(c(0, 150)) + scale_fill_manual(values = c("blue", 
        "red")) + scale_alpha_manual(values = c(0.6, 0.6))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

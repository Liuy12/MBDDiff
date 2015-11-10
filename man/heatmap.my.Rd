\name{heatmap.my}
\alias{heatmap.my}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
heatmap.my(Exprs, sel = F, thres_mean, thres_var, numbreaks = 100, col = c("blue", "white", "red"), breakratio = c(2, 1, 2), colsidebar, Colv = F, Rowv = T, scale = "row", labRow = F, labCol = F, dendrogram = "row")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Exprs}{
%%     ~~Describe \code{Exprs} here~~
}
  \item{sel}{
%%     ~~Describe \code{sel} here~~
}
  \item{thres_mean}{
%%     ~~Describe \code{thres_mean} here~~
}
  \item{thres_var}{
%%     ~~Describe \code{thres_var} here~~
}
  \item{numbreaks}{
%%     ~~Describe \code{numbreaks} here~~
}
  \item{col}{
%%     ~~Describe \code{col} here~~
}
  \item{breakratio}{
%%     ~~Describe \code{breakratio} here~~
}
  \item{colsidebar}{
%%     ~~Describe \code{colsidebar} here~~
}
  \item{Colv}{
%%     ~~Describe \code{Colv} here~~
}
  \item{Rowv}{
%%     ~~Describe \code{Rowv} here~~
}
  \item{scale}{
%%     ~~Describe \code{scale} here~~
}
  \item{labRow}{
%%     ~~Describe \code{labRow} here~~
}
  \item{labCol}{
%%     ~~Describe \code{labCol} here~~
}
  \item{dendrogram}{
%%     ~~Describe \code{dendrogram} here~~
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
function (Exprs, sel = F, thres_mean, thres_var, numbreaks = 100, 
    col = c("blue", "white", "red"), breakratio = c(2, 1, 2), 
    colsidebar, Colv = F, Rowv = T, scale = "row", labRow = F, 
    labCol = F, dendrogram = "row") 
{
    suppressPackageStartupMessages(invisible(require("gplots", 
        quietly = TRUE)))
    if (labRow) 
        labRow <- rownames(Exprs)
    if (labCol) 
        labCol <- colnames(Exprs)
    if (sel) {
        gene_mean <- apply(Exprs, 1, mean)
        gene_var <- apply(Exprs, 1, var)
        Exprs <- Exprs[gene_mean > thres_mean & gene_var > thres_var, 
            ]
    }
    if (scale == "row") 
        Exprs_scale <- t(scale(t(Exprs)))
    else Exprs_scale <- Exprs
    if (missing(colsidebar)) {
        lmat <- rbind(c(0, 4), c(0, 3), c(2, 1))
        lwid <- c(1, 4)
        lhei <- c(1, 0.1, 4)
        if (class(Colv) == "dendrogram") {
            lhei <- c(1, 1, 4)
            dendrogram <- "both"
        }
    }
    else {
        if (class(Colv) == "dendrogram") {
            lmat <- rbind(c(0, 5), c(0, 4), c(3, 2), c(0, 1))
            lwid <- c(1, 4)
            lhei <- c(1, 1, 4, 0.25)
            dendrogram <- "both"
        }
        else {
            if (Colv) {
                lmat <- rbind(c(0, 5), c(0, 4), c(3, 2), c(0, 
                  1))
                lwid <- c(1, 4)
                lhei <- c(1, 1, 4, 0.25)
                dendrogram <- "both"
            }
            lmat <- rbind(c(0, 5), c(0, 1), c(3, 2), c(0, 4))
            lwid <- c(1, 4)
            lhei <- c(1, 0.25, 4, 0.1)
        }
    }
    rg <- quantile(Exprs_scale, na.rm = T)
    rg_diff <- rg[4] - rg[2]
    rg_max <- max(abs(rg))
    Exprs_sd <- sd(Exprs_scale)
    Exprs_mean <- mean(Exprs_scale)
    if (rg_max > max(abs(c(Exprs_mean + 3 * Exprs_sd, Exprs_mean - 
        3 * Exprs_sd)))) {
        rg_iqr <- max(abs(c(rg[2], rg[4])))
        bp <- c((breakratio[1]/sum(breakratio)) * rg_diff - rg_iqr, 
            rg_iqr - (breakratio[3]/sum(breakratio)) * rg_diff)
        bk <- c(seq(-rg_max, -rg_iqr, length = numbreaks), seq(-rg_iqr, 
            bp[1], length = numbreaks), seq(bp[1], bp[2], length = numbreaks), 
            seq(bp[2], rg_iqr, length = numbreaks), seq(rg_iqr, 
                rg_max, length = numbreaks))
        hmcols <- colorRampPalette(col)(length(bk) - 1)
    }
    else {
        rg <- range(Exprs_scale, na.rm = T)
        bp <- c((breakratio[1]/sum(breakratio)) * diff(rg) - 
            rg_max, rg_max - (breakratio[3]/sum(breakratio)) * 
            diff(rg))
        bk <- c(seq(-rg_max, bp[1], length = numbreaks), seq(bp[1], 
            bp[2], length = numbreaks), seq(bp[2], rg_max, length = numbreaks))
        hmcols <- colorRampPalette(col)(length(bk) - 1)
    }
    heatmap.2(Exprs, Colv = Colv, Rowv = Rowv, dendrogram = dendrogram, 
        trace = "none", scale = scale, density.info = "none", 
        lmat = lmat, lwid = lwid, lhei = lhei, labRow = labRow, 
        labCol = labCol, col = hmcols, breaks = bk, ColSideColors = colsidebar)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

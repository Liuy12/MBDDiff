\name{pcaplot}
\alias{pcaplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pcaplot(x, subset = NULL, cv.Th = 0.1, var.Th = 0, mean.Th = 0, standardize = TRUE, method = c("cluster", "mds", "pca"), dimension = c(1, 2, 3), color = "black", princurve = F, lwd = 1, normals = NULL, col.curve = "red", text, main = NULL, psi = 4, type = "p", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{subset}{
%%     ~~Describe \code{subset} here~~
}
  \item{cv.Th}{
%%     ~~Describe \code{cv.Th} here~~
}
  \item{var.Th}{
%%     ~~Describe \code{var.Th} here~~
}
  \item{mean.Th}{
%%     ~~Describe \code{mean.Th} here~~
}
  \item{standardize}{
%%     ~~Describe \code{standardize} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{dimension}{
%%     ~~Describe \code{dimension} here~~
}
  \item{color}{
%%     ~~Describe \code{color} here~~
}
  \item{princurve}{
%%     ~~Describe \code{princurve} here~~
}
  \item{lwd}{
%%     ~~Describe \code{lwd} here~~
}
  \item{normals}{
%%     ~~Describe \code{normals} here~~
}
  \item{col.curve}{
%%     ~~Describe \code{col.curve} here~~
}
  \item{text}{
%%     ~~Describe \code{text} here~~
}
  \item{main}{
%%     ~~Describe \code{main} here~~
}
  \item{psi}{
%%     ~~Describe \code{psi} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (x, subset = NULL, cv.Th = 0.1, var.Th = 0, mean.Th = 0, 
    standardize = TRUE, method = c("cluster", "mds", "pca"), 
    dimension = c(1, 2, 3), color = "black", princurve = F, lwd = 1, 
    normals = NULL, col.curve = "red", text, main = NULL, psi = 4, 
    type = "p", ...) 
{
    if (is.matrix(x)) {
        dataMatrix <- x
    }
    else {
        stop("The class of \"x\" should be matrix!")
    }
    if (is.null(subset)) {
        if (!is.null(cv.Th)) {
            cv.gene <- apply(dataMatrix, 1, function(x) sd(x)/mean(x))
            dataMatrix <- dataMatrix[cv.gene > cv.Th, ]
            subset <- 1:nrow(dataMatrix)
        }
        else {
            if (!is.null(var.Th)) {
                var.gene <- apply(dataMatrix, 1, function(x) var(x))
                dataMatrix <- dataMatrix[var.gene > var.Th, ]
                subset <- 1:nrow(dataMatrix)
            }
            if (!is.null(mean.Th)) {
                mean.gene <- apply(dataMatrix, 1, function(x) mean(x))
                dataMatrix <- dataMatrix[mean.gene > mean.Th, 
                  ]
                subset <- 1:nrow(dataMatrix)
            }
        }
        if (is.null(main)) 
            main <- paste("Sample relations based on", length(subset), 
                "genes")
    }
    else {
        if (length(subset) == 1 && is.numeric(subset)) {
            subset <- sample(1:nrow(dataMatrix), min(subset, 
                nrow(dataMatrix)))
        }
        if (is.null(main)) 
            main <- paste("Sample relations based on", length(subset), 
                "selected genes")
    }
    if (standardize) 
        dataMatrix <- scale(dataMatrix)
    method <- match.arg(method)
    if (method == "cluster") {
        dd <- dist(t(dataMatrix[subset, ]))
        hc = hclust(dd, "ave")
        plot(hc, xlab = "Sample", main = main, ...)
        attr(hc, "geneNum") <- length(subset)
        return(invisible(hc))
    }
    if (method == "mds") {
        dd <- dist(t(dataMatrix[subset, ]))
        mds.result <- cmdscale(dd, k = max(dimension), eig = TRUE)
        ppoints <- mds.result$points
        eig <- mds.result$eig
        percent <- round(eig/sum(eig) * 100, 1)
        if (is.null(color)) {
            color <- 1
        }
        else {
            if (!is.numeric(color)) {
                allColor <- colors()
                if (!all(is.element(color, allColor))) {
                  color <- as.numeric(factor(color, levels = unique(color)))
                }
            }
        }
        require("rgl")
        plot3d(ppoints[, dimension[1]], ppoints[, dimension[2]], 
            ppoints[, dimension[3]], xlab = paste("Principal Component ", 
                dimension[1], " (", percent[dimension[1]], "\%)", 
                sep = ""), ylab = paste("Principal Component ", 
                dimension[2], " (", percent[dimension[2]], "\%)", 
                sep = ""), zlab = paste("Principal Component ", 
                dimension[3], " (", percent[dimension[3]], "\%)", 
                sep = ""), main = main, size = psi, col = color, 
            type = type)
        if (princurve) {
            start <- aggregate(ppoints[, 1:3], by = list(rank(!normals)), 
                FUN = mean)
            start <- as.matrix(start[, -1])
            fit <- principal.curve(ppoints[, 1:3], start = start, 
                plot.true = F)
            plot3d(fit$s[fit$tag, ], type = "l", add = T, col = col.curve, 
                lwd = lwd)
        }
        if (text) 
            text3d(ppoints[, dimension[1]], ppoints[, dimension[2]], 
                ppoints[, dimension[3]], col = color, texts = colnames(dataMatrix), 
                cex = 1)
        attr(ppoints, "geneNum") <- length(subset)
        pp <- ppoints
        return(pp)
    }
    if (method == "pca") {
        pca.result <- prcomp(t(dataMatrix))
        ppoints <- pca.result$x[, 1:max(dimension)]
        percent <- round((pca.result$sdev^2/sum(pca.result$sdev^2)) * 
            100, 1)
        if (is.null(color)) {
            color <- 1
        }
        else {
        }
        require("rgl")
        plot3d(ppoints[, dimension[1]], ppoints[, dimension[2]], 
            ppoints[, dimension[3]], xlab = paste("Principal Component ", 
                dimension[1], " (", percent[dimension[1]], "\%)", 
                sep = ""), ylab = paste("Principal Component ", 
                dimension[2], " (", percent[dimension[2]], "\%)", 
                sep = ""), zlab = paste("Principal Component ", 
                dimension[3], " (", percent[dimension[3]], "\%)", 
                sep = ""), main = main, size = psi, col = color, 
            type = type, pch = 19)
        if (princurve) {
            start <- aggregate(ppoints[, 1:3], by = list(rank(!normals)), 
                FUN = mean)
            start <- as.matrix(start[, -1])
            fit <- principal.curve(ppoints[, 1:3], start = start, 
                plot.true = F)
            plot3d(fit$s[fit$tag, ], type = "l", add = T, col = col.curve, 
                lwd = lwd)
        }
        if (text) {
            text3d(ppoints[, dimension[1]], ppoints[, dimension[2]], 
                ppoints[, dimension[3]], col = color, texts = colnames(dataMatrix), 
                cex = 1)
        }
        pp <- ppoints
        return(pp)
    }
    else {
        stop("the method has to be one of the three,'cluster','mds' or 'pca'")
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

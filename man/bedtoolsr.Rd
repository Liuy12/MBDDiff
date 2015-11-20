\name{bedtoolsr}
\alias{bedtoolsr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Counting reads mapped to a specified region using bedtools
}
\description{
   Given bam file or directories of bam files, counting reads mapped to genomic intervals specified by user. 
}
\usage{
   bedtoolsr(bam = NULL, bamdir = NULL, bed)
}
\arguments{
  \item{bam}{
   bam file, including absolute path of the file. This is usually the output of alignment programs like, bwa, bowtie, etc. 
}
  \item{bamdir}{
   Alternatively, you can also provide the path of a list of bam files to process all at once.
}
  \item{bed}{
   beb files, which include genomic intervels you want to measure.
}
}
\details{
   
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

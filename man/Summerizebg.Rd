\name{Summerizebg}
\alias{Summerizebg}
\title{
	Summerizing background noise for each promoter region
}
\description{
Summerizing background noise for each promoter region based on methylation levels}
\usage{
Summerizebg(dataMat_bg, gaplength)
}
\arguments{
  \item{dataMat_bg}{
  A matrix or data.frame containing background mapped reads for selected bins
  }
  \item{gaplength}{
  Gap length of bins
  }
}
\details{
	This function is generally called after \code{bedtoolsr} for counting reads mapped to selected bins for background annotation constructed by \code{IdentifyBackground}
}
\value{
A data.frame of the same dimension as promoter mapped reads
}
\references{
https://github.com/Liuy12/MBDDiff
}
\author{
	Yuanhang Liu
}
\seealso{
	\code{\link{IdentifyBackground}}, \code{\link{Summerizebg}}
}
\examples{
	\dontrun{
		Background_counts <- bedtoolsr(bamdir = '/path/to/bam/files', bed = '/path/to/bg_bed/file')
		Background_counts_summe <- Summerizebg(Background_counts[[1]], Background_counts[[2]])
		}
}
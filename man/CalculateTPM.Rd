\name{CalculateTPM}
\alias{CalculateTPM}
\title{Calculate transcript per million (TPM) for count matrix}
\description{
	Calculate TPM based on user-specifed count matrix and transcript length. This can be generalized to any certain genomic interval instead of limiting to transcript.
}
\usage{
CalculateTPM(dataMat, gaplength)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataMat}{
  a matrix or data.frame of count datasets
  }
  \item{gaplength}{
  a integer vector that contains length information for each row of dataMat.
  }
}
\value{
	A matrix or data.frame of TPM normalized values
}
\references{
	https://github.com/Liuy12/MBDDiff
}
\author{
	Yuanhang Liu
}
\seealso{
\code{\link{Summerizebg}}
}
\examples{
data(PromoterCount)
PromoterCount_TPM <- CalculateTPM(PromoterCount, Promoterlength)
}
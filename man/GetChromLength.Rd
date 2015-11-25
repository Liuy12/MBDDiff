\name{GetChromLength}
\alias{GetChromLength}
\title{Retrive chromsome length information from UCSC database}
\description{
Retrive chromsome length information from UCSC MySQL database 
}
\usage{GetChromLength(organism)}
\arguments{
	\item{organism}{
	The accession code of organisms to UCSC database.
	}
}
\details{
	To get a list of availabe organisms in UCSC database, see \code{\link{ShowAvailableOrg}}
}
\value{
	A data.frame that contains chromesome information:
	\item{chrom}{
	chromsome number
	}
	\item{size}{
	chromsome length
	}
}
\references{
	https://github.com/Liuy12/MBDDiff
}
\examples{
GetChromLength('hg19')
}
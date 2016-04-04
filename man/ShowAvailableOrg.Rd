\name{ShowAvailableOrg}
\alias{ShowAvailableOrg}
\alias{Organisms}
\title{
  Show a list of available organisms from UCSC database
}
\description{
    Show a list of available organisms from UCSC database
}
\usage{
  ShowAvailableOrg()
}
\value{
A data.frame of two columns:
\item{Organism name}{Description of organism name}
\item{Accession code}{Accession code to be used for each organism}
}
\references{
  https://github.com/Liuy12/MBDDiff
}
\author{Yuanhang Liu}
\seealso{\code{\link{GetPromoterAnno}}, \code{\link{GetChromLength}}, \code{\link{IdentifyBackground}}}
\examples{
ShowAvailableOrg()
}
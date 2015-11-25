\name{MethyEnrich}
\alias{MethyEnrich}
\title{
  Quality control of MBDcap-seq datasets
}
\description{
  Generate graph of enrichment of GC content grouped by methylation levels}
\usage{
MethyEnrich(bed_count, fa, interactive = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bed_count}{
  A vector containing number of reads mapped to each bin. See vignettes for details
  }
  \item{fa}{
  Absolute path of sequence file annotation for your interested organism
  }
  \item{interactive}{
  Logical variable indicating whether you want to generate a interactive version of the same graph
  }
}
\details{
  In order to use this function, you will firstly need to create a sliding window across the genome and measure number of reads mapped to each window. Please refer to vignettes for details
}
\references{
  https://github.com/Liuy12/MBDDiff
}
\author{
  Yuanhang Liu
}
\seealso{
  \code{\link{bedtoolsr}}, \code{\link{CreateWindows}}
}
\examples{
  \dontrun{
    whole_genome_bin_count <- bedtoolsr(bed, bam)
    whole_genome_bin_count_TPM <- CalculateTPM(whole_genome_bin_count[[1]], whole_genome_bin_count[[2]])
    MethyEnrich(whole_genome_bin_count_TPM, fa)
    }
}
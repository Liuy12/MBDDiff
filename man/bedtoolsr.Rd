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
   This function is basically a wrapper for bedtools, which means you still need to have bedtools installed on your computer and saved in your system path. 
}
\value{
A list of two elements:
\item{count}{
	a vector or matrix of values that indicates the number of reads that mapped to user-specified intervals
}
\item{gaplength}{
	a integer vector that indicates the gap length of user-specified intervals
}
}
\references{
https://github.com/Liuy12/MBDDiff
}
\author{
Yuanhang Liu
}
\examples{
\dontrun{
	# Process several bam files in a fold 
	Promoter_counts <- bedtoolsr(bamdir = '/path/to/bam/files', bed = '/path/to/bed/file')
	}
}
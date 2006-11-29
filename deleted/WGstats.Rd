\name{WGstats}
\alias{WGstats}
\title{Statistics for an object of class 'WGassociation'}

\description{
  Prints, for each SNP, a detailed table as the output of \code{\link{association}}}

\usage{
WGstats(object, ...)
}

\arguments{
  \item{object}{an object of class \code{\link{WGassociation}}}
}
\details{
  The function computes, for each chromosome (or gene), the number of SNPs,
  the percentage of SNPs with genotype error (e.g SNPs where the genotype rate is
  less than those indicated in the argument 'genotypingRate' from \code{\link{WGassociation}}
  function), the percentage of Monomorphic SNPs, and 
  the number of statistically significant associations at level given in 
  the 'pSig' argument.
}
\value{
  a data frame where rows are chromosomes (or genes) and the following columns:
  number of SNPs, percentage of SNPs with genotype rate less than 'genotypingRate',
  the percentage of Monomorphic SNPs, and the number of statistically significant 
  associations at level indicated in 'pSig'.
}

\seealso{\code{\link{WGassociation}}}

\examples{
data(HapMap)
WGstats(resHapMap)
}

\keyword{utilities}

\name{WGstats}
\alias{WGstats}
\title{Statistics for an object of class 'WGassociation'}

\description{
  Gives a sort summary for each chromosome (or gene) for an object of class 'WGassociation'}

\usage{
WGstats(x, pSig = 1e-06)
}

\arguments{
  \item{x}{an object of class \code{\link{WGassociation}}}
  \item{pSig}{level of significance. The default is 0.000001. See details}
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

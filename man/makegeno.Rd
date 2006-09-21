\name{make.geno}
\alias{make.geno}

\title{Create a group of locus objects from some SNPs, assign to 'model.matrix' class.
}

\description{
  This function prepares the CRITICAL element corresponding to matrix of genotypes necessary to be included in 'haplo.glm' function.
}

\usage{
make.geno(data, SNPs.sel, sep="/")
}

\arguments{
  \item{data}{data frame containing the the SNPs that will be used to estimate the haplotypes.}
  \item{SNPs.sel}{ a vector indicating the names of SNPs that are used to estimate the haplotypes  }
  \item{sep}{character separator used to divide alleles. See \code{snp} function for details.}
}


\value{
  the same as 'setupGeno' function, from 'haplo.stats' library, returns
}


\seealso{\code{\link{snp}}}

\examples{
data(SNPs)
geno<-make.geno(SNPs,c("snp10001","snp10002","snp10003"),sep="")
}


\keyword{utilities}
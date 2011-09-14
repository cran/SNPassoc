\name{dscore}
\alias{dscore}
\alias{dscore.default}
\alias{dscore.setupSNP}

\title{Genetic risk allele score}
\description{
  Computes the exact distribution of a genetic risk allele score 
}
\usage{
dscore(x, \dots)

\method{dscore}{default}(x, \dots)

\method{dscore}{setupSNP}(x, \dots)

}

\arguments{
  \item{x}{a vector of probabilities representing the MAF of each SNP or an object of class setupSNP.}

\item{\dots}{further arguments to be passed to or from methods.}

}

\references{
Lucas G, et al. (2011). Submitted to American Journal of Human Genetics
}

\value{
  A vector with probabilities for each allele count 
}

\seealso{
  \code{\link{snp}, \link{setupSNP} }
}


\examples{

# example with 4 SNPs - the user gives the probabilities
MAFs <- c(0.1, 0.07, 0.2, 0.4)
dscore(MAFs)

# example with 4 SNPs - using setupSNP
data(SNPs)
myDat<-setupSNP(SNPs,6:10,sep="")
dscore(myDat)

}


\keyword{utilities}

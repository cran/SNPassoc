\name{tableHWE}
\alias{tableHWE}
\alias{print.tableHWE}

\title{Test for Hardy-Weinberg Equilibrium}
\description{
  Test the null hypothesis that Hardy-Weinberg equilibrium holds in cases, 
  controls and both populations.
}
\usage{
tableHWE(x, strata)
print.tableHWE(x, digits=4, sig=0.05, na="-", ...) 
}

\arguments{
  \item{x}{an object of class 'setupSNP'}
  \item{strata}{a factor variable for a stratified analysis}
  \item{digits}{number of significant digits, see 'print.default'}
  \item{sig}{significance level in the flag column}
  \item{na}{character string (or 'NULL') indicating 'NA' values in
          printed output, see 'print.default}
  \item{...}{further arguments passed to or from other methods}
}
\details{
  This function calculates the HWE test for those variables of class 'snp' in the object x of 
  class 'setupSNP'
}
\value{
  A matrix with p values for Hardy-Weinberg Equilibrium 
}

\author{This function is based on an R function which computes an exact SNP test of Hardy-Weinberg Equilibrium 
written by Wigginton JE, Cutler DJ and Abecasis GR available at \url{http://www.sph.umich.edu/csg/abecasis/Exact/r_instruct.html}} 

\seealso{ \code{\link{setupSNP}}}

\references{Wigginton JE, Cutler DJ and Abecasis GR (2005). A note on exact tests of Hardy-Weinberg equilibrium. Am J Hum Genet 76:887-93 

}

\examples{
data(SNPs)
ans<-setupSNP(SNPs,6:40,sep="")
res<-tableHWE(ans)
print(res)
#change the significance level showed in the flag column
print(res,sig=0.001)


#stratified analysis
res<-tableHWE(ans,ans$sex)
print(res)


}
\keyword{utilities}


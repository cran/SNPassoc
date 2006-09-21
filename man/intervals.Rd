\name{intervals}
\alias{intervals}
\alias{intervals.haplo.glm}
\alias{print.intervals}
\alias{summary.haplo.glm}

\title{Print ORs and 95\% confidence intervals for an object of class 'haplo.glm'}

\description{Print ORs and confidence intervals for an object of class 'haplo.glm'}

\usage{
intervals(o, level=.95, ...)
intervals.haplo.glm(o, level=.95, sign = 1, FUN=exp, ...)
print.intervals(x, len = 6, d = 2, exclude.intercept=TRUE, pval=TRUE, ...)
summary.haplo.glm(object, ...)
}

\arguments{
  \item{o}{object of class 'haplo.glm'}
  \item{x}{object of class 'haplo.glm'}
  \item{object}{object of class 'haplo.glm'}
  \item{level}{significance level. Default is 95 percent}
  \item{sign}{sign}
  \item{FUN}{}
  \item{len}{length}
  \item{d}{number of decimals}
  \item{exclude.intercept}{is the intercept excluded?}
  \item{pval}{is p value printed?}
  \item{...}{other arguments}

}

\examples{
# Not Run
# data(SNPs)
# tag.SNPs<-c("snp100019","snp10001","snp100029")
# geno<-make.geno(SNPs,tag.SNPs)

# mod<-haplo.glm(casco~geno,data=SNPs, 
#      family=binomial,
#	locus.label=tag.SNPs,
#	allele.lev=attributes(geno)$unique.alleles,
#	control = haplo.glm.control(haplo.freq.min=0.05))

# intervals(mod)
# summary(mod)

}
\keyword{utilities}
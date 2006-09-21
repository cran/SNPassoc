\name{WGassociation}
\alias{WGassociation}
\alias{print.WGassociation}
\alias{summary.WGassociation}

\title{ Whole genome association analysis}

\description{
  This function carries out a whole genome association analysis between the SNPs and
  a dependent variable (phenotype) under five different genetic models (inheritance patterns): 
  codominant, dominant, recessive, overdominant and log-additive. The phenotype may be quantitative 
  or categorical. In the second case (e.g. case-control studies) this variable must be of class 'factor' 
  with two levels.

}
\usage{
   WGassociation(formula, data, model = c("all"), quantitative = is.quantitative(formula, data),
                 genotypingRate = 80, level = 0.95)
}

\arguments{
  \item{formula}{ a symbolic description of the model to be fited (a formula object) without the SNP. 
                It might have either a continuous variable (quantitative traits) or a 
                factor variable (case-control studies) as the response on the left of the \code{~} 
                operator and terms with additional covariates on the right of the \code{~} operator may be 
                added to fit an adjusted model (e.g., \code{~}var1+var2+...+varN+SNP). See details}
  \item{data}{ a required dataframe of class 'setupSNP' containing the variables in the model and the SNPs}
  \item{model}{a character string specifying the type of genetic model (mode of inheritance) for the SNP. 
               This indicates how the genotypes should be collapsed. Possible
               values are "codominant", "dominant", "recessive", "overdominant", "log-additive" or "all". The default
               is "all" that fits the 5 possible genetic models. Only the first words are required, e.g "co", "do", etc.}
  \item{quantitative}{logical value indicating whether the phenotype (that which is in the
                      left of the operator ~ in 'formula' argument) is quantitative. The function 
                      'is.quantitative' returns FALSE when the phenotype
                      is a variable with two categories (i.e. indicating case-control status). Thus,
                      it is not a required argument but it may be modified by the user.}
  \item{genotypingRate}{ minimum percentage of genotype rate for a given SNP to be included in the analysis.
                  Default is 80\%.}
  \item{level}{signification level for confidence intervals. Defaul 95\%.}
}

\details{
  This function assesses the association between the response variable included in the left side in 
  the 'formula' and the SNPs included in the 'data' argument adjusted by those variables included 
  in the right side of the 'formula'. Different genetic models may be analyzed using 'model' argument.
}


\value{
  An object of class 'WGassociation'. 

  The functions 'print' and 'summary' are used to obtain and print a summary of the results. The 
  p values are saved in the attribute 'pvalues' as a matrix. They may be obtained using attr(,"pvalues") 
  (see examples). The first column indicates whether a problem with genotyping is present.

  The function 'plot' is used to obtain a plot of p values in the -log scale. See plot.WGassociation for
  further details
}

\seealso{ \code{\link{getSignificantSNPs}} \code{\link{association}} 
     \code{\link{setupSNP}} \code{\link{plot.WGassociation}}}

\examples{
data(SNPs)
datSNP<-setupSNP(SNPs,6:40,sep="")
ansAll<-WGassociation(protein~1,data=datSNP,model="all")


#only codominant and log-additive
ansCoAd<-WGassociation(protein~1,data=datSNP,model=c("co","log-add"))

#for printing p values
print(ansAll)
print(ansCoAd)

#for obtaining a matrix with the p palues
pvalAll<-attr(ansAll,"pvalues")
pvalCoAd<-attr(ansCoAd,"pvalues")

#for a big summary
summary(ansAll)
summary(ansCoAd)

#for plotting the p values
plot(ansAll)


#
# Whole genome analysis
#
data(HapMap)


# Next steps may be very time consuming. So they are not executed

#myDat<-setupSNP(HapMap, colSNPs=3:9809, sort = TRUE,
#   info=HapMap.SNPs.pos, sep="")
#resHapMap<-WGassociation(group~1, data=myDat, model="log")


# However, the results are saved in the object "resHapMap"
# to illustrate print, summary and plot functions
print(resHapMap)
summary(resHapMap)
plot(resHapMap)


}

\keyword{utilities}

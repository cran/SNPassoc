% --- SNPstat-internal.Rd ---

\name{SNPassoc-internal}
\alias{SNPassoc-internal}
\alias{association.fit}
\alias{extractPval}
\alias{extractPval.i}
\alias{SNPHWE}
\alias{GenotypeRate}
\alias{intervals.dif}
\alias{intervals.or}
\alias{is.quantitative}
\alias{haplo.inter.fit}
\alias{table.corner}
\alias{table.interaction}
\alias{crea.lab}
\alias{orderChromosome}
\alias{togeno}
\alias{expandsetupSNP}
\alias{pvalTest}
\alias{modelTest}
\alias{assoc}


\title{Internal SNPstat functions}
\description{Internal SNPassoc functions}
\usage{
association.fit(var, dep, adj, quantitative, type, level, 
       genotypingRate = 0) 
extractPval(x)
extractPval.i(i,x,pos,models)
SNPHWE(x)
GenotypeRate(x)
haplo.inter.fit(geno, var2, dep, adj = NULL, fam, 
          haplo.freq.min, ...)
crea.lab(x,pos.ini,cex,dist)
orderChromosome(x)
togeno(f,sep=sep,lab=lab)
expandsetupSNP(o)
pvalTest(dataX,Y,quantitative,type,genotypingRate)
modelTest(X,Y,quantitative,type,genotypingRate)
assoc(y,x,test="lrt",quantitative)

}

\details{These are not to be called by the user}
\keyword{internal}
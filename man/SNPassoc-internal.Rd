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

\title{Internal SNPstat functions}
\description{Internal SNPassoc functions}
\usage{
association.fit(var, dep, adj, quantitative, type, level, 
       genotypingRate = 0) 
extractPval(x)
extractPval.i(i,x,pos,models)
SNPHWE(x)
GenotypeRate(x)
intervals.dif(o, level, x.b, var, pval=TRUE, ...)
intervals.or(o, level, x.b, var, ...)
haplo.inter.fit(geno, var2, dep, adj = NULL, fam, 
          haplo.freq.min, ...)
crea.lab(x,pos.ini,cex,dist)
orderChromosome(x)
togeno(f,sep=sep,lab=lab)
expandsetupSNP(o)
}

\details{These are not to be called by the user}
\keyword{internal}
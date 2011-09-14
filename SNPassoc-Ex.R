pkgname <- "SNPassoc"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('SNPassoc')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("BonferroniSig")
### * BonferroniSig

flush(stderr()); flush(stdout())

### Name: Bonferroni.sig
### Title: Bonferroni correction of p values
### Aliases: Bonferroni.sig
### Keywords: utilities

### ** Examples

data(SNPs)
datSNP<-setupSNP(SNPs,6:40,sep="")
ans<-WGassociation(protein~1,data=datSNP,model="all")
Bonferroni.sig(ans, model="codominant", alpha=0.05, include.all.SNPs=FALSE)




cleanEx()
nameEx("GenomicControl")
### * GenomicControl

flush(stderr()); flush(stdout())

### Name: GenomicControl
### Title: Population substructure
### Aliases: GenomicControl
### Keywords: utilities

### ** Examples


data(SNPs) 
datSNP<-setupSNP(SNPs,6:40,sep="")
res<-scanWGassociation(casco,datSNP,model=c("do","re","log-add"))

# Genomic Control 
resCorrected<-GenomicControl(res)
plot(resCorrected)




cleanEx()
nameEx("HapMap")
### * HapMap

flush(stderr()); flush(stdout())

### Name: HapMap
### Title: SNPs from HapMap project
### Aliases: HapMap HapMap.SNPs.pos resHapMap
### Keywords: datasets

### ** Examples

data(HapMap)



cleanEx()
nameEx("TableNPer")
### * TableNPer

flush(stderr()); flush(stdout())

### Name: Table.N.Per
### Title: Descriptive sample size and percentage
### Aliases: Table.N.Per
### Keywords: utilities

### ** Examples

data(SNPs)
#sample size and percentage of cases and controls for each genotype 
#  Table.N.Per(SNPs$snp10001,SNPs$casco)

# The same table for a subset (males)
#  Table.N.Per(SNPs$snp10001,SNPs$casco,SNPs$sex=="Male")

# The same table assuming a dominant model
#  Table.N.Per(dominant(snp(SNPs$snp10001,sep="")),SNPs$casco,SNPs$sex=="Male")





cleanEx()
nameEx("Tablemeanse")
### * Tablemeanse

flush(stderr()); flush(stdout())

### Name: Table.mean.se
### Title: Descriptive sample size, mean, and standard error
### Aliases: Table.mean.se
### Keywords: utilities

### ** Examples

data(SNPs)
# sample size, mean age and standard error for each genotype
#  Table.mean.se(SNPs$snp10001,SNPs$protein)

# The same table for a subset (males)
#  Table.mean.se(SNPs$snp10001,SNPs$protein,SNPs$sex=="Male")

# The same table assuming a dominant model
#  Table.mean.se(dominant(snp(SNPs$snp10001,sep="")),SNPs$protein,SNPs$sex=="Male")





cleanEx()
nameEx("WGassociation")
### * WGassociation

flush(stderr()); flush(stdout())

### Name: WGassociation
### Title: Whole genome association analysis
### Aliases: WGassociation WGstats print.WGassociation
###   summary.WGassociation labels.WGassociation pvalues
###   pvalues.WGassociation codominant.WGassociation dominant.WGassociation
###   recessive.WGassociation overdominant.WGassociation
###   additive.WGassociation [.WGassociation
### Keywords: utilities

### ** Examples

data(SNPs)
datSNP<-setupSNP(SNPs,6:40,sep="")
ansAll<-WGassociation(protein~1,data=datSNP,model="all")

# In that case the formula is not required. You can also write:
# ansAll<-WGassociation(protein,data=datSNP,model="all")


#only codominant and log-additive
ansCoAd<-WGassociation(protein~1,data=datSNP,model=c("co","log-add"))

#for printing p values
print(ansAll)
print(ansCoAd)

#for obtaining a matrix with the p palues
pvalAll<-pvalues(ansAll)
pvalCoAd<-pvalues(ansCoAd)

# when all models are fitted and we are interested in obtaining p values for different genetic models

# codominant model
pvalCod<-codominant(ansAll)

# recessive model
pvalRec<-recessive(ansAll)

# and the same for additive, dominant or overdominant


#summary
summary(ansAll)

#for a detailed report
WGstats(ansAll)

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
summary(resHapMap)
plot(resHapMap)
print(resHapMap)




cleanEx()
nameEx("association")
### * association

flush(stderr()); flush(stdout())

### Name: association
### Title: Association analysis between a single SNP and a given phenotype
### Aliases: association print.snpOut
### Keywords: utilities

### ** Examples

data(SNPs)

# first, we create an object of class 'setupSNP'
datSNP<-setupSNP(SNPs,6:40,sep="")

# case-control study, crude analysis
association(casco~snp10001, data=datSNP)

# case-control study, adjusted by sex and arterial blood pressure
association(casco~sex+snp10001+blood.pre, data=datSNP)


# quantitative trait, crude analysis
association(log(protein)~snp10001,data=datSNP)
# quantitative trait, adjusted by sex
association(log(protein)~snp10001+sex,data=datSNP)


#
# Interaction analysis
#

# Interaction SNP and factor
association(log(protein)~snp10001*sex+blood.pre, data=datSNP, 
          model="codominant")

# Interaction SNP and SNP (codominant and codominant)
association(log(protein)~snp10001*factor(snp10002)+blood.pre, 
          data=datSNP, model="codominant")

# Interaction SNP and SNP (dominant and recessive)
association(log(protein)~snp10001*factor(recessive(snp100019))+blood.pre, 
          data=datSNP, model="dominant")





cleanEx()
nameEx("dscore")
### * dscore

flush(stderr()); flush(stdout())

### Name: dscore
### Title: Genetic risk allele score
### Aliases: dscore dscore.default dscore.setupSNP
### Keywords: utilities

### ** Examples


# example with 4 SNPs - the user gives the probabilities
MAFs <- c(0.1, 0.07, 0.2, 0.4)
dscore(MAFs)

# example with 4 SNPs - using setupSNP
data(SNPs)
myDat<-setupSNP(SNPs,6:10,sep="")
dscore(myDat)




cleanEx()
nameEx("getSignificantSNPs")
### * getSignificantSNPs

flush(stderr()); flush(stdout())

### Name: getSignificantSNPs
### Title: Extract significant SNPs from an object of class 'WGassociation'
### Aliases: getSignificantSNPs
### Keywords: utilities

### ** Examples

data(HapMap)
# resHapMap contains the results for a log-additive genetic model

# to get the significant SNPs for chromosome 12
getSignificantSNPs(resHapMap,chromosome=12)
# to get the significant SNPs for chromosome 5
getSignificantSNPs(resHapMap,5)
# to get the significant SNPs for chromosome X at level 1e-8
getSignificantSNPs(resHapMap,5,sig=1e-8)




cleanEx()
nameEx("haplointeraction")
### * haplointeraction

flush(stderr()); flush(stdout())

### Name: haplo.interaction
### Title: Haplotype interaction with a covariate
### Aliases: haplo.interaction print.haploOut
### Keywords: utilities

### ** Examples

# not Run
# data(SNPs)
# datSNP<-setupSNP(SNPs,6:40,sep="")
# res<-haplo.interaction(log(protein)~int(sex), data=datSNP,
#      SNPs.sel=c("snp100019","snp10001","snp100029"))
# res



cleanEx()
nameEx("inheritance")
### * inheritance

flush(stderr()); flush(stdout())

### Name: inheritance
### Title: Collapsing (or recoding) genotypes into different categories
###   (generally two) depending on a given genetic mode of inheritance
### Aliases: inheritance geneticModel codominant dominant recessive
###   overdominant additive
### Keywords: utilities

### ** Examples

data(SNPs)
dominant(snp(SNPs$snp10001,sep=""))
overdominant(snp(SNPs$snp10001,sep=""))



cleanEx()
nameEx("int")
### * int

flush(stderr()); flush(stdout())

### Name: int
### Title: Identify interaction term
### Aliases: int
### Keywords: utilities

### ** Examples

# Not Run
# data(SNPs)
# mod <- haplo.interaction(casco~int(sex)+blood.pre, data=SNPs,
# SNPs.sel=c("snp10001","snp10004","snp10005"))
#



cleanEx()
nameEx("interactionPval")
### * interactionPval

flush(stderr()); flush(stdout())

### Name: interactionPval
### Title: Two-dimensional SNP analysis for association studies
### Aliases: interactionPval print.SNPinteraction plot.SNPinteraction
### Keywords: utilities

### ** Examples


data(SNPs)
datSNP<-setupSNP(SNPs,6:40,sep="")

ansCod<-interactionPval(log(protein)~sex,datSNP)
print(ansCod)
plot(ansCod)




cleanEx()
nameEx("intervals")
### * intervals

flush(stderr()); flush(stdout())

### Name: intervals
### Title: Print ORs and 95% confidence intervals for an object of class
###   'haplo.glm'
### Aliases: intervals intervals.haplo.glm print.intervals
###   summary.haplo.glm
### Keywords: utilities

### ** Examples

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




cleanEx()
nameEx("isMonomorphic")
### * isMonomorphic

flush(stderr()); flush(stdout())

### Name: is.Monomorphic
### Title: Check whether a SNP is Monomorphic
### Aliases: is.Monomorphic
### Keywords: internal

### ** Examples

#
# data(SNPs)
# is.Monomorphic(SNPs$snp10001)
# is.Monomorphic(SNPs$snp100020)
#



cleanEx()
nameEx("makegeno")
### * makegeno

flush(stderr()); flush(stdout())

### Name: make.geno
### Title: Create a group of locus objects from some SNPs, assign to
###   'model.matrix' class.
### Aliases: make.geno
### Keywords: utilities

### ** Examples

## Not run:
data(SNPs)
# first, we create an object of class 'setupSNP'
datSNP<-setupSNP(SNPs,6:40,sep="")
geno<-make.geno(datSNP,c("snp10001","snp10002","snp10003"))
## End(Not run)





cleanEx()
nameEx("maxstat")
### * maxstat

flush(stderr()); flush(stdout())

### Name: maxstat
### Title: max-statistic for a 2x3 table
### Aliases: maxstat maxstat.default maxstat.table maxstat.setupSNP
###   maxstat.matrix print.maxstat
### Keywords: utilities

### ** Examples


# example from Sladek et al. (2007) for the SNP rs1111875 
 tt<-matrix(c(77,298,310,122,316,231),nrow=2,ncol=3,byrow=TRUE)
 maxstat(tt)

 data(SNPs)
 maxstat(SNPs$casco,SNPs$snp10001) 
 myDat<-setupSNP(SNPs,6:40,sep="")
 maxstat(myDat,casco)

 




cleanEx()
nameEx("odds")
### * odds

flush(stderr()); flush(stdout())

### Name: odds
### Title: Extract odds ratios, 95% CI and pvalues
### Aliases: odds
### Keywords: utilities

### ** Examples

 data(SNPs)
 datSNP<-setupSNP(SNPs,6:40,sep="")
 ans<-WGassociation(casco~1,data=datSNP,model="all")
 odds(ans)



cleanEx()
nameEx("permTest")
### * permTest

flush(stderr()); flush(stdout())

### Name: permTest
### Title: Permutation test analysis
### Aliases: permTest print.permTest plot.permTest
### Keywords: utilities

### ** Examples


data(SNPs)
datSNP<-setupSNP(SNPs,6:40,sep="")
ans<-scanWGassociation(casco~1,data=datSNP,model="co",nperm=1000)

# pPerm<-permTest(ans)
# print(pPerm)
# plot(pPerm)





cleanEx()
nameEx("plotMissing")
### * plotMissing

flush(stderr()); flush(stdout())

### Name: plotMissing
### Title: Plot of missing genotypes
### Aliases: plotMissing
### Keywords: utilities

### ** Examples

 data(SNPs)
 data(SNPs.info.pos) 
 ans<-setupSNP(SNPs,colSNPs=6:40,sep="")
 plotMissing(ans)
 
 # The same plot with the SNPs sorted by genomic position and 
 # showing the information about chromosomes

 ans<-setupSNP(SNPs,colSNPs=6:40,sort=TRUE,SNPs.info.pos,sep="") 
 plotMissing(ans)



cleanEx()
nameEx("qqpval")
### * qqpval

flush(stderr()); flush(stdout())

### Name: qqpval
### Title: Functions for inspecting population substructure
### Aliases: qqpval
### Keywords: utilities

### ** Examples

data(SNPs)
datSNP<-setupSNP(SNPs,6:40,sep="")
res<-scanWGassociation(casco,datSNP,model=c("do","re","log-add"))

# observed vs expected p values for recessive model
qqpval(recessive(res))




cleanEx()
nameEx("scanWGassociation")
### * scanWGassociation

flush(stderr()); flush(stdout())

### Name: scanWGassociation
### Title: Whole genome association analysis
### Aliases: scanWGassociation
### Keywords: utilities

### ** Examples


# Next steps may be very time consuming. So they are not executed

#data(HapMap)
#myDat<-setupSNP(HapMap, colSNPs=3:9307, sort = TRUE,
#   info=HapMap.SNPs.pos, sep="")
#resHapMap<-scanWGassociation(group~1, data=myDat, model="log")




cleanEx()
nameEx("setupSNP")
### * setupSNP

flush(stderr()); flush(stdout())

### Name: setupSNP
### Title: Convert columns in a dataframe to class 'snp'
### Aliases: setupSNP summary.setupSNP plot.setupSNP [.setupSNP
###   [[<-.setupSNP [<-.setupSNP $<-.setupSNP labels.setupSNP
### Keywords: utilities

### ** Examples


 data(SNPs)
 myDat<-setupSNP(SNPs,6:40,sep="")


#sorted SNPs and having genomic information
 data(SNPs.info.pos)
 myDat.o<-setupSNP(SNPs,6:40,sep="",sort=TRUE, info=SNPs.info.pos)

# summary
 summary(myDat.o)

# plot one SNP
  plot(myDat,which=2)




cleanEx()
nameEx("snp")
### * snp

flush(stderr()); flush(stdout())

### Name: snp
### Title: SNP object
### Aliases: snp is.snp as.snp reorder.snp summary.snp plot.snp
###   dominant.snp codominant.snp recessive.snp additive.snp print.snp
###   [.snp print.summary.snp
### Keywords: utilities

### ** Examples

# some examples of snp data in different formats

dat1  <- c("21", "21", "11", "22", "21",
                    "22", "22", "11", "11", NA)
ans1  <- snp(dat1,sep="")
ans1

dat2 <- c("A/A","A/G","G/G","A/G","G/G",
                    "A/A","A/A","G/G",NA)
ans2  <- snp(dat2,sep="/")
ans2

dat3 <- c("C-C","C-T","C-C","T-T","C-C",
                    "C-C","C-C","C-C","T-T",NA)
ans3 <- snp(dat3,sep="-")
ans3


dat4 <- c("het","het","het","hom1","hom2",
                    "het","het","hom1","hom1",NA)
ans4 <- snp(dat4,name.genotypes=c("hom1","het","hom2"))
ans4


# summary 
summary(ans3)

# plots

plot(ans3)
plot(ans3,type=pie)
plot(ans3,type=pie,label="SNP 10045")




cleanEx()
nameEx("sortSNPs")
### * sortSNPs

flush(stderr()); flush(stdout())

### Name: sortSNPs
### Title: Sort a vector of SNPs by genomic position
### Aliases: sortSNPs
### Keywords: utilities

### ** Examples

#
# data(SNPs)
# data(SNPs.info.pos)
# colSNPs.order<-sortSNPs(SNPs,c(6:40),SNPs.info.pos)
#



cleanEx()
nameEx("tableHWE")
### * tableHWE

flush(stderr()); flush(stdout())

### Name: tableHWE
### Title: Test for Hardy-Weinberg Equilibrium
### Aliases: tableHWE print.tableHWE
### Keywords: utilities

### ** Examples

data(SNPs)
ans<-setupSNP(SNPs,6:40,sep="")
res<-tableHWE(ans)
print(res)
#change the significance level showed in the flag column
print(res,sig=0.001)


#stratified analysis
res<-tableHWE(ans,ans$sex)
print(res)





### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')

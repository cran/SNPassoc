GenomicControl<-function(x, snp.sel=NULL)
 {
  if(!inherits(x,"WGassociation"))
   stop("x must be an object of class 'WGassociation'")

  if ("codominant"%in%names(x))
   stop(" \n  Genomic control cannot be applied to codominant model. Only to 2x2 tables") 

  WGchisq<-function(x) {
  pok<-1-x
  chisq<-qchisq(pok,1)
  chisq
  } 

  if (!is.null(snp.sel))
   {
    chisq.obs.sel<-apply(pvalues(x)[snp.sel,-1],2,WGchisq)
    num<-apply(chisq.obs.sel,2,median,na.rm=TRUE)
    lambda<-num/0.456
   }
  else
   {
    chisq.obs<-apply(pvalues(x)[,-1],2,WGchisq)
    num<-apply(chisq.obs,2,median,na.rm=TRUE)
    lambda<-num/0.456
   }

  lambdaOK<-ifelse(lambda<1,1,lambda)
  chisq.corrrected<-sweep(chisq.obs, 2, lambdaOK,FUN="/")
  pOK<-1-pchisq(chisq.corrrected,1)
  pOK[pOK==0]<-NA

  k<-length(names(x))
  attr(x,"pvalues")[,2:k]<-pOK

  cat("\n lambda: ",lambda,"\n")
  x
}

`assoc` <-
function(y,x,test="lrt",quantitative)
 {
 lrt<-function(m) 
  {
   if (m$family$family=="gaussian") {
    df1<-m$df.null
    df2<-m$df.residual
    df<-df1-df2 
    ans<-1-pchisq(((m$null.deviance-m$deviance))/(m$deviance/df2),df)
   }
   else {
    ans<-1-pchisq(m$null.deviance-m$deviance,m$df.null-m$df.residual)
   }
   ans 
  }
  if (length(levels(x))==1) {
    pval<-NA
  }
  else {
    if (test=="lrt") {
     if (quantitative)  
      pval<-lrt(glm(y~x,family="gaussian"))
     else  
      pval<-lrt(glm(y~x,family="binomial"))
    }
  }
  pval
 }


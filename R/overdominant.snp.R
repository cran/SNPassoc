`overdominant.snp` <-
function (o) 
{
  o<-overdominant.default(o)
  class(o)<-c("snp","factor")
  o
}


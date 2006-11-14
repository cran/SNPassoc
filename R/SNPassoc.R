
snp<-function (x, sep = "/", name.genotypes, reorder="common", remove.spaces = TRUE, 
       allow.partial.missing = FALSE) 
{

if(missing(name.genotypes))
 {
    alleles<-NULL
    x.d <- dim(x)
    x <- as.character(x)
    dim(x) <- x.d
    x[is.na(x)] <- ""
   
    if (remove.spaces) {
        xdim <- dim(x)
        x <- gsub("[ \t]", "", x)
        dim(x) <- xdim
        
    }

    if (!is.null(dim(x)) && ncol(x) > 1) 
        parts <- x[, 1:2]
    else {
        if (sep == "") 
            sep <- 1
        if (is.character(sep)) {
            part.list <- strsplit(x, sep)
            part.list[sapply(part.list, length) == 0] <- NA
            half.empties <- lapply(part.list, length) == 1
            part.list[half.empties] <- lapply(part.list[half.empties], 
                c, NA)
            empties <- is.na(x) | lapply(part.list, length) == 
                0
            part.list[empties] <- list(c(NA, NA))
            parts <- matrix(unlist(part.list), ncol = 2, byrow = TRUE)
        }
        else if (is.numeric(sep)) 
            parts <- cbind(substring(x, 1, sep), substring(x, 
                sep + 1, 9999))
        else stop(paste("I don't know how to handle sep=", sep))
    }
    mode(parts) <- "character"
    temp <- grep("^[ \t]*$", parts)
    parts[temp] <- NA
    if (!allow.partial.missing) 
        parts[is.na(parts[, 1]) | is.na(parts[, 2]), ] <- c(NA, 
            NA)
    alleles <- unique(c(na.omit(parts)))

    if(length(alleles)>2)
       stop("SNP must have only two alleles")


    tmp <- ifelse(is.na(parts[, 1]) & is.na(parts[, 2]), NA, 
        apply(parts, 1, paste, collapse = "/"))
    object <- factor(tmp)
    
    ll <- levels(object) <- na.omit(levels(object))

    if (length(ll)==4)
     {
       object[object==ll[3]]<-ll[2]
       object<-factor(object) 
     }
   
    control <- paste(rep(alleles[1],2),collapse="/")%in%ll
 
    if (sum(control)==0 & length(ll)==3)
     {
       object[object==ll[2]]<-ll[1]
       object<-factor(object) 
     }


    control <- paste(rep(alleles[2],2),collapse="/")%in%ll
 
    if (sum(control)==0 & length(ll)==3)
     {
       object[object==ll[3]]<-ll[2]
       object<-factor(object) 
     }

    if (length(object)==sum(is.na(object)))
     stop("choose the correct character separator to divide alleles") 

    class(object) <- c("snp","factor")
    object<-reorder.snp(object, ref=reorder)
    attr(object, "allele.names") <- alleles
}

else
 {
   if (any(is.na(match(x[!is.na(x)],name.genotypes))))
    stop("'name.genotypes' must match with the observed genotypes")
   x[x==name.genotypes[1]]<-"A/A"
   x[x==name.genotypes[2]]<-"A/B"
   x[x==name.genotypes[3]]<-"B/B"
   object<-as.factor(x)
   attr(object, "allele.names") <- c("A","B")
   class(object) <- c("snp","factor")

 }

object
}





reorder.snp<-function(x, ref="common", ...)
{
s<-x
if(!inherits(s,"snp"))
    stop("object must be of class 'snp'")


type<-charmatch(ref,c("common","minor"))

if (is.na(type))
 stop("ref must be either 'common' or 'minor'")


if (type==1)
 {
  class(s)<-"factor"
  tt <- table(s)
  if (length(tt) == 3 & min(tt) > 0) 
   {
    if (tt[1] < tt[3]) 
     {
      s <- relevel(relevel(s, 2), 3)
     } 
   }
  else 
   {
    if (length(unique(unlist(strsplit(names(tt)[1], "/")))) == 2) 
      {
       s <- relevel(s, 2)
      } 
   }
 } 

else
 {
  class(s)<-"factor"
  tt <- table(s)
  if (length(tt) == 3 & min(tt) > 0) 
   {
    if (tt[3] < tt[1]) 
     {
      s <- relevel(relevel(s, 2), 3)
     } 
   }
  else 
   {
    if (length(unique(unlist(strsplit(names(tt)[1], "/")))) == 2) 
      {
       s <- relevel(s, 2)
      } 
   }
 } 

class(s)<-c("snp","factor")
s
}



summary.snp<-function (object, print.out=TRUE, ...) 
{
    n <- length(object)
    nas <- is.na(object)
    n.typed <- n - sum(nas)
    ll <- levels(object)
    tbl <- table(object)
    tt <- c(tbl)
    names(tt) <- dimnames(tbl)[[1]]
    if (any(nas)) 
      { 
        tt.g <- c(tt, "NA's" = sum(nas))
        missing.allele<-sum(nas)/(sum(tt)+sum(nas))
      }
    else
      {
       tt.g <- tt
       missing.allele<-0
      }
    tt.g.prop <- prop.table(tbl)
    if (any(nas)) 
        tt.g.prop <- c(tt.g.prop, NA)
    ans.g <- cbind(frequency = tt.g, percentage = tt.g.prop * 
        100)
    if(print.out)
     {
       cat("Genotypes: \n")
       print(ans.g)
     }

    alle <- attr(object, "allele.names")
    alle1 <- length(grep(paste(alle[1], "/", sep = ""), as.character(object))) + 
        length(grep(paste("/", alle[1], sep = ""), as.character(object)))
    alle2 <- length(grep(paste(alle[2], "/", sep = ""), as.character(object))) + 
        length(grep(paste("/", alle[2], sep = ""), as.character(object)))
    tt.a <- c(alle1, alle2)
    tt.a.prop <- prop.table(tt.a)
    ans.a <- cbind(frequency = tt.a, percentage = tt.a.prop * 
        100)
    if (length(alle) > 1) {
        dimnames(ans.a)[[1]] <- alle
        if (any(nas)) 
            ans.a <- rbind(ans.a, "NA's" = c(2 * sum(nas), NA))
        if (print.out)  
         {
          cat("\n")
          cat("Alleles: \n")
          print(ans.a)
         }
    }
    else {
       if(print.out)
        {
         cat("\n")
         cat("Alleles: \n")
         cat("   Monomorphic \n") 
        }
    }
    pvalueHWE <- SNPHWE(ans.g[, 1])
    if(print.out)
     {
      cat("\n")
      cat("HWE (p value):", pvalueHWE, "\n")
     }
    ans <- list(allele.names = alle, allele.freq = ans.a, genotype.freq = ans.g, 
        n = n, n.typed = n.typed, HWE = pvalueHWE, missing.allele=missing.allele)
    class(ans) <- "summary.snp"
    invisible(ans)
}





plot.snp<-function (x, type = barplot, label, ...) 
{
    if (!inherits(x, "snp")) 
        stop("snp must be an object of class 'WGassociation'")
    if (missing(label)) 
        label <- deparse(substitute(x))
    old.mar <- par("mar")
    old.mfrow <- par("mfrow")
    on.exit(par(mar = old.mar, mfrow = old.mfrow))
    m <- m <- matrix(c(1, 2), nrow = 2, ncol = 1, byrow = TRUE)
    layout(m, heights = c(1, 5.5))
    par(mar = c(0, 0, 0, 0))
    xx <- summary(x)
    plot(c(1:5), rep(1, 5), ylim=c(0.1,1.6), type = "n", axes = FALSE, xlab = "", 
        ylab = "")
    text(1, 1.5, label, font = 2, adj = 0)
    crea.lab(xx$allele.freq, 1.6, 0.8, 0.25)
    crea.lab(xx$genotype.freq, 2.8, 0.8, 0.25)
    text(4.5, 1, paste("HWE (pvalue):", round(xx$HWE, 6)), cex = 0.8)
    par(mar = old.mar)
    type(xx$genotype.freq[, 1], ...)
}


crea.lab<-function (x, pos.ini, cex, dist) 
{
    n <- nrow(x)
    nn <- dimnames(x)[[1]]
    text(pos.ini + 0.1, 1.2, "frequency", cex = cex, adj = 0)
    text(pos.ini + 0.5, 1.2, "percentage", cex = cex, adj = 0)
    for (i in 1:n) {
        control <- (i - 1) * dist
        text(pos.ini, 1 - control, nn[i], cex = cex)
        text(pos.ini + 0.4, 1 - control, x[i, 1], adj = 1, cex = cex)
        text(pos.ini + 0.8, 1 - control, formatC(x[i, 2], 2,2,format="f"), adj = 1, 
            cex = cex)
    }
}


codominant<- function(o)
{
 if (length(unique(o[!is.na(o)]))>3)
    stop("variable should have 3 levels max")
 else factor(o)
}


dominant<-function (o) 
{
  codominant(o)
  if(length(levels(o))==3)
   {
    o[o == levels(o)[3]] <- levels(o)[2]
    levels(o)[2] <- paste(levels(o)[2:3], collapse = "-")
   } 
   factor(o)
}


recessive<-function (o) 
{
  codominant(o)
  if(length(levels(o))==3)
   {
    o[o == levels(o)[1]] <- levels(o)[2]
    levels(o)[2] <- paste(levels(o)[1:2], collapse = "-")
   } 
  else
   {
    allele<-attr(o,"allele.names")
    if(sum(levels(o)%in%paste(allele,collapse="/"))>0)
      {
       o[o == levels(o)[2]] <- levels(o)[1]
       levels(o)[1] <- paste(levels(o)[1:2], collapse = "-")
      }
   } 

  factor(o)
}


overdominant<-function (o) 
{
  codominant(o)
  if(length(levels(o))==3)
   {
    o[o == levels(o)[3]] <- levels(o)[1]
    levels(o)[2] <- paste(levels(o)[c(1, 3)], collapse = "-")
   } 
  else
   {
    allele<-attr(o,"allele.names")
    if(sum(levels(o)%in%paste(allele,collapse="/"))==0)
      {
       o[o == levels(o)[2]] <- levels(o)[1]
       levels(o)[1] <- paste(levels(o)[1:2], collapse = "-")
      }
   } 

  factor(o)
}



additive<-function (o) 
{
  codominant(o)
  if(length(levels(o))==3)
    o<-as.numeric(o)-1
  else
   {
    allele<-attr(o,"allele.names")
    if(sum(levels(o)%in%paste(allele,collapse="/"))>0)
      {
       o<-as.numeric(o)-1     
      }
    else
      {
       o<-as.numeric(o)-1
       o[o==1]<-2
      }
   } 
  o 
}



as.snp<- function (x, ...) 
 {
   if (is.snp(x)) x else snp(x, ...)
 } 

is.snp<-function(x)
{
 inherits(x, "snp")
}


setupSNP<-function(data, colSNPs, sort=FALSE, info, sep="/", ...)

{
 x<-data
 if (missing(x)) 
        stop("Required argument x is missing")
 if (!is.data.frame(x) & !is.matrix(x)) 
        stop("Argument x is not a data.frame or matrix")


 if (sort)
   {
     temp <- sortSNPs(x, colSNPs, info)
     pos <- temp$pos 
     info <- temp$dataSorted
     temp <- x[, pos]
     dataSNPs <- lapply(temp, snp, sep = sep, ...)
   }
 else
   { 
    dataSNPs <- lapply(x[,colSNPs],snp,sep=sep, ...)
   }
 
 datPhen<-x[,-colSNPs]
 ans<-data.frame(datPhen,dataSNPs)
 if(is.null(dim(datPhen)))
  names(ans)[1]<-names(x)[-colSNPs]

 label.SNPs <- names(dataSNPs)
 class(ans)<-c("setupSNP","data.frame")
 attr(ans,"row.names")<-1:length(ans[[1]])
 attr(ans,"label.SNPs")<-label.SNPs
 attr(ans,"colSNPs")<-c((length(ans)-length(label.SNPs)+1):length(ans))
 if (sort)
   attr(ans,"gen.info")<-info
 ans

}


summary.setupSNP<-function(object,...)
 {
  if (!inherits(object, "setupSNP")) 
        stop("object must be an object of class 'setupSNP'")
  colSNPs<-attr(object,"colSNPs")
  temp<-lapply(object[,colSNPs],expandsetupSNP)

  ans<-temp[[1]] 
  i<-2
  while (i <= length(temp))
   { 
    ans<-rbind(ans,temp[[i]])
    i<-i+1 
   } 
  dimnames(ans)[[1]]<-attr(object,"label.SNPs")
  out<-as.matrix(ans)
  dimnames(out)[[2]][4]<-"missing (%)"
  print(out,quote=FALSE,na.print="-")
  invisible(ans)
 }



expandsetupSNP<-function(o)
  {
   x<-summary(o, print.out=FALSE)
   allele<-rbind(x$allele.names)
   alleles<-paste(allele,collapse="/")
   if (length(allele)==1)
     alleles<-c(allele) 
    

   out<-data.frame(alleles=alleles,major.allele.freq=round(max(x$allele.freq[,2],na.rm=TRUE),1),
          HWE=round(x$HWE,6),"missing"=round(x$missing*100,1))
   out

  }



is.quantitative <- function(formula, data){
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m0 <- match(c("formula", "data"), names(mf), 0)
    mf <- mf[c(1, m0)]
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    ans<-length(table(unique(mf[,1])))>2
    ans
}


intervals.dif <- function(o, level, x.b, var, pval=TRUE, ...)
{

 x<-o
 z<-abs(qnorm((1-level)/2))

 mat.coef <- merge(x$coef, summary(x)$coef, by=0, all.x=TRUE, sort=FALSE)
 nom.pos <- data.frame(names(x$coef), ordre=1:length(x$coef))
 mat.ordre <- merge(nom.pos, mat.coef, by.x=1, by.y=1, all.x=TRUE, sort=FALSE)
 mat.ordre <- mat.ordre[order(mat.ordre$ordre),]
 
 a <- as.matrix(mat.ordre[,c("Estimate")])
 se <- as.matrix(mat.ordre[,c("Std. Error")])
 
 li <- a - (z * se)
 ls <- a + (z * se)
 
 if (missing(var))
 {
  focus <- nrow(a);
  m <- cbind(a[focus,],li[focus,],ls[focus,])
 
  if (pval)
  {
   p.as <- anova(x, x.b, test = "Chi")$"P(>|Chi|)"[2]
   m <- cbind(m, p.as)
   colnames(m) <- c("dif","lo","up","pval")
  }
  else
  {
   colnames(m) <- c("dif","lower","uppper")
  }
 }
 else
 {
  focus <- nrow(a) - length(levels(var)) + 2:length(levels(var))
  m <- cbind(a[focus,],li[focus,],ls[focus,])
  m <- rbind(c(0,NA,NA),m)
 
  if (pval)
  {
   p.as <- anova(x, x.b, test = "Chi")$"P(>|Chi|)"[2]
   m <- cbind(m, c(p.as,rep(NA,times=length(levels(var))-1)))
   colnames(m) <- c("dif","lower","upper","pval")
  }
  else
  {
   colnames(m) <- c("dif","lower","upper")
  }
 }
 
 list(m=m);
}

intervals.or <- function(o, level, x.b, var, ...)
{
 x<-o
 z<-abs(qnorm((1-level)/2))

 mat.coef <- merge(x$coef, summary(x)$coef, by=0, all.x=TRUE, sort=FALSE)
 nom.pos <- data.frame(names(x$coef), ordre=1:length(x$coef))
 mat.ordre <- merge(nom.pos, mat.coef, by.x=1, by.y=1, all.x=TRUE, sort=FALSE)
 mat.ordre <- mat.ordre[order(mat.ordre$ordre),]
 
 a <- as.matrix(mat.ordre[,c("Estimate")])
 se <- as.matrix(mat.ordre[,c("Std. Error")])
    
 or <- exp(a)
 li <- exp(a - z * se)
 ls <- exp(a + z * se)
 if(missing(var))
 {
  focus <- dim(a)[1.]
  or.ic <- round(cbind(or[focus,  ], li[focus,  ], ls[focus,  ]), 2.)
  or.ic[or.ic > 999.] <- NA
  p.as <- anova(x, x.b, test = "Chi")$"P(>|Chi|)"[2]
  or.ic <- cbind(or.ic, p.as)
  dimnames(or.ic) <- NULL
 }
    else
 {
  focus <- dim(a)[1.] - length(levels(var)) + 2:length(levels(var))
  or.ic <- round(cbind(or[focus,  ], li[focus,  ], ls[focus,  ]), 2)
  or.ic[or.ic > 999.] <- NA
  or.ic <- round(rbind(c(1, NA, NA), or.ic), 2)
  p.as <- anova(x, x.b, test = "Chi")$"P(>|Chi|)"[2]
  or.ic <- cbind(or.ic, c(p.as, rep(NA, times = length(levels(var)) - 1)))
  dimnames(or.ic) <- list(levels(var), c("   OR ", "lower", "upper", "p-value"))
 }
 list(or.ic = or.ic)
}




Table.mean.se <- function(var, dep, subset = !is.na(var))
{
    var <- as.factor(var)
    n <- ifelse(is.na(tapply(dep[subset],var[subset],FUN=length)),0,tapply(dep[subset],var[subset],FUN=length))
    me <- ifelse(is.na(tapply(dep[subset],var[subset],FUN=mean)),0,tapply(dep[subset],var[subset],FUN=mean))
    se <- ifelse(is.na(tapply(dep[subset],var[subset],FUN=function(x){sd(x)/sqrt(length(x))})),0,tapply(dep[subset],var[subset],FUN=function(x){sd(x)/sqrt(length(x))}))
    ta <- cbind(n=n,me=me,se=se)
    rownames(ta) <- names(n)
    list(tp = ta)
}



Table.N.Per <- function(var, dep, subset = !is.na(var))
    {
        var <- as.factor(var)
        dep <- as.factor(dep)
        ta <- table(var[subset], dep[subset])
        dimnames(ta) <- list(levels(var), levels(dep))
        per <- matrix(nrow = dim(ta)[1], ncol = dim(ta)[2])
        dimnames(per) <- list(levels(var[subset]), c("%", "%"))
        for(i in 1.:dim(ta)[1.]) {
                for(j in 1.:dim(ta)[2]) {
                        per[i, j] <- cbind(as.matrix(round((ta[i, j] * 100)/sum(ta[, j]), 1)))
                }
        }
        tp <- cbind(ta, per)
        tp <- tp[, order(rep(1:2, 2))]
        list(tp = tp)
    }



plotMissing<- function (x, print.labels.SNPs = TRUE, 
    main = "Genotype missing data", ...) 
{ 
   if(!inherits(x,"setupSNP"))
     stop("x must be an object of class 'setupSNP'")

   colSNPs<-attr(x,"colSNPs")
   data.SNPs <- t(x[colSNPs])
   label.SNPs<- attr(x,"label.SNPs")
   genInfo<-attr(x,"gen.info")

   data.Missing <- is.na(data.SNPs)
   old.xpd <- par("xpd")
   old.las <- par("las")
   par(xpd = TRUE)
   on.exit(par(xpd = old.xpd, las = old.las))
   image(1:nrow(data.Missing), 1:ncol(data.Missing), data.Missing, 
        col = c("white", "black"), ylab = "Individuals", xlab = ifelse(print.labels.SNPs, 
            "", "SNPs"), axes = !print.labels.SNPs)
   if (print.labels.SNPs) {
        axis(1, at = c(1:length(label.SNPs)), label = label.SNPs, 
            las = 3, cex.axis = 0.7)
        axis(2)
    }
    title(main, line = 3)
    if (!is.null(genInfo)) 
        n.snps <- table(genInfo[, 2])
    else n.snps <- length(label.SNPs)
    a <- c(0.5, cumsum(n.snps) + 0.5)
    b <- par("usr")
    if (!is.null(genInfo)) 
        col.ok <- c("black", rep("red", length(a) - 1))
    else col.ok <- c("black", rep("black", length(a) - 1))
    segments(a, b[3], a, b[4] + diff(b[3:4]) * 0.02, col = col.ok)
    abline(h = 0.5 + c(0, ncol(data.Missing)), xpd = FALSE)
    a <- par("usr")
    wh <- cumsum(c(0.5, n.snps))
    if (!is.null(genInfo)) {
        segments(a, b[3], a, b[4] + diff(b[3:4]) * 0.02, col = c("black", 
            rep("red", length(a) - 1)))
        names.geno <- unique(genInfo[, 2])
        n.gen <- length(names.geno)
        for (i in 1:n.gen) text(mean(wh[i + c(0, 1)]), a[4] + 
            (a[4] - a[3]) * 0.025, names.geno[i], srt = 45, cex = 0.8, 
            adj = 0.2)
    }
}







association<-function (formula, data, model = c("all"), model.interaction = c("codominant"), 
    subset, name.snp = NULL, quantitative = is.quantitative(formula, 
        data),  genotypingRate= 0, level = 0.95) 
{
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m0 <- match(c("formula", "data", "subset"), names(mf), 0)
    mf <- mf[c(1, m0)]
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    control <- unlist(lapply(mf, FUN = is.snp))
    if (sum(control) == 0) {
        stop("a variable of class 'snp' should be included in the model")
    }
    special <- c("strata")
    Terms <- if (missing(data)) 
        terms(formula, special)
    else terms(formula, special, data = data)
    strats <- attr(Terms, "specials")$strata
    ord <- attr(Terms, "order")
    if (any(ord > 1)) {
        varPos <- c(1:ncol(mf))[control][1]
        var <- mf[, varPos]
        dep <- mf[, 1]
        aux0 <- apply(attr(Terms, "factors"), 1, FUN = sum)
        aux <- names(aux0[aux0 == 2])
        control2 <- aux[aux != names(control)[control]]
        intPos <- c(1:ncol(mf))[names(mf) == control2]
        int <- mf[, intPos]
        if (!length(levels(int))) 
            stop("interaction term must be a factor")
        if (!is.null(strats)) 
            stop("interaction analysis does not support 'strata'")
        if (ncol(mf) > 3 & is.null(strats)) {
            adj <- data.frame(mf[, -c(1, varPos, intPos)])
            rownames(adj) <- 1:nrow(mf)
            variables <- attr(mt, "term.labels")
            varAdj <- variables[-c(varPos - 1, intPos - 1, length(variables))]
        }
        else {
            adj <- NULL
            varAdj <- NULL
            variables <- attr(mt, "term.labels")
        }
        int.nom <- variables[intPos - 1]
        if (is.null(name.snp)) 
            var.nom <- variables[varPos - 1]
        else var.nom <- name.snp
        model.type <- c("codominant", "dominant", "recessive", 
            "overdominant")
        m <- charmatch(model.interaction, model.type, nomatch = 0)
        if (m == 0) 
            stop("interaction analysis need a pre-defined model: codominant, dominant, recessive, or overdominant")
        if (length(table(var)) == 1) 
            stop("Monomorphic SNP")
        if (length(table(var)) == 2 & m > 1) 
            stop("SNP with only two genotypes. Codominant model is the only model that can be fitted")
        mod.inher <- switch(m, codominant, dominant, recessive, 
            overdominant)
        var <- mod.inher(var)
        res.corner <- table.corner(var, dep, adj, int, num.status = ifelse(quantitative, 
            1, 0), level)
        temp0 <- table.interaction(var, dep, adj, int, num.status = ifelse(quantitative, 
            1, 0), level)
        temp <- temp0$table
        p.interaction <- temp0$pval
        p.trend1 <- temp0$trend
        control.etiq <- ifelse(quantitative, 6, 5)
        etiq1 <- dimnames(temp)[[1]]
        aux0 <- dimnames(temp)[[2]]
        etiq2 <- aux0[seq(3, length(aux0), control.etiq)]
        ans <- list(NA)
        for (i in 1:nrow(temp)) {
            ans.i <- matrix(temp[i, ], nrow = length(etiq2), 
                ncol = control.etiq, byrow = TRUE)
            ans[[i]] <- data.frame(ans.i)
            dimnames(ans[[i]])[[1]] <- etiq2
            if (!quantitative) 
                names(ans[[i]]) <- c(aux0[1:2], "OR", "lower", 
                  "upper")
            else names(ans[[i]]) <- c(aux0[1:2], "se", "dif", 
                "lower", "upper")
        }
        names(ans) <- etiq1
        res.int1 <- ans
        temp0 <- table.interaction(int, dep, adj, var, num.status = ifelse(quantitative, 
            1, 0), level)
        temp <- temp0$table
        p.trend2 <- temp0$trend
        etiq1 <- dimnames(temp)[[1]]
        aux0 <- dimnames(temp)[[2]]
        etiq2 <- aux0[seq(3, length(aux0), control.etiq)]
        ans2 <- list(NA)
        for (i in 1:nrow(temp)) {
            ans.i <- matrix(temp[i, ], nrow = length(etiq2), 
                ncol = control.etiq, byrow = TRUE)
            ans2[[i]] <- data.frame(ans.i)
            dimnames(ans2[[i]])[[1]] <- etiq2
            if (!quantitative) 
                names(ans2[[i]]) <- c(aux0[1:2], "OR", "lower", 
                  "upper")
            else names(ans2[[i]]) <- c(aux0[1:2], "se", "dif", 
                "lower", "upper")
        }
        names(ans2) <- etiq1
        res.int2 <- ans2
        res <- list(res.corner, res.int1, res.int2, p.interaction, 
            p.trend1, p.trend2)
        interaction <- TRUE
    }
    else {
        type <- charmatch(model, c("codominant", "dominant", 
            "recessive", "overdominant", "log-additive", "all"))
        if (any(is.na(type))) 
            stop("model must be 'codominant','dominant','recessive','overdominant', \n                         'log-additive', 'all' or any combination of them")
        varPos <- c(1:ncol(mf))[control][1]
        dep <- mf[, 1]
        if (quantitative & !is.numeric(dep)) 
            stop("dependent variable should be numeric. It has more than two categories")
        var <- mf[, varPos]
        if (ncol(mf) > 2 & is.null(strats) | ncol(mf) > 3 & !is.null(strats)) {
            adj <- data.frame(mf[, -c(1, varPos, strats)])
            rownames(adj) <- 1:nrow(mf)
            variables <- attr(mt, "term.labels")
            varAdj <- variables[-c(varPos - 1, strats - 1)]
        }
        else {
            adj <- NULL
            varAdj <- NULL
            variables <- attr(mt, "term.labels")
        }
        if (is.null(name.snp)) 
            var.nom <- variables[varPos - 1]
        else var.nom <- name.snp
        dropx <- NULL
        if (length(strats)) {
            temp <- untangle.specials(Terms, "strata", 1)
            dropx <- c(dropx, temp$terms)
            if (length(temp$vars) == 1) 
                strata.keep <- mf[[temp$vars]]
            else strata.keep <- strata(mf[, temp$vars], shortlabel = TRUE)
            strats <- as.numeric(strata.keep)
            nstrats <- length(table(strats))
            res <- list()
            if (is.null(adj)) {
                for (i in 1:nstrats) {
                  res[[i]] <- association.fit(var[strats == i], 
                    dep[strats == i], adj, quantitative, type, 
                    level, genotypingRate)
                }
            }
            else {
                for (i in 1:nstrats) {
                  res[[i]] <- association.fit(var[strats == i], 
                    dep[strats == i], data.frame(adj[strats == 
                      i, ]), quantitative, type, level, genotypingRate)
                }
            }
            attr(res, "strata") <- levels(strata.keep)
        }
        else {
            res <- association.fit(var, dep, adj, quantitative, 
                type, level, genotypingRate)
        }
        interaction <- FALSE
    }
    class(res) <- "snpOut"
    attr(res, "varAdj") <- varAdj
    attr(res, "label.snp") <- var.nom
    if (interaction) 
        attr(res, "label.int") <- int.nom
    attr(res, "BigTable") <- FALSE
    attr(res, "Interaction") <- interaction
    res
}






association.fit<-function (var, dep, adj, quantitative, type, level, genotypingRate=0) 
{
    if (!quantitative) {
      if (length(unique(dep))==1)
      {
       res <- "Genot error"
      } 
      else 
       {
        co<-dom<-co<-dom<-rec<-over<-Ad<-NULL
        var <- as.factor(var)
        dep <- as.factor(dep)
        controlGeno <- GenotypeRate(var)
        if (genotypingRate > controlGeno)
          {
            res <- c(paste("Genot ", round(controlGeno, 1), "\\%", sep = ""))
          }

        else if (length(table(var)) == 1 | (length(table(var)) > 1 & 
            min(table(var)) == 0)) {
            res <- "Monomorphic"
        }
        else {
         if (length(levels(var)) == 3) {
                var.co <- codominant(var)
              if (any(type%in%6) | any(type%in%2))
                var.dom <- dominant(var)
              if (any(type%in%6) | any(type%in%3)) 
                var.rec <- recessive(var)
              if (any(type%in%6) | any(type%in%4))
                var.over <- overdominant(var)
            if (is.null(adj)) {
                  m.co <- glm(dep ~ var.co, family = binomial)
                  subset <- 1:length(var) %in% as.numeric(rownames(m.co$model))
                  m.b <- glm(dep ~ NULL, subset = subset, family = binomial)
              if (any(type%in%6) | any(type%in%2))
                  m.dom <- glm(dep ~ var.dom, subset = subset, 
                    family = binomial)
              if (any(type%in%6) | any(type%in%3))
                  m.rec <- glm(dep ~ var.rec, subset = subset, 
                    family = binomial)
              if (any(type%in%6) | any(type%in%4))
                  m.over <- glm(dep ~ var.over, subset = subset, 
                    family = binomial)
              if (any(type%in%6) | any(type%in%5))
                  m.ad <- glm(dep ~ as.numeric(var.co), subset = subset, 
                    family = binomial)
                }
            else {
                  m.co <- glm(dep ~ . + var.co, family = binomial, 
                    data = adj)
                  subset <- 1:length(var) %in% as.numeric(rownames(m.co$model))
                  m.b <- glm(dep ~ ., subset = subset, family = binomial, 
                    data = adj)
              if (any(type%in%6) | any(type%in%2))
                  m.dom <- glm(dep ~ . + var.dom, subset = subset, 
                    family = binomial, data = adj)
              if (any(type%in%6) | any(type%in%3))
                  m.rec <- glm(dep ~ . + var.rec, subset = subset, 
                    family = binomial, data = adj)
              if (any(type%in%6) | any(type%in%4))
                  m.over <- glm(dep ~ . + var.over, subset = subset, 
                    family = binomial, data = adj)
              if (any(type%in%6) | any(type%in%5))
                  m.ad <- glm(dep ~ . + as.numeric(var.co), subset = subset, 
                    family = binomial, data = adj)
                }


              if (any(type%in%6) | any(type%in%1))
               {
                temptp<-Table.N.Per(var.co, dep, subset)$tp
                co <- cbind(temptp, 
                  intervals.or(m.co, level, m.b, var)$or.ic, 
                  c(round(AIC(m.co), 1), NA, NA))
                if (any(temptp == 0) & is.null(adj)) 
                    {

                      pp<-fisher.test(dep,var.co)$p
                      co[1, 8] <- pp
                    }
               }
              if (any(type%in%6) | any(type%in%2))
               {
                temptp<-Table.N.Per(var.dom, dep, subset)$tp
                dom <- cbind(temptp, 
                  intervals.or(m.dom, level, m.b, var.dom)$or.ic, 
                  c(round(AIC(m.dom), 1), NA))
                if (any(temptp == 0) & is.null(adj)) 
                    {
                      pp<-fisher.test(dep,var.dom)$p
                      dom[1, 8] <- pp
                    }
               } 
              if (any(type%in%6) | any(type%in%3))
               {
                temptp<-Table.N.Per(var.rec, dep, subset)$tp
                rec <- cbind(temptp, 
                  intervals.or(m.rec, level, m.b, var.rec)$or.ic, 
                  c(round(AIC(m.rec), 1), NA))
                if (any(temptp == 0) & is.null(adj)) 
                    {
                      pp<-fisher.test(dep,var.rec)$p
                      rec[1, 8] <- pp
                    }
               }
              if (any(type%in%6) | any(type%in%4))
               {
                temptp<-Table.N.Per(var.over, dep, subset)$tp
                over <- cbind(temptp, 
                  intervals.or(m.over, level, m.b, var.over)$or.ic, 
                  c(round(AIC(m.over), 1), NA))
                if (any(temptp == 0) & is.null(adj)) 
                    {
                      pp<-fisher.test(dep,var.over)$p
                      over[1, 8] <- pp
                    }
               }
              if (any(type%in%6) | any(type%in%5))
               {
                temptp<-Table.N.Per(var.co, dep, subset)$tp
                totals<-round(table(dep),1)
                prop.totals<-round(100*prop.table(totals),1)
                ansTot<-c(totals[1],prop.totals[1], totals[2],prop.totals[2])
                Ad <- c(ansTot, intervals.or(m.ad, 
                  level, m.b)$or.ic, round(AIC(m.ad), 1))
                if (any(temptp == 0) & is.null(adj)) 
                    {
                      pp<-fisher.test(dep,var.co)$p
                      Ad[8] <- pp
                    }
               } 
             res<-NULL 
             if(!is.null(co))
               res<-rbind(Codominant=rep(NA,9),co)
             if(!is.null(dom))
               res<-rbind(res,Dominant=rep(NA,9),dom)  
             if(!is.null(rec))
               res<-rbind(res,Recessive=rep(NA,9),rec)
             if(!is.null(over))
               res<-rbind(res,Overdominant=rep(NA,9),over)
             if(!is.null(Ad))
               res<-rbind(res,"log-Additive"=rep(NA,9),"0,1,2"=Ad)

             dimnames(res)[[2]][5:9] <- c("OR","lower","upper","p-value","AIC")
               
            }
            else if (length(levels(var)) == 2) {
                var.co <- codominant(var)
                if (is.null(adj)) {
                  m.co <- glm(dep ~ var.co, family = binomial)
                  subset <- 1:length(var) %in% as.numeric(rownames(m.co$model))
                  m.b <- glm(dep ~ NULL, subset = subset, family = binomial)
#                  m.ad <- glm(dep ~ as.numeric(var.co), subset = subset, 
#                    family = binomial)
                }
                else {
                  m.co <- glm(dep ~ . + var.co, family = binomial, 
                    data = adj)
                  subset <- 1:length(var) %in% as.numeric(rownames(m.co$model))
                  m.b <- glm(dep ~ ., subset = subset, family = binomial, 
                    data = adj)
#                  m.ad <- glm(dep ~ . + as.numeric(var.co), subset = subset, 
#                    family = binomial, data = adj)
                }
                co <- cbind(Table.N.Per(var.co, dep, subset)$tp, 
                  intervals.or(m.co, level, m.b, var)$or.ic, 
                  c(round(AIC(m.co), 1), NA))
#                Ad <- c(rep(NA, times = 4), intervals.or(m.ad, level, m.b)$or.ic, round(AIC(m.ad), 1))

                totals<-table(dep)
                prop.totals<-round(100*prop.table(totals),1)
                ansTot<-c(totals[1],prop.totals[1], totals[2],prop.totals[2])
                Ad <- c(ansTot, intervals.or(m.co, level, m.b)$or.ic, round(AIC(m.co), 1))
                  
                Ad[8]<-NA 

               if(any(Table.N.Per(var.co, dep, subset)$tp==0) & is.null(adj))
                {
                 pp<-fisher.test(dep,var.co)$p
                 Ad[8]<-pp
                 co[1,8]<-pp 
                }

                res <- rbind(Codominant=rep(NA,9),co,"log-Additive"=rep(NA,9), "0,1,2"=Ad)
                dimnames(res)[[2]][5:9] <- c("OR","lower","upper","p-value","AIC")
            }
        }
       }
    }
    else {    # quantitative trait
        co<-dom<-co<-dom<-rec<-over<-Ad<-NULL 
        var <- as.factor(var)
        controlGeno <- GenotypeRate(var)
        if (genotypingRate > controlGeno)
          {
            res <- c(paste("Genot ", round(controlGeno, 1), "\\%", sep = ""))
          }

        else if (length(table(var)) == 1 | (length(table(var)) > 1 & 
            min(table(var)) == 0)) {
            res <- "Monomorphic"
        }
        else {
            if (length(levels(var)) == 3) {
              var.co <- codominant(var)
              if (any(type%in%6) | any(type%in%2))
                var.dom <- dominant(var)
              if (any(type%in%6) | any(type%in%3))
                var.rec <- recessive(var)
              if (any(type%in%6) | any(type%in%4))
                var.over <- overdominant(var)
                if (is.null(adj)) {
                  m.co <- glm(dep ~ var.co, family = gaussian)
                  subset <- 1:length(var) %in% as.numeric(rownames(m.co$model))
                  m.b <- glm(dep ~ NULL, subset = subset, family = gaussian)
              if (any(type%in%6) | any(type%in%2))
                  m.dom <- glm(dep ~ var.dom, subset = subset, 
                    family = gaussian)
              if (any(type%in%6) | any(type%in%3))
                  m.rec <- glm(dep ~ var.rec, subset = subset, 
                    family = gaussian)
              if (any(type%in%6) | any(type%in%4))
                  m.over <- glm(dep ~ var.over, subset = subset, 
                    family = gaussian)
              if (any(type%in%6) | any(type%in%5))
                  m.ad <- glm(dep ~ as.numeric(var.co), subset = subset, 
                    family = gaussian)
                }
                else {
                  m.co <- glm(dep ~ . + var.co, family = gaussian, 
                    data = adj)
                  subset <- 1:length(var) %in% as.numeric(rownames(m.co$model))
                  m.b <- glm(dep ~ ., subset = subset, family = gaussian, 
                    data = adj)
              if (any(type%in%6) | any(type%in%2))
                  m.dom <- glm(dep ~ . + var.dom, subset = subset, 
                    family = gaussian, data = adj)
              if (any(type%in%6) | any(type%in%3))
                  m.rec <- glm(dep ~ . + var.rec, subset = subset, 
                    family = gaussian, data = adj)
              if (any(type%in%6) | any(type%in%4))
                  m.over <- glm(dep ~ . + var.over, subset = subset, 
                    family = gaussian, data = adj)
              if (any(type%in%6) | any(type%in%5))
                  m.ad <- glm(dep ~ . + as.numeric(var.co), subset = subset, 
                    family = gaussian, data = adj)
                }
              if (any(type%in%6) | any(type%in%1))
                co <- cbind(Table.mean.se(var.co, dep, subset)$tp, 
                  intervals.dif(m.co, level, m.b, var)$m, AIC = c(round(AIC(m.co), 
                    1), NA, NA))
              if (any(type%in%6) | any(type%in%2))
                dom <- cbind(Table.mean.se(var.dom, dep, subset)$tp, 
                  intervals.dif(m.dom, level, m.b, var.dom)$m, 
                  AIC = c(round(AIC(m.dom), 1), NA))
              if (any(type%in%6) | any(type%in%3))
                rec <- cbind(Table.mean.se(var.rec, dep, subset)$tp, 
                  intervals.dif(m.rec, level, m.b, var.rec)$m, 
                  AIC = c(round(AIC(m.rec), 1), NA))
              if (any(type%in%6) | any(type%in%4))
                over <- cbind(Table.mean.se(var.over, dep, subset)$tp, 
                  intervals.dif(m.over, level, m.b, var.over)$m, 
                  AIC = c(round(AIC(m.over), 1), NA))
              if (any(type%in%6) | any(type%in%5))
                Ad <- c(rep(NA, 3), intervals.dif(m.ad, level, 
                  m.b)$m, AIC(m.ad))

             res<-NULL 
             if(!is.null(co))
               res<-rbind(Codominant=rep(NA,8),co)
             if(!is.null(dom))
               res<-rbind(res,Dominant=rep(NA,8),dom)  
             if(!is.null(rec))
               res<-rbind(res,Recessive=rep(NA,8),rec)
             if(!is.null(over))
               res<-rbind(res,Overdominant=rep(NA,8),over)
             if(!is.null(Ad))
               res<-rbind(res,"log-Additive"=rep(NA,8),"0,1,2"=Ad)
              dimnames(res)[[2]][4:8] <- c("dif","lower","upper","p-value","AIC")
            }


           else if (length(levels(var)) == 2) {
                var.co <- codominant(var)
                if (is.null(adj)) {
                  m.co <- glm(dep ~ var.co, family = gaussian)
                  subset <- 1:length(var) %in% as.numeric(rownames(m.co$model))
                  m.b <- glm(dep ~ NULL, subset = subset, family = gaussian)
#                  m.ad <- glm(dep ~ as.numeric(var.co), subset = subset, 
#                    family = gaussian)
                }
                else {
                  m.co <- glm(dep ~ . + var.co, family = gaussian, 
                    data = adj)
                  subset <- 1:length(var) %in% as.numeric(rownames(m.co$model))
                  m.b <- glm(dep ~ ., subset = subset, family = gaussian, 
                    data = adj)
#                  m.ad <- glm(dep ~ . + as.numeric(var.co), subset = subset, 
#                    family = gaussian, data = adj)
                }
                co <- cbind(Table.mean.se(var.co, dep, subset)$tp, 
                  intervals.dif(m.co, level, m.b, var)$m, AIC = c(AIC(m.co), 
                    NA))
#                Ad <- c(rep(NA, 3), intervals.dif(m.ad, level, m.b)$m, round(AIC(m.ad), 1))
                Ad <- c(rep(NA, 3), intervals.dif(m.co, level,m.b)$m, round(AIC(m.co), 1)) 

                Ad[7]<-NA 

                if(any(Table.mean.se(var.co, dep, subset)$tp==0) & is.null(adj))
                 {
                  pp<-fisher.test(dep,var.co)$p
                  Ad[7]<-pp
                  co[1,7]<-pp 
                 }

                res <- rbind(Codominant=rep(NA,8),co,"log-Additive"=rep(NA,8), "0,1,2"=Ad)
                dimnames(res)[[2]][4:8] <- c("dif","lower","upper","p-value","AIC")
            }
        }
    }
    res
}


print.snpOut<-function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    temp <- attr(x, "label.snp")
    temp2 <- gsub("(snp\\()", "", temp)
    attr(x, "label.snp") <- gsub("\\)", "", temp2)
    cat("\n")
    if (!attr(x, "Interaction")) {
        if (is.null(attr(x, "strata"))) {
            cat("SNP:", attr(x, "label.snp"), " adjusted by:", 
                attr(x, "varAdj"), "\n")
            class(x) <- NULL
            attr(x, "varAdj") <- attr(x, "label.snp") <- attr(x, 
                "BigTable") <- attr(x, "Interaction") <- NULL
            if (!is.Monomorphic(x)) {
                print(x, na.print = "", digits = digits, quote = FALSE)
            }
            else cat("Monomorphic\n")
        }
        else {
            nstrat <- length(attr(x, "strata"))
            for (i in 1:nstrat) {
                cat("       strata:", attr(x, "strata")[i], "\n")
                cat("SNP:", attr(x, "label.snp"), " adjusted by:", 
                  attr(x, "varAdj"), "\n")
                class(x) <- NULL
                attr(x, "varAdj") <- attr(x, "label.snp") <- attr(x, 
                  "BigTable") <- attr(x, "Interaction") <- NULL
                if (!is.Monomorphic(x)) {
                  print(x[[i]], na.print = "", digits = digits, 
                    quote = FALSE)
                }
                else cat("Monomorphic\n")
                cat("\n")
            }
        }
    }
    else {
        cat("      SNP:", attr(x, "label.snp"), " adjusted by:", 
            attr(x, "varAdj"), "\n")
        cat(" Interaction \n")
        cat("---------------------\n")
        print(x[[1]], digits = digits)
        cat("\n")
        cat("p interaction:", x[[4]], "\n")
        cat("\n", paste(attr(x, "label.int"), "within",attr(x, "label.snp")), "\n")
        cat("---------------------\n")
        for (i in 1:length(x[[2]])) {
            cat(names(x[[2]])[i], "\n")
            print(x[[2]][[i]], digits = digits)
            cat("\n")
        }
        cat("p trend:", x[[5]], "\n")
        cat("\n", paste(attr(x, "label.snp"),"within",attr(x, "label.int")), "\n")
        cat("---------------------\n")
        for (i in 1:length(x[[3]])) {
            cat(names(x[[3]])[i], "\n")
            print(x[[3]][[i]], digits = digits)
            cat("\n")
        }
        cat("p trend:", x[[6]], "\n")
    }
}



GenotypeRate<-function(x)
{
 temp<-sum(!is.na(x))/length(x)
 ans<-temp*100
 ans
}




is.Monomorphic<- function (x) 
{
   ans<-FALSE   
   if (length(x)==1)
    {
    if (x[1] == "Monomorphic") 
        ans <- TRUE
    }
   else
      ans <- length(table(x)[table(x) > 0]) == 1
   return(ans)
}




WGassociation<-function (formula, data, model=c("all"), quantitative = is.quantitative(formula, 
    data), genotypingRate = 80, level=0.95)
 {
  
    if(!inherits(data,"setupSNP"))
     stop("data must be an object of class 'setupSNP'")

    if (length(attr(data,"colSNPs")) > 2000 & (length(model) > 1 | any(model%in%"all"))) 
        stop("Select only one genetic model when more than 2000 SNPs are analyzed \n or use 'fastWGassociation' function")
 
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m0 <- match(c("formula", "data"), names(mf), 0)
    mf <- mf[c(1, m0)]

#
# aceptar respuesta sin formula
#    
	if( length(grep("~",mf[[2]]))==0){
		formula<-as.formula(paste(mf[[2]],"~1",sep=""))
		formula.1<- list(formula)
		mode(formula.1)<-"call"
		mf[2]<-formula.1
	}

    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    temp0 <- as.character(mt)
    adj <- paste(temp0[2], temp0[1], temp0[3])

    Terms <- if (missing(data)) 
        terms(formula)
    else terms(formula, data = data)
    ord <- attr(Terms, "order")
    if (any(ord > 1)) 
     stop("interaction term is not implemented")


    association.i<-function(snp.i,adj,data,model,quantitative,genotypingRate,level)
     {
       association(as.formula(paste(adj,"+",snp.i)), data=data,
          model=model, quantitative=quantitative, genotypingRate=
          genotypingRate, level=level)
     }
 

    colSNPs<-attr(data,"colSNPs")
    if (is.vector(colSNPs) & length(colSNPs) > 1) 
        dataSNPs <- data[, colSNPs]
    else stop("data should have an attribute called 'colSNPs'. Try again 'setupSNP' function")
    
    type<-charmatch(model,c("codominant","dominant","recessive","overdominant","log-additive","all"))
    type<-sort(type)

    if (any(type%in%6))
      type<-1:6
  
    if(any(is.na(type)))
     stop("model must be 'codominant','dominant','recessive','overdominant', 
                     'log-additive', 'all' or any combination of them")

    SNPs<-attr(data,"label.SNPs")

    out<-lapply(SNPs, association.i, adj=adj, data=data, model=model,
            quantitative=quantitative, genotypingRate=
            genotypingRate, level=level) 

    class(out)<-"WGassociation"
    attr(out,"label.SNPs")<-attr(data,"label.SNPs")
    attr(out,"models")<-type
    attr(out,"quantitative")<-quantitative
    attr(out,"pvalues")<-extractPval(out)
    attr(out,"gen.info")<-attr(data,"gen.info")
    attr(out,"whole")<-attr(data,"whole")
    attr(out,"colSNPs")<-attr(data,"colSNPs")
    out
}




extractPval<-function(x)
 {
   
   models<-attr(x,"models")
   if(length(models)==6)
    models<-c(1:5) 
   quantitative<-attr(x,"quantitative")
   pos<-ifelse(quantitative,7,8)

   ans<-t(data.frame(lapply(1:length(x),extractPval.i,x=x,pos=pos,models=models)))

   ans<-data.frame(ans)
   for (i in 2:ncol(ans))
    ans[,i]<-as.numeric(as.character(ans[,i]))

   dimnames(ans)[[1]]<-attr(x,"label.SNPs")
   dimnames(ans)[[2]]<-c("comments",c("codominant","dominant","recessive","overdominant","log-additive")[models])
   ans

 }



extractPval.i<-function (i, x, pos, models) 
{
    tt<-x[[i]]
    if (is.null(nrow(tt)))
      control<-1000
    else
      control<-nrow(tt)     

    if(length(models)*2>control)
     {
       tt <- tt[, pos]
       ans <- tt[!is.na(tt)][1]
       ans <- c(NA,ans, rep(NA, length(models) - 1))
     }
    else
     {
      if (!is.null(dim(tt))) {
        if (length(models) == 1) {
            tt <- tt[, pos]
            ans <- c(NA, tt[!is.na(tt)][1])
        }
        else {
            tt <- tt[, pos]
            ans <- c(NA, tt[!is.na(tt)])
            if ((length(ans) - 1) < length(models)) 
#                ans <- c(ans, rep(NA, length(models) - 1))
                ans <- c(ans, rep(NA, length(models) - (length(ans) - 1)))
          }
       }
      else if (!is.na(charmatch("Geno", tt))) 
        ans <- c(tt[1], rep(NA, length(models)))
      else ans <- c("Monomorphic", rep(NA, length(models)))
    }
    ans
}


getSignificantSNPs<-function (x, chromosome, model, sig = 1e-15) 
{
    if (!inherits(x, "WGassociation")) 
        stop("x must be an object of class 'WGassociation'")

    if(is.null(attr(x, "gen.info")))
     {
        pvalues <- attr(x, "pvalues")
        mm<-charmatch(model,dimnames(pvalues)[[2]])
        if (is.na(mm))
          stop("model selected is not correct")  
        sel2<-pvalues[,mm]<=sig
        SNPs.sel <- pvalues[sel2, ]
        pos.sel <- attr(x, "colSNPs")[sel2]
        out <- list(names = dimnames(SNPs.sel)[[1]], column = pos.sel)
     }

    else
     {
      if (!chromosome %in% c(1:22) & chromosome != "X") 
        stop("chromosome should be either a number between 1 and 22 or X")
      if (chromosome == "X") 
        chromosome <- 23
      gen.info <- attr(x, "gen.info")
      pvalues <- attr(x, "pvalues")
      chrs <- gen.info[, 2]
      chr.l <- unique(chrs)
      chr <- chr.l[orderChromosome(chr.l)]
      sel <- chr[chromosome]
      SNPs <- gen.info[gen.info[, 2] %in% sel, 1]
      sel2 <- dimnames(pvalues)[[1]] %in% SNPs & !is.na(pvalues[, 
        2]) & pvalues[, 2] <= sig
      SNPs.sel <- pvalues[sel2, ]
      pos.sel <- attr(x, "colSNPs")[sel2]
      out <- list(names = dimnames(SNPs.sel)[[1]], column = pos.sel)
     }
    out
}



print.WGassociation<- function (x, digits = 5, ...) 
{
    if (!inherits(x, "WGassociation")) 
        stop("x must be an object of class 'WGassociation'")
    ans <- attr(x, "pvalues")
    ans[, -1] <- round(ans[, -1], digits)
    out <- as.matrix(ans)
    out[,1]<-gsub("\\\\","",out[,1])
    print(out, quote = FALSE, na.print = "-", ...)
    invisible(ans)
}




summary.WGassociation<-function (object, ...) 
{
    x <- object
    if (!inherits(x, "WGassociation")) 
        stop("x must be an object of class 'WGassociation'")

    if (!is.null(attr(x,"fast")))
       stop("\n summary is implemented only for 'WGassociation' function")
    nn <- names(x)
    attr(x, "label.SNPs") <- attr(x, "models") <- attr(x, "quantitative") <- attr(x, 
        "pvalues") <- attr(x, "gen.info") <- attr(x, "whole") <- attr(x, 
        "colSNPs") <- NULL
    class(x) <- NULL
    print(x, na.print = "", ...)
    invisible(x)
}



"plot.WGassociation" <-
function (x, alpha = 0.05, plot.all.SNPs = FALSE, print.label.SNPs = TRUE, 
    cutPval = c(0, 1e-10, 1), whole, ylim.sup = ifelse(is.null(attr(x,"fast")),1e-40,
    1e-30), col.legend = c("red","gray60"), sort.chromosome=TRUE, centromere, ...) 
{
    if (!inherits(x, "WGassociation")) 
        stop("x must be an object of class 'WGassociation'")

    if (missing(whole))
     {
       if (ncol(attr(x, "pvalues"))==2 & length(unique(attr(x,"gen.info")[,2]))>10)
         whole <- TRUE
       else 
         whole <- FALSE 
     }
    if (!whole) {
        x <- attr(x, "pvalues")
        ylims<-range(-log(x[,-1]),na.rm=TRUE)
        control.Mono <- grep("Monomorphic", as.character(x[, 
            1]))
        control.Geno <- grep("Genot", as.character(x[, 1]))

        if (plot.all.SNPs | length(c(control.Mono, control.Geno))==0) 
            ans <- x
        else ans <- x[-c(control.Mono, control.Geno), ]
        old.mar <- par("mar")
        old.mfrow <- par("mfrow")
        on.exit(par(mar = old.mar, mfrow = old.mfrow))
        m <- matrix(c(1, 2, 1, 3, 1, 4, 1, 5, 1, 6, 1, 7, 1, 
            8, 1, 9, 1, 10, 1, 11, 1, 12), nrow = 11, ncol = 2, 
            byrow = TRUE)
        layout(m, heights = c(0.3, 1, 0.2, 1, 0.2, 1, 0.2, 1, 
            0.2, 1, 0.7), widths = c(0.05, 1))
        par(mar = c(0, 0, 0, 0))
        n.models <- ncol(ans) - 1
        if (n.models>2)
         control.labelY <- 6.5 - ceiling(n.models/2)
        else if (n.models==2)
         control.labelY <- 4.5
        else
          control.labelY <- 5
        plot(rep(1, 5), 1:5, type = "n", axes = FALSE, xlab = "", 
            ylab = "")
        text(1, control.labelY, "-log (p value)", font = 2, srt = 90, 
            cex = 2, adj=1)
        models <- c("codominant", "dominant", "recessive", "overdominant", 
            "log-additive")
        ok <- 2
        for (i in 2:6) {
            par(mar = c(0, 0, 0, 0))
            if (any(models[i - 1] %in% names(ans))) {
                plot(1:3, rep(1, 3), type = "n", axes = FALSE, 
                  xlab = "", ylab = "")
                if (i == 2) 
                  legend(2.5, 1.5, c("Nominal p value", "Bonferroni correction"), 
                    lty = c(2, 2), col = c("pink1", "red"), cex = 1, 
                    bty = "n", y.intersp = 0.8)
                text(2, 1, models[i - 1], col = "blue", cex = 1.5, 
                  font = 2)
                pval <- ans[, ok]
                xx <- -log(pval)

                par(mar = c(0.3, 0, 0, 0))
                plot(c(1:length(xx)), xx, type = "b", axes = FALSE, 
                  ylab = "", xlab = "", ylim=ylims)
                axis(1, at = c(1:length(xx)), label = rep("", 
                  length(xx)), pos = 0)
                control.y <- ceiling(seq(0,ceiling(ylims[2]),length=6))
                axis(2, at = control.y, label = control.y, 
                  pos = 1)
                cut <- -log(alpha)
                segments(1, cut, length(xx), cut, col = "pink2", 
                  lty = 2)
                pvalues <- xx[!is.na(xx)]
                cut.p <- alpha/length(pvalues)
                cut.trans <- -log(cut.p)
                if (cut.trans > max(xx, na.rm = TRUE)) {
                  cat("Warning: Any SNP is statistically significant after \n                     Bonferroni Correction under", 
                    models[i - 1], "model \n")
                }
#                else {
                  segments(1, cut.trans, length(xx), cut.trans, 
                    col = "red", lty = 2)
#                }
                ok <- ok + 1
            }
        }
        if (ok < 6) 
            plot(1:3, rep(1, 3), type = "n", axes = FALSE, xlab = "", 
                ylab = "")
        plot(c(1:nrow(ans)), rep(0, nrow(ans)), type = "n", axes = FALSE)
        if (print.label.SNPs) {
            text(c(1:nrow(ans)), rep(1, nrow(ans)), dimnames(ans)[[1]], 
                srt = 90, adj = 1, xpd = TRUE, ...)
        }
        else {
            text(nrow(ans)/2, 0, "SNPs", cex = 2)
        }
    }
    else {
        gen.info <- attr(x, "gen.info")
        chrs <- gen.info[, 2]
        chr.l <- unique(chrs)

        if (sort.chromosome)
          chr <- chr.l[orderChromosome(chr.l)]
        else
         chr<-unique(chrs)

        pval <- attr(x, "pvalues")
        limits <- c(min(gen.info[, 3]), max(gen.info[, 3]))
        if (missing(centromere)) {
            centro <- c(1.23e+08, 93500000, 91500000, 5.1e+07, 
                47700000, 60500000, 58900000, 4.5e+07, 4.9e+07, 
                4e+07, 5.3e+07, 35400000, 1.5e+07, 15600000, 
                1.7e+07, 3.9e+07, 2.4e+07, 1.6e+07, 28500000, 
                27800000, 12200000, 11900000, 58500000, 1e+07)
        }
        else {
            centro <- centromere
        }
        n.chr <- length(chr)
        old.mfrow <- par("mfrow")
        old.mar <- par("mar")
        old.xpd <- par("xpd")
        on.exit(par(mfrow = old.mfrow, mar = old.mar))
        par(mfrow = c(26, 1))
        par(mar = c(0.5, 5, 0, 3))
        par(xpd = TRUE)
        plot(c(1:3), rep(1, 3), axes = FALSE, xlab = "", ylab = "", 
            type = "n")
        text(2, 1, paste("Genetic model:", names(pval)[2]), adj = 0.5, 
            font = 2, cex = 2)
        plot(c(1:3), rep(1, 2, 3), axes = FALSE, xlab = "", ylab = "", 
            type = "n")
        text(1, 1.2, "p value", adj = 0)
        legend(1.2, 1.2, levels(cut(runif(100), cutPval)), col = col.legend, 
            pt.bg = col.legend, pch = rep(22, 4), horiz = TRUE, 
            cex = 1, bty = "n", pt.cex = 1.6, yjust = 0.5)
        max.y <- min(pval[pval[,2]>0, 2], na.rm = TRUE)
        if (ylim.sup > max.y) 
            ylim.sup <- max.y
        for (i in 1:n.chr) {
            snp <- gen.info[gen.info[, 2] == chr[i], 1]
            pos <- gen.info[gen.info[, 2] == chr[i], 3]
            dat <- pval[dimnames(pval)[[1]] %in% snp, -1]
            col <- as.character(cut(dat, cutPval, labels = col.legend))
            col[is.na(col)] <- "white"
            dat[is.na(dat)] <- 1
            plot(pos, -log(dat), axes = FALSE, xlab = "", ylab = "", 
                type = "h", col = col, xlim = limits, ylim = c(0, 
                  -log(ylim.sup)))
            segments(limits[1], 0, max(pos), 0)
            segments(centro[i], -1, centro[i], -log(max.y), lwd = 3, 
                col = "darkblue", xpd = TRUE)
            text(par("usr")[1], 0, chr[i], cex = 1, adj = 0.5)
        }
        plot(limits, c(1, 1), axes = FALSE, xlab = "", ylab = "", 
            type = "n", col = col, xlim = limits)
        text(limits[1], 1, limits[1], cex = 1, adj = 0)
        text(limits[2], 1, limits[2], cex = 1, adj = 1)
        text((limits[2] - limits[1])/2, 1, "Genomic Position", 
            cex = 1, adj = 0, font = 2)
    }
}






Bonferroni.sig<-
function (x, model = "codominant", alpha = 0.05, include.all.SNPs = FALSE) 
{
    if (!inherits(x, "WGassociation")) 
        stop("x must be a 'WGassociation' object")
    x <- attr(x, "pvalues")
    model.type <- names(x)
    m <- charmatch(model, model.type, nomatch = 0)
    if (m == 0) 
        stop("this model is was not fitted")
    temp1 <- grep("no", as.character(x[, m ]))
    temp2 <- c(1:nrow(x))[is.na(x[, m ])]
    temp <- c(temp1, temp2)
    if (!(include.all.SNPs)) {
        x <- x[-temp, c(1, (m ))]
        cut.p <- alpha/nrow(x)
    }
    else cut.p <- alpha/nrow(x)
    cat("number of tests: ", nrow(x), "\n")
    cat("alpha:", alpha, "\n")
    cat("corrected alpha:", cut.p, "\n")
    significant <- x[as.numeric(x[, 2]) <= cut.p, ]
    if (dim(significant)[1] == 0) {
        cat("   No significant SNPs after Bonferroni correction \n")
        ans <- NULL
    }
    else {
        ans <- significant
        print(as.matrix(ans), na.print = "-", quote = FALSE)
    }
    invisible(ans)
}




sortSNPs<- function (data, colSNPs, info) 
{
    o <- order(info[, 2], info[, 3])
    label.SNPs.o <- info[o, 1]
    label.SNPs <- names(data[, colSNPs])

#control
    ans <- match(label.SNPs, label.SNPs.o)
    if (sum(is.na(ans)) > 0) {
        cat("Warning: ")
        cat("the SNPs: ", as.character(label.SNPs[is.na(ans)]), 
            "\n  are not included in the file with the genomic positions and they are discarded \n")
    }

    ans <- match(label.SNPs.o, label.SNPs)


    out <- colSNPs[ans[!is.na(ans)]]
    out <- out[!is.na(out)]
    res <- list(pos=out, dataSorted=info[o,])
    res 
}






interactionPval<-function (formula, data, quantitative = is.quantitative(formula, 
    data), model="codominant") 
  {
 
    if(!inherits(data,"setupSNP"))
     stop("data must be an object of class 'setupSNP'")

    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m0 <- match(c("formula", "data"), names(mf), 0)
    mf <- mf[c(1, m0)]
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    temp0 <- as.character(mt)
    adj <- paste(temp0[2], temp0[1], temp0[3])

    fam <- ifelse(quantitative,"gaussian","binomial")   

    model.type <- c("codominant", "dominant", "recessive", 
                "overdominant","log-additive")
    m <- charmatch(model, model.type, nomatch = 0)
    if (m == 0) 
          stop("model must be codominant dominant recessive overdominant or log-additive")

    modelOK<-switch(m,codominant,dominant,recessive,overdominant,additive)


    colSNPs<-attr(data,"colSNPs")
    if (is.vector(colSNPs) & length(colSNPs) > 1) 
        dataSNPs.sel <- data[, colSNPs]
    else stop("data should have an attribute called 'colSNPs'. Try again 'setupsNP' function")
    
    dataSNPs <- data.frame(lapply(dataSNPs.sel,function(x,model.sel) model.sel(x),model.sel=modelOK)) 

    SNPs.label <- names(dataSNPs)

    dimnames(data)[[1]]<-1:nrow(data)

    i<-1
    n<-ncol(dataSNPs)

    pval<-matrix(NA,nrow=n,ncol=n)


    while (i<=n)
     {
         nas <- sum(!is.na(dataSNPs[, i]))
         n.nas <- length(dataSNPs[, i])

         if (is.Monomorphic(dataSNPs[, i]))
          {
            pval[i,]<-rep(NA,n)
          }

         else if (nas/n.nas<0.80)
          {
            pval[i,]<-rep(NA,n)
          }
         
         else if (length(table(dataSNPs[, i]))==1)
          {
            pval[i,]<-rep(NA,n)
          }

         else
         {
          j<-i+1
          while (j<=n)
           {
            nas <- sum(!is.na(dataSNPs[, j]))
            n.nas <- length(dataSNPs[, j])

            if (is.Monomorphic(dataSNPs[, j]))
             {
                pval[i,j]<-NA
             }
            else if (nas/n.nas<0.80)
             {
              pval[i,j]<-NA
             }

            else if (length(table(dataSNPs[, j]))==1)
             {
              pval[i,j]<-NA
             }


            else
             {
              mod.i <- glm(as.formula(paste(adj, "+ dataSNPs[, i]*dataSNPs[, j]")),
                        data = data, family=fam)
              subset <- 1:nrow(data) %in% as.numeric(rownames(mod.i$model))
              mod.a <- glm(as.formula(paste(adj, "+ dataSNPs[, i]+dataSNPs[, j]")),
                         data = data, family=fam, subset=subset)

              mod.b1 <- glm(as.formula(paste(adj, "+ dataSNPs[, i]")), data = data,
                         family=fam,subset=subset)
              mod.b2 <- glm(as.formula(paste(adj, "+ dataSNPs[, j]")), data = data,
                         family=fam,subset=subset)

              pval[i,j]<-anova(mod.a,mod.i,test="Chisq")$"P(>|Chi|)"[2]

              if(mod.b1$aic<=mod.b2$aic) 
                pval[j,i]<-anova(mod.b1,mod.a,test="Chisq")$"P(>|Chi|)"[2]
              else  
                pval[j,i]<-anova(mod.b2,mod.a,test="Chisq")$"P(>|Chi|)"[2]

             }
             j<-j+1
           }

          mod.0 <- glm(as.formula(paste(adj, "+ dataSNPs[, i]")), data = data,
                         family="gaussian")
          subset <- 1:nrow(data) %in% as.numeric(rownames(mod.0$model))
          mod.b <- glm(as.formula(paste(adj)), data = data,
                         family="gaussian",subset=subset)
          pval[i,i] <- anova(mod.b,mod.0,test="Chisq")$"P(>|Chi|)"[2]
         }

          i<-i+1
       }

 dimnames(pval)[[2]]<-SNPs.label
 dimnames(pval)[[1]]<-SNPs.label

 class(pval)<-"SNPinteraction" 
 attr(pval,"model") <- model.type[m]
 attr(pval,"gen.info")<-attr(data,"gen.info")
 pval
}



print.SNPinteraction<-function(x, ...)
{
 print.table(x,na.print="-")
}




plot.SNPinteraction<-function(x, main.tit, ...)
{
 
 control<-apply(x,1,function(x) sum(is.na(x))==length(x))
 x.OK<-x[!control,!control]
 
 if (!is.null(attr(x,"gen.info"))){
        genInfo<-attr(x,"gen.info")
        o <- order(genInfo[, 2], genInfo[, 3])
        label.SNPs <- as.character(genInfo[o, 1])
        label.SNPs <- label.SNPs[label.SNPs%in%dimnames(x.OK)[[1]]]
        orderSNPs.ok<-match(label.SNPs, dimnames(x.OK)[[1]])
        x.OK <- x.OK[orderSNPs.ok,orderSNPs.ok ]
        genInfo <- genInfo[genInfo[,1]%in%label.SNPs,]
    }
    else {
        label.SNPs <- dimnames(x.OK)[[1]]
    }
 

 old.xpd <- par("xpd")
 old.las <- par("las")
 old.mfrow <- par("mfrow")
 
 par(xpd=NA)
 
 m <- matrix(1:2, 1, 2)
 layout(m, widths=c(4.5, 1))

 on.exit(par(xpd = old.xpd, mfrow = old.mfrow, las = old.las))

# Other palettes:
# mypaletteOld<-brewer.pal(9,"Greens")
# mypaletteOld<-c("#F7FCF5", "#E5F5E0", "#C7E9C0","#A1D99B", "#74C476", "#41AB5D", "#238B45", "#006D2C", "#00441B")
# mypaletteOld<-brewer.pal(9,"Reds")
#  mypaletteOld<- c("#FFF5F0", "#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D", "#A50F15", "#67000D")
# This is used:  mypaletteOld<-brewer.pal(9,"YlGn")

 mypaletteOld <- c("#FFFFE5", "#F7FCB9", "#D9F0A3", "#ADDD8E", "#78C679", "#41AB5D", "#238443", "#006837", "#004529")

 mypalette<-mypaletteOld[c(9,6,4,3,3,2,2,1,1)]

 pvalCut<-c(0,0.001,0.01,0.05,0.1,0.2,0.3,0.5,0.7,1)

 image(1:nrow(x.OK),1:ncol(x.OK),x.OK,col=mypalette,breaks=pvalCut,
     axes=FALSE,xlab="",ylab="")
 
 axis(1,at=c(1:nrow(x.OK)),label=label.SNPs,las=3,cex.axis=0.7,col="darkgreen")
 axis(2,at=c(1:nrow(x.OK)),label=label.SNPs,las=1,cex.axis=0.7,col="darkgreen")
 

 if (missing(main.tit))
  main.tit<-paste("SNPs interactions --",attr(x,"model"),"model")

 title(main.tit,line=3)

 if (!is.null(attr(x,"gen.info"))) 
        n.snps <- table(genInfo[, 2])
 else n.snps <- nrow(x.OK)

 
 a <- c(0.5, cumsum(n.snps) + 0.5)

 b <- par("usr")
 segments(a, b[3], a, b[4] + diff(b[3:4]) * 0.02, col="darkblue",lwd=2)
 segments(b[3], a, b[4]+diff(b[3:4]) * 0.02, a, col="darkblue",lwd=2)

 abline(coef=c(0,1),xpd=FALSE,col="yellow")

 if(!is.null(attr(x,"gen.info")))
  {
   a <- par("usr")
   wh <- cumsum(c(0.5, n.snps))
   names.geno<-unique(genInfo[,2])
   n.gen<-length(names.geno)

   for (i in 1:n.gen)
    { 
      text(mean(wh[i + c(0, 1)]), a[4] + (a[4] - a[3]) * 0.025, names.geno[i],srt=45,cex=0.8,adj=0.2)
      text(a[4] + (a[4] - a[3]) * 0.025, mean(wh[i + c(0, 1)]), names.geno[i],srt=45,cex=0.8,adj=0.5)
    }
  }  

 
 image(0.5,1:10,matrix(pvalCut,nrow=1,ncol=10),col=rev(mypalette),breaks=pvalCut,axes=FALSE,
          xlab="",ylab="")
 marcas<-c(0.5,3.5,4.5,5.5,7.5,8.5,9.5,10.5)
 axis(2,marcas,rev(c(0,0.001,0.01,0.05,0.1,0.2,0.3,1)),pos=0.5)
 text(30,5.5,"pvalues",srt=90)

}


WGstats<-function(x,pSig=0.000001)
{

 if (!inherits(x, "WGassociation")) 
      stop("x must be an object of class 'WGassociation'")

 SNPs<-attr(x,"label.SNPs")
 genes<-attr(x,"gen.info")
 if (is.null(genes))
      stop("only implemented when chromosome information is available") 
 pvalues<-attr(x,"pvalues")
 nSNPs<-table(genes[,2])
 chr.l <- names(nSNPs)
 o<-orderChromosome(chr.l)
 chr <- chr.l[o]
 nSNPs.o<-nSNPs[o]

 info<-matrix(NA,nrow=length(chr),ncol=5)

 for (i in 1:length(chr))
 { 
  info0<-nSNPs.o[i]
  temp<-genes[genes[,2]==chr[i],] 
  aux<-pvalues[dimnames(pvalues)[[1]]%in%temp[,1],]
  info1<-round((table(aux[,1])/nrow(aux))*100,1)
  nSig<-sum(aux[,2]<=pSig,na.rm=TRUE)
  info2<-c(nSig,round((nSig/nrow(aux))*100,1))
  info[i,]<-c(info0,info1,info2)
 }

 ans<-data.frame(info)
 names(ans)<-c("SNPs (n)","Genot error (%)","Monomorphic (%)",
   "Significant* (n)","(%)")
 dimnames(ans)[[1]]<-chr
 
 print(ans)
 
 cat("\n *Number of statistically significant associations at level", pSig)
 cat("\n")
 invisible(ans)
}






haplo.interaction <- function(formula, data, SNPs.sel, quantitative = is.quantitative(formula, data), 
             haplo.freq.min=0.05,...)

{
   if (!inherits(data, "setupSNP")) 
        stop("data must be an object of class 'setupSNP'")
 
   control.SNPs<-sum(!is.na(match(names(data),SNPs.sel)))
   if (control.SNPs!=length(SNPs.sel))
     stop("Some of the SNPs selected are not in the data set")
    
   cl <- match.call()
   mf <- match.call(expand.dots = FALSE)
   m0 <- match(c("formula", "data", "subset"), names(mf), 0)
   mf <- mf[c(1, m0)]
   mf[[1]] <- as.name("model.frame")
   mf <- eval(mf, parent.frame())
   mt <- attr(mf, "terms")
   
   special <- c("int")
   Terms <- if (missing(data)) 
        terms(formula, special)
    else terms(formula, special, data = data)
   posInt <- attr(Terms, "specials")$int
   
   if(length(posInt))
    {
      var2<-mf[,posInt]
      if(!length(levels(var2))) 
       {
        stop("interaction variable must be a factor")
       }
    }  
   else 
      stop("formula needs an 'interaction' term")   
   
   control.missing<-dimnames(mf)[[1]]
   geno <- make.geno(data[dimnames(data)[[1]]%in%control.missing,], 
                   SNPs.sel)

   dep <- mf[, 1]
   if (ncol(mf) > 2)
    adj <- data.frame(mf[, -c(1,posInt)])
   else
    adj <- NULL 
   varAdj <- attr(mt, "term.labels")[-(posInt-1)]
   varInt <- attr(mt, "term.labels")[posInt-1]
   varInt <-gsub("int\\(","",varInt)
   varInt <-gsub("\\)","",varInt)

   out<-haplo.inter.fit(geno, var2, dep, adj , ifelse(quantitative,"gaussian","binomial"), haplo.freq.min, ...)
   
   res.corner<-out[[1]]
   xx<-dimnames(res.corner)[[2]] 
   xx[xx=="li"]<-"lower"
   xx[xx=="ls"]<-"upper"
   dimnames(res.corner)[[2]]<-xx

   temp<-out[[2]]   
   etiq1<-dimnames(temp)[[1]]
   aux0<-dimnames(temp)[[2]]
   etiq2<-aux0[seq(2,length(aux0),3)]
      
   ans<-list(NA)
   for (i in 1:length(etiq2))
       {
         ans[[i]]<-temp[,c(1,(2+3*(i-1)):(4+3*(i-1)))]
         ans[[i]][1,2]<-ifelse(quantitative,0,1)  
         ans[[i]]<-ans[[i]][,-1]  
         if (!quantitative)
           dimnames(ans[[i]])[[2]]<-c("OR","lower","upper")               
         else
            dimnames(ans[[i]])[[2]]<-c("diff","lower","upper")               
       }
   names(ans)<-etiq2      
   res.int1<-ans      

   temp<-out[[3]]   
   etiq1<-dimnames(temp)[[1]]
   aux0<-dimnames(temp)[[2]]
   etiq2<-aux0[seq(2,length(aux0),3)]
      
   ans2<-list(NA)
   for (i in 1:length(etiq1))
       {
         ans.i<- matrix(temp[i,][-1],nrow=length(etiq2) ,ncol=3,byrow=TRUE)
         ans2[[i]]<-data.frame(ans.i)
         dimnames(ans2[[i]])[[1]]<-etiq2 
         ans2[[i]][1,1]<-ifelse(quantitative,0,1)  
         if (!quantitative)
           dimnames(ans2[[i]])[[2]]<-c("OR","lower","upper")               
         else
           dimnames(ans2[[i]])[[2]]<-c("diff","lower","upper")               
       }
   names(ans2)<-etiq1      
   res.int2<-ans2      
  
   res<-list(res.corner,res.int1,res.int2,out$pval)
   attr(res,"label.snp")<-SNPs.sel
   attr(res,"varAdj")<-varAdj 
   attr(res,"varInt")<-varInt
   attr(res,"quantitative")<-quantitative
   class(res)<-"haploOut"
   res

}



haplo.inter.fit <- function(geno, var2, dep, adj = NULL, fam, haplo.freq.min, ...)
{
	param <- "var2"; #This is the name of the formal parameter that contains the interaction variable

	# haplo.inter(genotip, Datos$sexo, Datos$grupo)
	# haplo.inter(genotip, Datos$sexo, Datos$grupo, Datos[,c("rcal.dia", "n.edad")],0.005,1)

	# Funció que retorna les tres taules d'interaccions a partir dels coef i se de l'ajust d'un model haplo.glm

	# Li passo els seguents paràmetres per crear un model haplo.glm amb interaccions

	# La variable geno és del tipus: (genotip <- setupGeno(cbind(allele(g.81,1:2), allele(g.82,1:2), allele(g.83,1:2))) on g.81<-genotype(XRCC1.81)  )
	# var2 (Variable que interacciona amb l'haplotip)
	# dep (Variable depenent)
	# adj (Variables d'ajust)
	# fam (funció pel glm "binomial"/"gaussian")
	# haplo.freq.min (haplo.freq.min = 0.01) Frequència mínima haplotípica.
	# num.status (variable resposta categòrica o numèrica?)

	# Models
	if (is.null(adj))
	{
		mod <- haplo.glm(dep~geno*var2, family=fam, allele.lev=attributes(geno)$unique.alleles, control=haplo.glm.control(haplo.freq.min=haplo.freq.min,...))
		mod.b <- haplo.glm(dep~geno+var2, family=fam, allele.lev=attributes(geno)$unique.alleles, control=haplo.glm.control(haplo.freq.min=haplo.freq.min,...))
	}
	else
	{
		mod <- haplo.glm(dep~.+geno*var2, family=fam, allele.lev=attributes(geno)$unique.alleles, control=haplo.glm.control(haplo.freq.min=haplo.freq.min,...), data=adj)
		mod.b <- haplo.glm(dep~.+geno+var2, family=fam, allele.lev=attributes(geno)$unique.alleles, control=haplo.glm.control(haplo.freq.min=haplo.freq.min,...), data=adj)
	}

	# Calculo la p d'interacció entre l'haplotip i la variable
	pval.haplo <- 1-pchisq(mod$lrt$lrt-mod.b$lrt$lrt, mod$lrt$df-mod.b$lrt$df)

	#  Guardo els coeficients en el vector co
	co <- mod$coef
	noms.coef <- names(co)

	# Matriu de covariàncies
	# Selecciono les mateixes files i columnes que a la matriu de coef, ja que la matriu q retorna mod$var.mat no té dimnames i estan tots els haplotips sense agrupar els rare. 

	mat.cov <- mod$var.mat[1:length(mod$coef), 1:length(mod$coef)]
	dimnames(mat.cov)<- list(names(mod$coef), names(mod$coef))

	# Creo la nova matriu (mat.coef) amb els resultats dels coeficients, per files var2 i columnes geno, amb el format de taula final.
	mat.coef <- matrix(nrow=1+length(mod$haplo.names), ncol=length(levels(var2)))
	dimnames(mat.coef) <- list(c(mod$haplo.base, mod$haplo.names), levels(var2))

	# Creo una matriu amb les variàncies del mateix format que la de coeficients.
	mat.var <- mat.coef
	mat.coef[1,1] <- 0 # Referència

	# creo matrius per taules marginals de row i col
	mat.coef.c <- mat.coef
	mat.var.c <- mat.var
	mat.coef.r <- mat.coef
	mat.var.r <- mat.var
 
	# Inicialitzo i, j
	i <- 1
	j <- 1

	# Selecciono la posició dels efectes principals (de tot el conjunt elimino les interaccions)
	# Omplo la primera columna (correspon als efectes principals de la geno)
	for (i in 1:(dim(mat.coef)[1]-1))
	{
		pos <- rownames(mat.coef)[i+1];

		mat.coef[i+1,1] <- co[pos];
		mat.var[i+1,1] <- mat.cov[pos,pos];
		mat.coef.c[i+1,1] <- co[pos]
		mat.var.c[i+1,1] <- mat.cov[pos,pos]
		mat.coef.r[i+1,1] <- NA
		mat.var.r[i+1,1] <- NA
	}

	# Omplo la primera fila (correspon als efectes principals de la var2)
	for (j in 1:(dim(mat.coef)[2]-1))
	{
		pos <- paste(param,colnames(mat.coef)[j+1],sep="");

		mat.coef[1,j+1] <- co[pos];
		mat.var[1,j+1] <- mat.cov[pos,pos];

		mat.coef.c[1,j+1] <- NA
		mat.var.c[1,j+1] <- NA

		mat.coef.r[1,j+1] <- co[pos]
		mat.var.r[1,j+1] <- mat.cov[pos,pos]
	}

	# Omplo les interaccions

	for (j in 1:(dim(mat.coef)[2]-1))
	{
		for (i in 1:(dim(mat.coef)[1]-1))
		{

			pos.x <- rownames(mat.coef)[i+1];
			pos.y <- paste(param,colnames(mat.coef)[j+1],sep="");
			pos.xy <- paste(rownames(mat.coef)[i+1],":",param,colnames(mat.coef)[j+1],sep="");

			mat.coef[i+1,j+1] <- mat.coef[i+1,1] + mat.coef[1,j+1] + co[pos.xy]
			mat.var[i+1,j+1]  <- mat.var[i+1,1] + mat.var[1,j+1] + mat.cov[pos.xy, pos.xy] + 2*(mat.cov[pos.x, pos.y] + mat.cov[pos.x, pos.xy] + mat.cov[pos.y, pos.xy])

			mat.coef.c[i+1,j+1] <- mat.coef.c[i+1,1] + co[pos.xy]
			mat.var.c[i+1,j+1]  <- mat.var.c[i+1,1] +  mat.cov[pos.xy, pos.xy] + 2*mat.cov[pos.x, pos.xy]

			mat.coef.r[i+1,j+1] <- mat.coef.r[1,j+1] + co[pos.xy]
			mat.var.r[i+1,j+1]  <- mat.var.r[1,j+1] + mat.cov[pos.xy, pos.xy] + 2*mat.cov[pos.y, pos.xy]
		}
	}

	# Calculo la matriu d'ORs/coefs i ICs
	mat.or <- mat.coef
	mat.li <- mat.coef - (1.96 * sqrt(mat.var))
	mat.ls <- mat.coef + (1.96 * sqrt(mat.var))

	mat.or.c <- mat.coef.c
	mat.li.c <- mat.coef.c - (1.96 * sqrt(mat.var.c))
	mat.ls.c <- mat.coef.c + (1.96 * sqrt(mat.var.c))

	mat.or.r <- mat.coef.r
	mat.li.r <- mat.coef.r - (1.96 * sqrt(mat.var.r))
	mat.ls.r <- mat.coef.r + (1.96 * sqrt(mat.var.r))

	#En cas de regressió logística transformem els coeficients a OR's
	if (fam=="binomial")
	{
		mat.or <- exp(mat.or)
		mat.li <- exp(mat.li)
		mat.ls <- exp(mat.ls)

		mat.or.c <- exp(mat.or.c)
		mat.li.c <- exp(mat.li.c)
		mat.ls.c <- exp(mat.ls.c)

		mat.or.r <- exp(mat.coef.r)
		mat.li.r <- exp(mat.coef.r - 1.96 * sqrt(mat.var.r))
		mat.ls.r <- exp(mat.coef.r + 1.96 * sqrt(mat.var.r))
	}

	# Ordeno les tres matrius d'ORs i IC en una taula final (mat.fi)
	mat.fi <- NULL
	mat.fi.c <- NULL
	mat.fi.r <- NULL
	i <- 1
	j <- 1

	for(i in 1:length(levels(var2)))
	{
		mat.fi <- cbind(mat.fi, mat.or[,i], mat.li[,i], mat.ls[,i])
		dimnames(mat.fi)[[2]][j:(j+2)] <- c(dimnames(mat.or)[[2]][i], "li", "ls")

		mat.fi.c <- cbind(mat.fi.c, mat.or.c[,i], mat.li.c[,i], mat.ls.c[,i])
		dimnames(mat.fi.c)[[2]][j:(j+2)] <- c(dimnames(mat.or.c)[[2]][i], "li", "ls")
		mat.fi.r <- cbind(mat.fi.r, mat.or.r[,i], mat.li.r[,i], mat.ls.r[,i])
		dimnames(mat.fi.r)[[2]][j:(j+2)] <- c(dimnames(mat.or.r)[[2]][i], "li", "ls")

		j <- j + 3
	}

	# Afegeixo una primera columna a mat.fi amb la frequència dels haplotips
	HapFreq <- NA
	mat.fi <- cbind(HapFreq, mat.fi)
	mat.fi.c <- cbind(HapFreq, mat.fi.c)
	mat.fi.r <- cbind(HapFreq, mat.fi.r)

	# Arreglo les etiquetes de les files de la matriu corresponent als haplotips.
	# Elimino de les etiquetes el "geno.", de forma que em quedi només el número (geno.1 per 1),
	# per això selecciono a partir de la posició 6 i fins la 9 pensant que el màxim de llargada són els rare.
	# I selecciono a partir de la segona fila pq a la primera hi ha el base que ja és sempre un número.

	dimnames(mat.fi)[[1]][2:nrow(mat.fi)] <- substr(dimnames(mat.fi)[[1]][2:nrow(mat.fi)], 6,9)
	dimnames(mat.fi.c)[[1]][2:nrow(mat.fi.c)] <- substr(dimnames(mat.fi.c)[[1]][2:nrow(mat.fi.c)], 6,9)
	dimnames(mat.fi.r)[[1]][2:nrow(mat.fi.r)] <- substr(dimnames(mat.fi.r)[[1]][2:nrow(mat.fi.r)], 6,9)

	# Un cop tinc els dimnames amb números (3,4,...), li ajunto les etiquetes i freq. haplotípiques que li corresponen.
	# Ex: geno.2, ara és 2, i correspon a CGA, freq 0.34.
	i <- 1

	for(i in 1:nrow(mat.fi))
	{
		if (dimnames(mat.fi)[[1]][i]!="rare")
		{
			num <- as.numeric(dimnames(mat.fi)[[1]][i]) 
			mat.fi[i,c("HapFreq")] <- mod$haplo.freq[num]
			dimnames(mat.fi)[[1]][i] <- paste(mod$haplo.unique[num,], collapse="")

			mat.fi.c[i,c("HapFreq")] <- mod$haplo.freq[num]
			dimnames(mat.fi.c)[[1]][i] <- paste(mod$haplo.unique[num,], collapse="")
			mat.fi.r[i,c("HapFreq")] <- mod$haplo.freq[num]
			dimnames(mat.fi.r)[[1]][i] <- paste(mod$haplo.unique[num,], collapse="")
		}
 		else if (dimnames(mat.fi)[[1]][i]=="rare")
 		{
			mat.fi[i,c("HapFreq")] <- 1-sum(mat.fi[1:(nrow(mat.fi)-1),1])
			mat.fi.c[i,c("HapFreq")] <- 1-sum(mat.fi.c[1:(nrow(mat.fi.c)-1),1])
			mat.fi.r[i,c("HapFreq")] <- 1-sum(mat.fi.r[1:(nrow(mat.fi.r)-1),1])
		}
	}

	# Ordeno la taula final per frequència haplotípica (de major a menor), deixant els rare sempre al final.
	# I suposant que el haplobase sempre serà el més frequent i quedarà a la primera categoria ( en cas que no fos així no passaria res, només que 
	# la categoria de referència no es mostraria a la primera fila.

	mat.fi <- mat.fi[c(order(mat.fi[1:(nrow(mat.fi)-1),c("HapFreq")], decreasing=TRUE), nrow(mat.fi)),]
	mat.fi.c <- mat.fi.c[c(order(mat.fi.c[1:(nrow(mat.fi.c)-1),c("HapFreq")], decreasing=TRUE), nrow(mat.fi.c)),]
	mat.fi.r <- mat.fi.r[c(order(mat.fi.r[1:(nrow(mat.fi.r)-1),c("HapFreq")], decreasing=TRUE), nrow(mat.fi.r)),]

	if (fam=="binomial")
	{
		# Si hi ha algun valor molt gran a la taula final li assigno el valor Inf.
		mat.fi[mat.fi>999] <- Inf
		mat.fi.c[mat.fi.c>999] <- Inf
		mat.fi.r[mat.fi.r>999] <- Inf
	}

	# Arrodoneixo a 4 decimals la frequència haplotípica i la resta a 2.
	mat.fi[, c("HapFreq")] <- round(mat.fi[, c("HapFreq")], 4)
	mat.fi[, 2:ncol(mat.fi)] <- round(mat.fi[, 2:ncol(mat.fi)], 2)

	mat.fi.c[, c("HapFreq")] <- round(mat.fi.c[, c("HapFreq")], 4)
	mat.fi.c[, 2:ncol(mat.fi.c)] <- round(mat.fi.c[, 2:ncol(mat.fi.c)], 2)

	mat.fi.r[, c("HapFreq")] <- round(mat.fi.r[, c("HapFreq")], 4)
	mat.fi.r[, 2:ncol(mat.fi.r)] <- round(mat.fi.r[, 2:ncol(mat.fi.r)], 2)

	list(mat.fi=mat.fi, mat.fi.c=mat.fi.c, mat.fi.r=mat.fi.r, pval=pval.haplo) 
	
}


make.geno<- function (data, SNPs.sel) 
{
    if (!inherits(data, "setupSNP")) 
        stop("data must be an object of class 'setupSNP'")

    ans<-togeno(data[,SNPs.sel[1]],sep="/",lab=SNPs.sel[1])
    for (i in 2:length(SNPs.sel))
     {
      ans.i<-togeno(data[,SNPs.sel[i]],sep="/",lab=SNPs.sel[i])
      ans<-cbind(ans,ans.i)
     }
    geno<-setupGeno(ans)
    geno
}


togeno<-function(f,sep=sep,lab=lab)
{
  nam<-paste(lab,c("1","2"),sep=".")
  f<-as.character(factor(f))
  f[is.na(f)]<-paste("0",sep,"0",sep="")
  g<-as.data.frame(t(matrix(unlist(strsplit(f,sep)),2,length(f))))
  names(g)<-nam
  g
}


print.haploOut<-function(x, digits = max(3, getOption("digits") - 3), ...)
{
   cat("\n")
   cat("     Haplotype using SNPs:", attr(x, "label.snp"), " adjusted by:", attr(x, 
            "varAdj"), "\n")
   cat(" Interaction \n")
   cat("-------------------------\n")
   etiq<-dimnames(x[[1]])[[2]]

   if(attr(x,"quantitative"))
    {  
      etiq[2]<-paste(etiq[2],"(dif)")
      etiq[5]<-paste(etiq[5],"(dif)")
    }
   else
    {  
      etiq[2]<-paste(etiq[2],"(OR)")
      etiq[5]<-paste(etiq[5],"(OR)")
    }

   dimnames(x[[1]])[[2]]<-etiq
   print(x[[1]])
   cat("\n")
   cat("p interaction:",x[[4]],"\n")
   cat("\n",paste(attr(x,"varInt"),"within haplotype"), "\n")
   cat("-------------------------\n")
   print(x[[2]])
   cat(paste("haplotype within",attr(x,"varInt")), "\n")
   cat("-------------------------\n")
   print(x[[3]])

}



table.interaction <- function(var, dep, adj = NULL, int, num.status, level)
{       
	# taula.int(Datos$XRCC1.81, Datos$grupo, Datos[,c("sexo", "rcal.dia")], Datos$n.edad)
	# taula.int(Datos$XRCC1.81, Datos$grupo, NULL, Datos$n.edad)

	if (num.status==0) #Categorical response variable
	{
	    var <- as.factor(var)
	    dep <- as.factor(dep)

	    if (is.null(adj))
	    {
	       
	        m.t  <- glm(dep~ as.numeric(var) + int, family = binomial)
	
	        subset <- 1:length(var)%in%as.numeric(rownames(m.t$model));
	
	        m.b   <- glm(dep~ var + int,         subset = subset, family = binomial)
	        m.int <- glm(dep~ var/int,         subset = subset, family = binomial)
	        m.t.int <- glm(dep~ as.numeric(var) * int, subset = subset, family = binomial)
	
	    }
	    else
	    {
	        m.t  <- glm(dep~. + as.numeric(var) + int, family = binomial, data=adj)
	
	        subset <- 1:length(var)%in%as.numeric(rownames(m.t$model));
	
	        m.b   <- glm(dep~. + var + int,         subset = subset, family = binomial, data=adj)
	        m.int <- glm(dep~. + var/int,         subset = subset, family = binomial, data=adj)
	        m.t.int <- glm(dep~. + as.numeric(var) * int, subset = subset, family = binomial, data=adj)
	        
	    }
	       
		var.int <- factor(paste(levels(var)[var], levels(int)[int]), levels = outer(levels(var), levels(int), paste),
						  exclude = c(paste(levels(var), ""), paste("", levels(int)), paste(" ")))
	
		ta <- table(var.int[subset], dep[subset])

		# Matriu de coeficients i cov
	
		mat.coef <- merge(m.int$coef, summary(m.int)$coef, by=0, all.x=TRUE, sort=FALSE)
		nom.pos <- data.frame(names(m.int$coef), ordre=1:length(m.int$coef))
		mat.ordre <- merge(nom.pos, mat.coef, by.x=1, by.y=1, all.x=TRUE, sort=FALSE)
		mat.ordre <- mat.ordre[order(mat.ordre$ordre),]
	
		a <- as.matrix(mat.ordre[,c("Estimate")])
		se <- as.matrix(mat.ordre[,c("Std. Error")])
		mat <- cbind(a, se)
		selec <- dim(mat)[1.] - (length(levels(int)) - 1.) * length(levels(var))
		o <- (selec + 1.):dim(mat)[1.]
		k <- matrix(nrow = length(levels(var)), ncol = 3.)
		k[, 1.] <- 1.
		taula <- cbind(exp(mat[o, 1.]), exp(mat[o, 1.] - 1.96 * mat[o, 2.]), exp(mat[o, 1.] + 1.96 * mat[o, 2.]))
		taula[taula > 999.] <- NA
		ktaula <- rbind(k, round(taula, 2.))

		ktaula <- cbind(ta, ktaula)
	
		i <- 1;
		j <- 1;
		step <- length(levels(var));
		taula.int <- NULL;
		while (i <= nrow(ktaula))
	    {
			aux <- ktaula[i:(i+step-1),];
			colnames(aux)[3] <- levels(int)[j];
			taula.int <- cbind(taula.int, aux);
			i <- i + step;
			j <- j + 1;
		}

		#Check if interaction pvalues are NA
		pval <- anova(m.b, m.int, test = "Chi")$"P(>|Chi|)"[2];
		if (is.na(pval))
		{
			pval <- "NA";
		}
		else
		{
			pval <- format.pval(pval);
		}
		
		pval.trend <- anova(m.t, m.t.int, test = "Chi")$"P(>|Chi|)"[2];
		if (is.na(pval.trend))
		{
			pval.trend <- "NA";
		}
		else
		{
			pval.trend <- format.pval(pval.trend);
		}
		rownames(taula.int) <- levels(var);
		list(table=taula.int,pval=pval,trend=pval.trend);
	}
	else #Continuous response variable
	{
	    var <- as.factor(var)

	    if (is.null(adj))
	    {
	        m.t  <- glm(dep~as.numeric(var) + int, family = gaussian)

	        subset <- 1:length(var)%in%as.numeric(rownames(m.t$model));

	        m.b   <-   glm(dep~ var + int,             subset = subset, family = gaussian)
	        m.int <-   glm(dep~ var/int,               subset = subset, family = gaussian)
	        m.t.int <- glm(dep~ as.numeric(var) * int, subset = subset, family = gaussian)
	    }
	    else
	    {
	        m.t  <- glm(dep~. + as.numeric(var) + int, family = gaussian, data=adj)

	        subset <- 1:length(var)%in%as.numeric(rownames(m.t$model));

	        m.b <-     glm(dep~. + var + int,             subset = subset, family = gaussian, data=adj)
	        m.int <-   glm(dep~. + var/int,               subset = subset, family = gaussian, data=adj)
	        m.t.int <- glm(dep~. + as.numeric(var) * int, subset = subset, family = gaussian, data=adj)
	    }
		var.int <- factor(paste(levels(var)[var], levels(int)[int]), levels = outer(levels(var), levels(int), paste),
						  exclude = c(paste(levels(var), ""), paste("", levels(int)), paste(" ")))

		# Matriu de coeficients i cov

		mat.coef <- merge(m.int$coef, summary(m.int)$coef, by=0, all.x=TRUE, sort=FALSE)
		nom.pos <- data.frame(names(m.int$coef), ordre=1:length(m.int$coef))
		mat.ordre <- merge(nom.pos, mat.coef, by.x=1, by.y=1, all.x=TRUE, sort=FALSE)
		mat.ordre <- mat.ordre[order(mat.ordre$ordre),]
	
		a <- as.matrix(mat.ordre[,c("Estimate")])
		se <- as.matrix(mat.ordre[,c("Std. Error")])
		mat <- cbind(dif=a, lo=a-(1.96*se), up=a+(1.96*se))
		selec <- dim(mat)[1] - (length(levels(int)) - 1.) * length(levels(var))
		o <- (selec + 1):dim(mat)[1]
		mat <- mat[o,];
		
		i <- 1;
		while (i <= length(levels(var)))
		{
			mat <- rbind(c(0,NA,NA),mat);
			i <- i + 1;
		}
        
	    res <- cbind(Table.mean.se(var.int, dep, subset)$tp, mat);

	    i <- 1;
	    j <- 1;
	    step <- length(levels(var));
	    taula.int <- NULL;
	    while (i <= nrow(res))
	    {
	        aux <- res[i:(i+step-1),];
	        colnames(aux)[3] <- levels(int)[j];
	        taula.int <- cbind(taula.int, aux);
	        i <- i + step;
	        j <- j + 1;
	    }

		pval <- anova(m.b, m.int, test = "Chi")$"P(>|Chi|)"[2];
		if (is.na(pval))
		{
			pval <- "NA";
		}
		else
		{
			pval <- format.pval(pval);
		}
		
		pval.trend <- anova(m.t, m.t.int, test = "Chi")$"P(>|Chi|)"[2];
		if (is.na(pval.trend))
		{
			pval.trend <- "NA";
		}
		else
		{
			pval.trend <- format.pval(pval.trend);
		}
	    rownames(taula.int) <- levels(var);
		list(table=taula.int,pval=pval,trend=pval.trend);
	}
}

table.corner <- function (var, dep, adj = NULL, int = NULL, num.status, level) 
{


glm2<-function (formula, family = gaussian, data, weights, subset, 
    na.action, start = NULL, etastart, mustart, offset, control = glm.control(...), 
    model = TRUE, method = "glm.fit", x = FALSE, y = TRUE, contrasts = NULL, 
    ...) 
{
    call <- match.call()
    if (is.character(family)) 
        family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family)) 
        family <- family()
    if (is.null(family$family)) {
        print(family)
        stop("'family' not recognized")
    }
    if (missing(data)) 
        data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action", 
        "etastart", "mustart", "offset"), names(mf), 0)
    mf <- mf[c(1, m)]
#    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    switch(method, model.frame = return(mf), glm.fit = 1, glm.fit.null = 1, 
        stop("invalid 'method': ", method))
    mt <- attr(mf, "terms")
    Y <- model.response(mf, "numeric")
    if (length(dim(Y)) == 1) {
        nm <- rownames(Y)
        dim(Y) <- NULL
        if (!is.null(nm)) 
            names(Y) <- nm
    }
    X <- if (!is.empty.model(mt)) 
        model.matrix(mt, mf, contrasts)
    else matrix(, NROW(Y), 0)
    weights <- model.weights(mf)
    offset <- model.offset(mf)
    if (!is.null(weights) && any(weights < 0)) 
        stop("negative weights not allowed")
    if (!is.null(offset) && length(offset) != NROW(Y)) 
        stop(gettextf("number of offsets is %d should equal %d (number of observations)", 
            length(offset), NROW(Y)), domain = NA)
    mustart <- model.extract(mf, "mustart")
    etastart <- model.extract(mf, "etastart")
    fit <- glm.fit(x = X, y = Y, weights = weights, start = start, 
        etastart = etastart, mustart = mustart, offset = offset, 
        family = family, control = control, intercept = attr(mt, 
            "intercept") > 0)
    if (any(offset) && attr(mt, "intercept") > 0) {
        fit$null.deviance <- glm.fit(x = X[, "(Intercept)", drop = FALSE], 
            y = Y, weights = weights, offset = offset, family = family, 
            control = control, intercept = TRUE)$deviance
    }
    if (model) 
        fit$model <- mf
    fit$na.action <- attr(mf, "na.action")
    if (x) 
        fit$x <- X
    if (!y) 
        fit$y <- NULL
    fit <- c(fit, list(call = call, formula = formula, terms = mt, 
        data = data, offset = offset, control = control, method = method, 
        contrasts = attr(X, "contrasts"), xlevels = .getXlevels(mt, 
            mf)))
    class(fit) <- c("glm", "lm")
    fit
}




    if (num.status == 0) {
        var <- as.factor(var)
        dep <- as.factor(dep)
        var.int <- factor(paste(levels(var)[var], levels(int)[int]), 
            levels = outer(levels(var), levels(int), paste), 
            exclude = c(paste(levels(var), ""), paste("", levels(int)), 
                paste(" ")))
        if (is.null(adj)) {
            m.var.int <- glm2(dep ~ var.int, family = binomial)
            subset <- 1:length(var.int) %in% as.numeric(rownames(m.var.int$model))
            m.b <- glm2(dep ~ NULL, subset = subset, family = binomial)
        }
        else {
            m.var.int <- glm2(dep ~ . + var.int, family = binomial, 
                data = adj)
            subset <- 1:length(var.int) %in% as.numeric(rownames(m.var.int$model))
            m.b <- glm2(dep ~ ., subset = subset, family = binomial, 
                data = adj)
        }
        res <- cbind(table(var.int[subset], dep[subset]), intervals.or(m.var.int, 
            level, m.b, var.int)$or.ic)[, 1:5]
        i <- 1
        j <- 1
        step <- length(levels(var))
        taula.int <- NULL
        while (i <= nrow(res)) {
            aux <- res[i:(i + step - 1), ]
            colnames(aux)[3] <- levels(int)[j]
            taula.int <- cbind(taula.int, aux)
            i <- i + step
            j <- j + 1
        }
        rownames(taula.int) <- levels(var)
        colnames(taula.int)[1] <- colnames(taula.int)[3]
        colnames(taula.int)[2]<-c("")
        colnames(taula.int)[3]<-c("OR")
        colnames(taula.int)[6] <- colnames(taula.int)[8]
        colnames(taula.int)[7]<-c("")
        colnames(taula.int)[8]<-c("OR")

        taula.int
    }
    else {
        var <- as.factor(var)
        var.int <- factor(paste(levels(var)[var], levels(int)[int]), 
            levels = outer(levels(var), levels(int), paste), 
            exclude = c(paste(levels(var), ""), paste("", levels(int)), 
                paste(" ")))
        if (is.null(adj)) {
            m.var.int <- glm2(dep ~ var.int, family = gaussian)
            subset <- 1:length(var.int) %in% as.numeric(rownames(m.var.int$model))
            m.b <- glm2(dep ~ NULL, subset = subset, family = gaussian)
        }
        else {
            m.var.int <- glm2(dep ~ . + var.int, family = gaussian, 
                data = adj)
            subset <- 1:length(var.int) %in% as.numeric(rownames(m.var.int$model))
            m.b <- glm2(dep ~ ., subset = subset, family = gaussian, 
                data = adj)
        }
        res <- cbind(Table.mean.se(var.int, dep, subset)$tp, 
            intervals.dif(m.var.int, level, m.b, var.int, pval = FALSE)$m)
        i <- 1
        j <- 1
        step <- length(levels(var))
        taula.int <- NULL
        while (i <= nrow(res)) {
            aux <- res[i:(i + step - 1), ]
            colnames(aux)[3] <- levels(int)[j]
            taula.int <- cbind(taula.int, aux)
            i <- i + step
            j <- j + 1
        }
        rownames(taula.int) <- levels(var)
        colnames(taula.int)[2] <- colnames(taula.int)[3]
        colnames(taula.int)[c(1,3)]<-c("","")
        colnames(taula.int)[8] <- colnames(taula.int)[9]
        colnames(taula.int)[c(7,9)]<-c("","")

        taula.int
    }
}





intervals <-
function(o, level = .95, ...)
UseMethod("intervals")


intervals.haplo.glm <-
function (o, level = 0.95, sign = 1, FUN = exp, ...) 
{
    if (o$family$family != "binomial") 
        FUN = function(x) x
    z <- abs(qnorm((1 - level)/2))
    co <- summary(o)$coef
    control0 <- gsub("geno.", "", dimnames(summary(o)$coef)[[1]])
    control.geno<-grep("geno.", dimnames(summary(o)$coef)[[1]])
    control<-control0[c(1,control.geno)]
    n.control<-length(control)
    nombres <- rep(NA, n.control)
    freqs <- rep(NA, n.control)
    for (i in 1:n.control) {
        if (control[i] != "rare" & control[i] != "(Intercept)") {
            nombres[i] <- paste(o$haplo.unique[as.numeric(control[i]), 
                ], collapse = "")
            freqs[i] <- o$haplo.freq[as.numeric(control[i])]
        }
        else if (control[i] == "(Intercept)") {
            nombres[i] <- "(Intercept)"
            freqs[i] <- -1
        }
        else {
            nombres[i] <- "rare"
            freqs[i] <- sum(o$haplo.freq[o$haplo.rare])
        }
    }

    or <- FUN(co[, 1] * sign)
    li <- FUN(co[, 1] * sign - z * co[, 2])
    ls <- FUN(co[, 1] * sign + z * co[, 2])


# Add the reference haplotype (modified JRG 12-Nov-06
    if (o$family$family != "binomial") 
      or<-c(or[1],or[1],or[-1])
    else
      or<-c(or[1],1,or[-1])
    li<-c(li[1],NA,li[-1])
    ls<-c(ls[1],NA,ls[-1])
    pvals<-co[,4]
    pvals<-c(pvals[1],NA,pvals[-1])
    nombre.ref<-paste(o$haplo.unique[o$haplo.base,], collapse = "")
    nombre.cov<-dimnames(summary(o)$coef)[[1]][-c(1:n.control)]
    nombres<-c(nombres[1],nombre.ref,nombres[-1],nombre.cov)
    ncov<-length(nombre.cov)
    freqs<-c(freqs[1],o$haplo.freq[o$haplo.base],freqs[-1],rep(NA,ncov)) 
    names(freqs)<-names(or)
#
    
    r <- cbind(freqs, or, li, ls, pvals)

    if (o$family$family != "binomial") 
        dimnames(r) <- list(nombres, c("freq", "diff", paste(level * 
            100, "%", sep = ""), "C.I.", "      P-val"))
    else dimnames(r) <- list(nombres, c("freq", "or", paste(level * 
        100, "%", sep = ""), "C.I.", "      P-val"))

    class(r) <- "intervals"
    r
}


print.intervals <-
function (x, len = 6, d = 2, exclude.intercept = TRUE, pval = TRUE, 
    ...) 
{
    n <- x
    dd <- dim(n)
    mx <- 10^(len - (d + 1))
    n[n > mx] <- Inf
    a <- formatC(n, d, len, format = "f")
    dim(a) <- dd
    if (length(dd) == 1) {
        dd <- c(1, dd)
        dim(a) <- dd
        lab <- " "
    }
    else lab <- dimnames(n)[[1]]
    if (!pval) {
        mx <- max(nchar(lab)) + 1
        cat(paste(rep(" ", mx), collapse = ""), paste(" ", dimnames(n)[[2]]), 
            "\n")
        for (i in (1 + exclude.intercept):dd[1]) {
            lab[i] <- paste(c(rep(" ", mx - nchar(lab[i])), lab[i]), 
                collapse = "")
            if (i == (1 + exclude.intercept)) 
                cat(lab[i], formatC(n[i, 1], 4, 6, format = "f"), 
                  a[i, 2], "Reference haplotype", "\n")
            else cat(lab[i], ifelse(is.na(n[i, 1]),"      ",formatC(n[i, 1], 4, 6, format = "f")), 
                a[i, 2], "(", a[i, 3], "-", a[i, 4], ") \n")
        }
    }
    else {
        mx <- max(nchar(lab)) + 1
        cat(paste(rep(" ", mx), collapse = ""), paste(" ", dimnames(n)[[2]]), 
            "\n")
        for (i in (1 + exclude.intercept):dd[1]) {
            lab[i] <- paste(c(rep(" ", mx - nchar(lab[i])), lab[i]), 
                collapse = "")
            if (i == (1 + exclude.intercept)) 
                cat(lab[i], formatC(n[i, 1], 4, 6, format = "f"), 
                  a[i, 2], "Reference haplotype", "\n")
            else cat(lab[i], ifelse(is.na(n[i, 1]),"      ",formatC(n[i, 1], 4, 6, format = "f")), 
                a[i, 2], "(", a[i, 3], "-", a[i, 4], ") ", formatC(n[i, 
                  5], 4, 6, format = "f"), "\n")
        }
    }
}


summary.haplo.glm <- function(object, ...)
 {
   o <- object
   coe<-o$coe
   se<-sqrt(diag(o$var.mat)[1:length(coe)])
   z<-coe/se
   p<-2-2*pnorm(abs(z))
   list(coeficients=cbind(coe,se,z,p))
 }



tableHWE<-function (x, strata)
{
    if (!inherits(x, "setupSNP"))
        stop("x must be an object of class 'setupSNP'")
    colSNPs <- attr(x, "colSNPs")
    data.SNPs <- x[colSNPs]
    tt <- lapply(data.SNPs, table)
    ans <- cbind("HWE (p value)"=unlist(lapply(tt, SNPHWE)))
    if (!missing(strata)) {
        if (length(table(strata))>5) stop("strata looks numeric")
        strates <- names(table(strata) > 0)
        n.strata <- length(strates)
        i <- 1
        while (i <= n.strata) {
            data.SNPs.temp <- subset(data.SNPs, strata == strates[i])
            tt <- apply(data.SNPs.temp, 2, table)
            temp <- cbind("HWE (p value)" = unlist(lapply(tt, SNPHWE)))
            ans <- cbind(ans, temp)
            i <- i + 1
        }
        dimnames(ans)[[2]] <- c("all groups", strates)
    }
    class(ans)<-c("tableHWE","matrix")
    ans
}

print.tableHWE<-function(x, digits=4, sig=0.05, na="-", ...) 
 {
  x<-round(x,digits)
  x<-data.frame(x)
  if (ncol(x)<3)
   {
    names(x)[1]<-"HWE (p value)" 
    x$flag<- apply(x,1,function(x)ifelse(any(x<sig) & !is.na(x),"<-",""))
   }
  print(as.matrix(x), na.print=na, quote=FALSE, ...)
}



SNPHWE<-function(x)
{
 if (length(x)<3)
  {
   p<-NA
  }
 
 else
  {
   obs_hom1 <- x[1]
   obs_hets <- x[2]
   obs_hom2 <- x[3]

   if (obs_hom1 < 0 || obs_hom2 < 0 || obs_hets < 0)
      return(-1.0)

   # total number of genotypes
   N <- obs_hom1 + obs_hom2 + obs_hets
   
   # rare homozygotes, common homozygotes
   obs_homr <- min(obs_hom1, obs_hom2)
   obs_homc <- max(obs_hom1, obs_hom2)

   # number of rare allele copies
   rare  <- obs_homr * 2 + obs_hets

   # Initialize probability array
   probs <- rep(0, 1 + rare)

   # Find midpoint of the distribution
   mid <- floor(rare * ( 2 * N - rare) / (2 * N))
   if ( (mid %% 2) != (rare %% 2) ) mid <- mid + 1

   probs[mid + 1] <- 1.0
   mysum <- 1.0

   # Calculate probablities from midpoint down 
   curr_hets <- mid
   curr_homr <- (rare - mid) / 2
   curr_homc <- N - curr_hets - curr_homr

   while ( curr_hets >=  2)
      {
      probs[curr_hets - 1]  <- probs[curr_hets + 1] * curr_hets * (curr_hets - 1.0) / (4.0 * (curr_homr + 1.0)  * (curr_homc + 1.0))
      mysum <- mysum + probs[curr_hets - 1]

      # 2 fewer heterozygotes -> add 1 rare homozygote, 1 common homozygote
      curr_hets <- curr_hets - 2
      curr_homr <- curr_homr + 1
      curr_homc <- curr_homc + 1
      }    

   # Calculate probabilities from midpoint up
   curr_hets <- mid
   curr_homr <- (rare - mid) / 2
   curr_homc <- N - curr_hets - curr_homr
   
   while ( curr_hets <= rare - 2)
      {
      probs[curr_hets + 3] <- probs[curr_hets + 1] * 4.0 * curr_homr * curr_homc / ((curr_hets + 2.0) * (curr_hets + 1.0))
      mysum <- mysum + probs[curr_hets + 3]
         
      # add 2 heterozygotes -> subtract 1 rare homozygtote, 1 common homozygote
      curr_hets <- curr_hets + 2
      curr_homr <- curr_homr - 1
      curr_homc <- curr_homc - 1
      }    
 
    # P-value calculation
    target <- probs[obs_hets + 1]
    p <- min(1.0, sum(probs[probs <= target])/ mysum)

    #plo <- min(1.0, sum(probs[1:obs_hets + 1]) / mysum)

    #phi <- min(1.0, sum(probs[obs_hets + 1: rare + 1]) / mysum)

  }

 return(p)   
}



int<-function(x)
 {
  x
 }

orderChromosome<-function(x)
 {

   temp<-rep(NA,13)
   for (i in 10:22)
    {
     temp[i-9]<-grep(i,x)
    }

   temp2<-rep(NA,9)
   for (i in 1:9)
    {
      aux<-grep(i,x)
      temp2[i]<-aux[!aux%in%temp]
       
    }

   temp3<-grep("X",x)
   res<-c(temp2,temp,temp3)
   res
   

  }


####################################################################
####################################################################
#
#   WGassociation
#
####################################################################
####################################################################


scanWGassociation<-function (formula, data, model = c("all"), quantitative = is.quantitative(formula, 
    data), genotypingRate = 80) 
{
    if (!inherits(data, "setupSNP")) 
        stop("data must be an object of class 'setupSNP'")
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m0 <- match(c("formula", "data"), names(mf), 0)
    mf <- mf[c(1, m0)]
    if (length(grep("~", mf[[2]])) == 0) {
        formula <- as.formula(paste(mf[[2]], "~1", sep = ""))
        formula.1 <- list(formula)
        mode(formula.1) <- "call"
        mf[2] <- formula.1
    }
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    varDep<-mf[,1]
    control.missing<-dimnames(mf)[[1]]
    mt <- attr(mf, "terms")
    temp0 <- as.character(mt)

    ptm<-proc.time()

# Por ahora no
#    adj <- paste(temp0[2], temp0[1], temp0[3])
    if (temp0[3]!=1)
      stop("adjusted analysis is not yet implemented. Try 'WGassociation'")

    cat("Be patient. The program is computing ... \n")

    Terms <- if (missing(data)) 
        terms(formula)
    else terms(formula, data = data)
    ord <- attr(Terms, "order")
    if (any(ord > 1)) 
        stop("interaction term is not implemented")

    colSNPs <- attr(data, "colSNPs")
    if (is.vector(colSNPs) & length(colSNPs) > 1) 
        dataSNPs <- data[control.missing, colSNPs]
    else stop("data should have an attribute called 'colSNPs'. Try again 'setupSNP' function")

    type <- charmatch(model, c("codominant", "dominant", "recessive", 
        "overdominant", "log-additive", "all"))
    type <- sort(type)
    if (any(type %in% 6)) 
        type <- 1:6
    if (length(type)==0) 
        stop("model must be 'codominant','dominant','recessive','overdominant', \n                     'log-additive', 'all' or any combination of them")

    SNPs <- attr(data, "label.SNPs")

    out <-data.frame(pvalTest(dataSNPs, Y=varDep, quantitative=quantitative, type=type,
                      genotypingRate = genotypingRate ))
    lab.model <- c("codominant","dominant","recessive","overdominant","log-additive")
    if (max(type)==6)
     names(out)<-c("comments",lab.model)
    else
     names(out)<-c("comments",lab.model[type])
 
    for (i in 2:ncol(out)) out[, i] <- as.numeric(as.character(out[,i]))

    cost<-proc.time()-ptm
    cat("The program took", round(cost[3],2), "seconds \n")
     
    attr(out, "label.SNPs") <- attr(data, "label.SNPs")
    attr(out, "models") <- type
    attr(out, "quantitative") <- quantitative
    attr(out, "pvalues") <- out
    attr(out, "gen.info") <- attr(data, "gen.info")
    attr(out, "whole") <- attr(data, "whole")
    attr(out, "colSNPs") <- attr(data, "colSNPs")
    attr(out, "fast")<-TRUE
    class(out) <- c("WGassociation","data.frame")
    out
}


pvalTest<-function(dataX,Y,quantitative,type,genotypingRate)
 {
  pvalues<-t(data.frame(lapply(dataX,FUN=modelTest,Y=Y,
               quantitative=quantitative,type=type,
               genotypingRate = genotypingRate )))
  pvalues
 }


modelTest<-function(X,Y,quantitative,type,genotypingRate)
{
  control<-ifelse(6%in%type,5,length(type)) 
  controlGeno <- GenotypeRate(X)
  if (genotypingRate > controlGeno)
   {
    ans<-c("Genot error",rep(NA,control))
   }
  else
   {
    if (is.Monomorphic(X))
      ans<-c("Monomorphic",rep(NA,control))
    else { 
     ans<-NA  
     if (1%in%type | 6%in%type) { 
       mco<-assoc(Y,codominant(X),quantitative=quantitative) 
       ans<-c(ans,mco)
     }
     if (2%in%type | 6%in%type) {
       mdo<-assoc(Y,dominant(X),quantitative=quantitative) 
       ans<-c(ans,mdo)
     }
     if (3%in%type | 6%in%type) {
       mre<-assoc(Y,recessive(X),quantitative=quantitative) 
       ans<-c(ans,mre)
     }
     if (4%in%type | 6%in%type) {
       mov<-assoc(Y,overdominant(X),quantitative=quantitative) 
       ans<-c(ans,mov)
     }
     if (5%in%type | 6%in%type) {
      mad<-assoc(Y,additive(X),quantitative=quantitative) 
      ans<-c(ans,mad)
     }
    }
   }
 ans  
}



assoc<-function(y,x,test="lrt",quantitative)
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
 





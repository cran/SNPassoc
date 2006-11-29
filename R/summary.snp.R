`summary.snp` <-
function (object, print.out=TRUE, ...) 
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
        pvalueHWE <- SNPHWE(ans.g[, 1])
        dimnames(ans.a)[[1]] <- alle
        if (any(nas)) 
            ans.a <- rbind(ans.a, "NA's" = c(2 * sum(nas), NA))
        if (print.out) {
            cat("\n")
            cat("Alleles: \n")
            print(ans.a)
        }
     }
     else {
        pvalueHWE <- NA
        if (print.out) {
            cat("\n")
            cat("Alleles: \n")
            cat("   Monomorphic \n")
        }
     }
    
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


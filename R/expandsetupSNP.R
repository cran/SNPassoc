`expandsetupSNP` <-
function (o) 
{
    x <- summary(o, print.out = FALSE)
    allele <- rbind(x$allele.names)
    alleles <- paste(allele, collapse = "/")
    if (length(allele) == 1) 
        alleles <- c(allele)
    aux<-ifelse(any(!is.na(x$allele.freq[, 2])),
             max(x$allele.freq[, 2], na.rm = TRUE),NA)
    out <- data.frame(alleles = alleles, 
         major.allele.freq = ifelse(is.na(aux),NA,round(aux, 1)), 
         HWE = round(x$HWE, 6), missing = round(x$missing * 100, 1))
    out
}


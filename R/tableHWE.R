`tableHWE` <-
function (x, strata)
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


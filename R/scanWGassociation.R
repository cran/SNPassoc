`scanWGassociation` <-
function (formula, data, model = c("all"), quantitative = is.quantitative(formula, 
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


approx.bf.estimates.ave <- function (z, V, type, suffix=NULL, sdY=1) {
  listVec <- list(lABF.fn(z, V, sd.prior=sqrt(0.01)), lABF.fn(z, V, sd.prior=sqrt(0.1)), lABF.fn(z, V, sd.prior=sqrt(0.5)))
  m <- do.call(cbind, listVec)
  #lABF <- rowMeans(m)
  #lABF <- apply(m, 1, logsum-log(3))
  lABF <- apply(m, 1, function(x) logsum(x) -log(3))
  #ret <- data.frame(V, z, r, lABF)
  ret <- data.frame(V, z, lABF)
  if(!is.null(suffix))
    colnames(ret) <- paste(colnames(ret), suffix, sep = ".")
  return(ret)
}

lABF.fn <- function (z, V, sd.prior=0.15) {
  r <- sd.prior^2/(sd.prior^2 + V)
  lABF = 0.5 * (log(1 - r) + (r * z^2))
  return(lABF)
  #return(list("lABF" = lABF, "r" = r))
}

# changed the process.dataset in the coloc package
process.dataset <- function(d, suffix) {
  message('Processing dataset')

  nd <- names(d)
  if (! 'type' %in% nd)
    stop('The variable type must be set, otherwise the Bayes factors cannot be computed')

  if("beta" %in% nd && "varbeta" %in% nd && ("MAF" %in% nd || "sdY" %in% nd)) {
    if(length(d$beta) != length(d$varbeta))
      stop("Length of the beta vectors and variance vectors must match")
    if(!("snp" %in% nd))
      d$snp <- sprintf("SNP.%s",1:length(d$beta))
    if(length(d$snp) != length(d$beta))
      stop("Length of snp names and beta vectors must match")

    if(d$type == 'quant' & !('sdY' %in% nd))
      d$sdY <- sdY.est(d$varbeta, d$MAF, d$N)

    # is the BETA a log OR or OR/Effect size (not logged?)
    # if there are negative value, then it is a logOR?
    if (length(d$beta[d$beta<0])>0) log=TRUE  else log=FALSE
    #df <- approx.bf.estimates(z=d$beta/sqrt(d$varbeta),
    #                          V=d$varbeta, type=d$type, suffix=suffix, sdY=d$sdY)
    df <- approx.bf.estimates.ave(z=d$beta/sqrt(d$varbeta),
                              V=d$varbeta, type=d$type, suffix=suffix, sdY=d$sdY)
    if (type=="cc" & !log)  {
    #df <- approx.bf.estimates(z=log(d$beta)/sqrt(d$varbeta),
    #                          V=d$varbeta, type=d$type, suffix=suffix, sdY=d$sdY)
    df <- approx.bf.estimates.ave(z=log(d$beta)/sqrt(d$varbeta),
                              V=d$varbeta, type=d$type, suffix=suffix, sdY=d$sdY)
    }
    df$snp <- as.character(d$snp)
    return(df)
  }
    if("pvalues" %in% nd & "MAF" %in% nd & "N" %in% nd) {
    if (length(d$pvalues) != length(d$MAF))
      stop('Length of the P-value vectors and MAF vector must match')
    if(d$type=="cc" & !("s" %in% nd))
      stop("Must specify s if type=='cc' and you want to use approximate Bayes Factors")
    if(!("snp" %in% nd))
      d$snp <- sprintf("SNP.%s",1:length(d$pvalues))
    df <- data.frame(pvalues = d$pvalues,
                     MAF = d$MAF,
                     snp=as.character(d$snp))
    colnames(df)[-3] <- paste(colnames(df)[-3], suffix, sep=".")
    df <- subset(df, df$MAF>0 & df$pvalues>0) # all p values and MAF > 0
    abf <- approx.bf.p(p=df$pvalues, f=df$MAF, type=d$type, N=d$N, s=d$s, suffix=suffix)
    df <- cbind(df, abf)
    return(df)
  }

  stop("Must give, as a minimum, either (beta, varbeta, type) or (pvalues, MAF, N, type)")
}



# Apply this function to every row of data (every row is a locus)
combine.abf.locus <- function(l0, l1, l2, l3, l4, a0, a1, a2, a3, a4) {

  lH0.abf  <- log(a0) + l0
  lH1.abf  <- log(a1) + l1
  lH2.abf  <- log(a2) + l2
  lH3.abf  <- log(a3) + l3
  lH4.abf  <- log(a4) + l4

  all.abf <- c(lH0.abf, lH1.abf, lH2.abf, lH3.abf, lH4.abf)
  my.denom.log.abf <- logsum(all.abf)
  pp.abf <- exp(all.abf - my.denom.log.abf)
  names(pp.abf) <- paste("PP.H", (1:length(pp.abf)) - 1, ".abf", sep = "")
  print(signif(pp.abf,3))
  print(paste("PP abf for shared variant: ", signif(pp.abf["PP.H4.abf"],3)*100 , '%', sep=''))
  return(pp.abf)
}

fn.pw.gwas = function(p, data) {
  a0 = p[1]
  a1 = p[2]
  a2 = p[3]
  a3 = p[4]
  a4 = p[5]
  #print(nrow(data))
  suma = sum(exp(c(a0,a1,a2,a3,a4)))
  #print(paste("Alphas:" , exp(a0)/suma,exp(a1)/suma,exp(a2)/suma,exp(a3)/suma,exp(a4)/suma), sep=" ")
  lkl.frame.temp <- as.matrix(data)
  lkl.frame.temp[,1] <- log(exp(a0)/suma)
  lkl.frame.temp[,2] <- lkl.frame.temp[,2] + log(exp(a1)/suma)
  lkl.frame.temp[,3] <- lkl.frame.temp[,3] + log(exp(a2)/suma)
  lkl.frame.temp[,4] <- lkl.frame.temp[,4] + log(exp(a3)/suma)
  lkl.frame.temp[,5] <- lkl.frame.temp[,5] + log(exp(a4)/suma)
  #print(lkl.frame.temp[1,])
  #print(apply(lkl.frame.temp, MAR = 1, FUN = logsum))
  #print(log(sum(exp(lkl.frame.temp[1,]))))
  sumlkl = sum(apply(lkl.frame.temp, MAR = 1, FUN = logsum))
  #print(sumlkl)
  return(sumlkl)
}


##' Estimate the variance of Y in a processed dataset 
##'
##' Internal function
##' @title Var.data
##' @param f minor allele freq
##' @param N sample number
##' @return variance of MLE beta
##' @author Claudia Giambartolomei

estimate_vary = function(data_set){
  vars  =data_set$var *(data_set$n) * data_set$se ^2 * (data_set$n -1)  +  data_set$var * (data_set$n) * data_set$b^2
  print(paste("Var y = ",median(vars/(data_set$n-1))))
  return(median(vars/(data_set$n-1)))
}

############################################################################################################
# This program calculates the variance of a mean when the units that are included in it are chosen at random 
# from a finite population of $1, ..., N$ units. 
#
# The program is used to calculate the expresion:  \frac{1}{m} \cdot \sum_{i=1}^N Y_i(0) \cdot T_i
# The program is written such that X is equal to Y(0), however we can equally define X to be Y(1) and then we will need to switch "tr" with "1-tr"
############################################################################################################

meanBernulli = function(sigma0,mean0,n, p){
  
  stopifnot(  length(c(sigma0,mean0,n, p))==4  )
  
  var1 =  n * p * (1-p) * ( sigma0^2 + mean0^2  )
  var2 =  n * p * (1-p) 
  covariance = p * (1-p) * n * mean0
  Sigma = matrix(c(var1, covariance, covariance, var2), byrow=TRUE, ncol=2 )
  A = matrix(c( 1/(n*p) , - p * n * mean0/(n*p)^2  ), byrow=TRUE, ncol=1 )
  
  variance.of.mean = t(A) %*% Sigma %*% A
  return(variance.of.mean)
}









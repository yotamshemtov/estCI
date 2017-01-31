#' Confidence intervals for SATT and SATC
#' @return The function returns XXX 
# \itemize{
#    \item cohort
#    \item min.cohort
#    \item max.cohort
#    \item name.first - the first name that was inserted to the function.
#    \item percent.female - the percent of females with this first name.
#    \item frequency - the count of times the name appeared in this cohort year.   
# }
#' @examples
#' 
#' ### Example 1:
#' # Estimate the average treatment effect on the treated
#' y0 = rnorm(n,mean=10,sd=1)
#' y1 = rnorm(n,mean=13,sd=3)
#' tau = y1-y0
#' tr = rep(0,n)
#' tr.index = sample(c(1:n),size=n/2,replace=FALSE)
#' tr[tr.index]=1
#' y = tr*y1 + y0*(1-tr)
#' results = aveCI(y,tr, print=TRUE)
#' 
#' ### Example 2:
#' # Real data example using data from Tunca and Egeli (1996) that was also used by Rosenbaum (2001)
#' data(tunca.and.egeli.1996)
#' head(tunca.and.egeli.1996)
#' results = aveCI(outcome=tunca.and.egeli.1996$y,treatment=tunca.and.egeli.1996$tr, print=TRUE)
#' 
#' ### Example 3:
#' # Real data examle using data from Lalonde (1986)
#' data(lalonde1986)
#' head(lalonde1986)
#' y = lalonde1986$re78
#' tr= lalonde1986$treat
#' # A confidence interval at the 5 level - default.
#' results = aveCI(y,tr, print=TRUE, Sharp.CI = TRUE)
#' # A confidence interval at the 1 level.
#' results = aveCI(y,tr, print=TRUE, alpha=0.01, Sharp.CI = TRUE)
#' 
#' 
#' @param outcome  Outcome of interest \cr \cr
#' @param treatment  Treatmnet indicator that can take only two levels 0 or 1. \cr \cr
#' @param alpha The Type-I error rate (size) of the confidence interval.
#' @param adjust Indicator whether to adjust the outcome with respect to the design matrix X.
#' @param X the design matrix of pre-treatment variables that the outcome should be adjusted prior to estimating the CI. 

#' @description The funciton calculates confidence intervals (CI) for the sample average treatment effect on the treated (SATT), the sample average treatment effect on the controls (SATC), and the average treatment effect (SATE).
# The standard CI for the sample average treatment effect is based on Neyman's variance estimator,   
#' \deqn{ \text{Var}_{Neyman} = \left(  \frac{s}{d} \right) }
#' 
#' @export 

aveCI = function(outcome,treatment,X=NULL,alpha=0.05, print=TRUE, Sharp.CI=TRUE, Bernulli = FALSE){
  
  tr=as.numeric(treatment)
  y=outcome
  n = length(y)
  m = sum(tr)
  
  if(length(unique(tr))!=2){
    stop("Treatment variable has more than 2 levels - only two levels of treatment are allowed")
  }
  
  if(is.numeric(y)==FALSE){
    stop("Only numeric outcomes are allowed")
  }
  
  # pre-treatment covariate adjustment
  #if (adjust==TRUE){
  # y = lm(y~(.),data=as.data.frame(X))$res
  #}
  
  # The difference in means
  difference.means = mean(y[tr==1],na.rm=TRUE) - mean(y[tr==0],na.rm=TRUE)
  
  # Variance in control and treatment:
  var1 = var(y[tr==1],na.rm=TRUE)
  var0 = var(y[tr==0],na.rm=TRUE)
  var.pool = var(y,na.rm=TRUE)
  
  # Confidence interval based on Neyman's variance estimator:
  sd.neyman = sqrt(var1/m + var0/(n-m))
  ci.upper.neyman = difference.means + qnorm(1-alpha/2) * sd.neyman
  ci.lower.neyman = difference.means - qnorm(1-alpha/2) * sd.neyman
  #cat("Neyman's CI:  [",ci.lower.neyman, ci.upper.neyman,"]","\n")
  
  # Confidence interval based on Neyman's variance estimator with pooled variance:
  sd.pool = sqrt(var.pool/m + var.pool/(n-m))
  ci.upper.pool = difference.means + qnorm(1-alpha/2) * sd.pool
  ci.lower.pool = difference.means - qnorm(1-alpha/2) * sd.pool
  
  # Confidence interval using control units only:
  k.n.m = ( n/(n-m) )^2 * ( (n-m)/( (n-1)*m ) )
  sd.adj.control = sqrt( k.n.m * var0 ) 
  ci.upper.control = difference.means + qnorm(1-alpha/2) * sd.adj.control
  ci.lower.control = difference.means - qnorm(1-alpha/2) * sd.adj.control
  #cat("Our CI:  [",ci.low, ci.upper,"]","\n")
  
  # Confidence interval using treated units only:
  sd.adj.treated = sqrt( k.n.m * var1 ) 
  ci.upper.treated = difference.means + qnorm(1-alpha/2) * sd.adj.treated
  ci.lower.treated = difference.means - qnorm(1-alpha/2) * sd.adj.treated
  
  # CI using corr=1 
  sd.corr1 = sqrt(  var1/m + var0/(n-m) - ( sqrt(var1) - sqrt(var0) )^2/n  )
  ci.upper.corr1 = difference.means + qnorm(1-alpha/2) * sd.corr1
  ci.lower.corr1 = difference.means - qnorm(1-alpha/2) * sd.corr1
  
  
  ### Shortest CI:
  length.satc = ci.upper.treated-ci.lower.treated
  length.satt = ci.upper.control-ci.lower.control
  length.sate = ci.upper.neyman-ci.lower.neyman
  #length.sharp = ci.upper.sharp - ci.upper.sharp
  
  length.vec = c(length.sate,length.satc,length.satt)
  parameter.vec = c("SATE","SATC","SATT")
  
  length.shortestCI = length.vec[which.min(length.vec)]
  parameter.shortestCI = parameter.vec[which.min(length.vec)]
  
  if (parameter.shortestCI=="SATE"){
    ci.upper.shortestCI = ci.upper.neyman
    ci.lower.shortestCI = ci.lower.neyman
  } else if (parameter.shortestCI=="SATC"){
    ci.upper.shortestCI = ci.upper.treated
    ci.lower.shortestCI = ci.lower.treated
  } else if (parameter.shortestCI=="SATT"){
    ci.upper.shortestCI = ci.upper.control
    ci.lower.shortestCI = ci.lower.control
  } 
  
  length.gain.shortestCI = 1 - (ci.upper.shortestCI-ci.lower.shortestCI)/(ci.upper.neyman - ci.lower.neyman)
  
  
  if (Bernulli==TRUE){
    mean0 = mean(y[tr==0])
    mean1 = mean(y[tr==1])
    
    # Prediction interval for SATT
    var.mean.controls.under.treatment = meanBernulli(sigma0 = sqrt(var0), mean0 = mean0, n=n, p = m/n)
    var.difference.in.means = meanBernulli(sigma0 = sqrt(var1), mean0 = mean1, n=n, p = m/n) + meanBernulli(sigma0 = sqrt(var0), mean0 = mean0, n=n, p = 1-m/n)
    
    # Variance of the adjusted difference in means
    sd.adj.control.bernulli = sqrt( n/(n-m) * var.mean.controls.under.treatment ) 
    ci.upper.control.bernulli = difference.means + qnorm(1-alpha/2) * sd.adj.control.bernulli
    ci.lower.control.bernulli = difference.means - qnorm(1-alpha/2) * sd.adj.control.bernulli
    
    # Variance of the difference in means under Bernulli assignment mechanism
    sd.adj.diff.bernulli = sqrt(var.difference.in.means) 
    ci.upper.diff.bernulli = difference.means + qnorm(1-alpha/2) * sd.adj.diff.bernulli
    ci.lower.diff.bernulli = difference.means - qnorm(1-alpha/2) * sd.adj.diff.bernulli
  }
  
  ### results to return
  results = list(
    
    neymanCI = list(ci.lower.neyman = ci.lower.neyman,
                    ci.upper.neyman = ci.upper.neyman
                    ),
    
    sattCI = list(ci.lower.satt = ci.lower.control,
                  ci.upper.satt = ci.upper.control),
    
    satcCI = list(ci.lower.satc = ci.lower.treated,
                  ci.upper.satc = ci.upper.treated),
    
    shortestCI = list(ci.lower = ci.lower.shortestCI,
                      ci.upper = ci.upper.shortestCI,
                      parameter = parameter.shortestCI,
                      length.gain.shortestCI = length.gain.shortestCI)
    )
  
  if (Sharp.CI==TRUE){
    # CI sharp 
    sd.sharp = sqrt( sharpVar(y[tr==1],y[tr==0])  ) 
    ci.upper.sharp = round(difference.means + qnorm(1-alpha/2) * sd.sharp,dig=2)
    ci.lower.sharp = round(difference.means - qnorm(1-alpha/2) * sd.sharp,dig=2)
    
    length.gain.satt = 1 - (ci.upper.control - ci.lower.control)/(ci.upper.sharp - ci.lower.sharp)
    length.gain.satc = 1 - (ci.upper.treated - ci.lower.treated)/(ci.upper.sharp - ci.lower.sharp)
    
    results$sharpCI = list(ci.lower.sharp = ci.lower.sharp,
                   ci.upper.sharp = ci.upper.sharp,
                   length.gain.satt = length.gain.satt*100,
                   length.gain.satc = length.gain.satc*100)
  }
  
  if (Bernulli){
    
    results$Bernulli.SATT = list(ci.lower = ci.lower.control.bernulli,
                         ci.upper = ci.upper.control.bernulli,
                         length = ci.upper.control.bernulli - ci.lower.control.bernulli)
    
    results$Bernulli.diff.in.means = list(ci.lower = ci.lower.diff.bernulli,
                                  ci.upper = ci.upper.diff.bernulli,
                                  length = ci.upper.diff.bernulli - ci.lower.diff.bernulli)
    
    results$Bernulli.gain.SATT = 1 - (  ci.upper.control.bernulli - ci.lower.control.bernulli  )/( ci.upper.diff.bernulli - ci.lower.diff.bernulli )
  }
  
  
  # Print results
  if(print==TRUE){
    cat(
      " SATC - ", (1-alpha)*100," percent confidence interval: \n", 
        "      [",ci.lower.treated,", ",ci.upper.treated,"]",
      "\n \n",
      " SATT - ", (1-alpha)*100," percent confidence interval: \n", 
      "      [",ci.lower.control,", ",ci.upper.control,"]",
      "\n \n",
      " SATE (Neyman's CI) - ", (1-alpha)*100," percent confidence interval: \n", 
      "      [",ci.lower.neyman,", ",ci.upper.neyman,"]",
      "\n \n",
      "Shortest CI is for the parameter:  ",parameter.shortestCI,
      "\n",
      "The Shortest CI is smaller than Neyman's CI by:  ",round(length.gain.shortestCI,dig=2)*100,"%",
      "\n \n",
      " Variance of control units: ",round(var(y[tr==0],na.rm=TRUE),dig=2),
      "\n",
      " Variance of treated units: ",round(var(y[tr==1],na.rm=TRUE),dig=2),
      "\n",
      " F Test to Compare Two Variances:  P-value = ",var.test(y[tr==1],y[tr==0])$p.value,
      "\n \n \n",
      "Additional confidence intervals: \n \n",
      " SATE (Corr=1 CI) - ", (1-alpha)*100," percent confidence interval: \n", 
      "      [",ci.lower.corr1,", ",ci.upper.corr1,"]", 
      "\n",
        sep="")
    
    if (Sharp.CI==TRUE){
      cat(
      " SATE (Sharp CI) - ", (1-alpha)*100," percent confidence interval: \n", 
      "      [",ci.lower.sharp,", ",ci.upper.sharp,"]",
      "\n \n",
      sep="")
    }
    
  }
  
  # return results
  return(results)
}































#' Confidence intervals for SATT and SATC
#' @return The function returns TO-ADD help manual description 
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
#' n=1000
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

#' @description The funciton calculates confidence/prediction intervals (CI and PI) for sample average treatment effects. It calculates a prediction intervals for the sample average treatment effect on the treated (SATT), the sample average treatment effect on the controls (SATC), and a confidence intervals for the average treatment effect (SATE).
#' 
#' @export 

aveCI = function(outcome=NULL,treatment=NULL,alpha=0.05, print=TRUE, Sharp.CI=TRUE, variance.test=TRUE, stats.only=FALSE,
                 var1=NULL, var0=NULL, mu1=NULL, mu0=NULL, n=NULL, m=NULL, var.pool=NULL, rho = NULL){
  
  if(stats.only==TRUE){
    
    tmp = is.null(var1) | is.null(var0) |  is.null(mu0) | is.null(mu1) | is.null(m) |  is.null(n)
    if(tmp){
      stop("Error: One of the following moments is missing: var1, var0, mu1, mu0, n, or m.
           It is necesary to supply all the above moments when the stats.only option is equal TRUE")
    }
    
    variance.test=FALSE
    p = m/n
    k.n.m = 1/(p*(1-p)*n)
    
  }
  
  if(stats.only==FALSE){
    tr=as.numeric(treatment)
    y=outcome
    
    # General parameters and expressions
    n = length(y)
    m = sum(tr)
    p = m/n
    k.n.m = 1/(p*(1-p)*n)
    
    if(length(unique(tr))!=2){
      stop("Treatment variable has more than 2 levels - only two levels of treatment are allowed")
    }
    
    if(is.numeric(y)==FALSE){
      stop("Only numeric outcomes are allowed")
    }
    
    # Variance in control and treatment:
    var1 = var(y[tr==1],na.rm=TRUE)
    var0 = var(y[tr==0],na.rm=TRUE)
    #var.pool = var(y,na.rm=TRUE)
    
    # Means in treatment and control
    mu1 = mean(y[tr==1],na.rm=TRUE)
    mu0 = mean(y[tr==0],na.rm=TRUE)
  }
  
  # The difference in means
  difference.means =  mu1 - mu0
  
  # Confidence interval based on Neyman's variance estimator:
  sd.neyman = sqrt(var1/m + var0/(n-m))
  ci.upper.neyman = difference.means + qnorm(1-alpha/2) * sd.neyman
  ci.lower.neyman = difference.means - qnorm(1-alpha/2) * sd.neyman
  
  # Confidence interval based on Neyman's variance estimator with pooled variance:
  if( !is.null(var.pool) ){
    sd.pool = sqrt(var.pool/m + var.pool/(n-m))
    ci.upper.pool = difference.means + qnorm(1-alpha/2) * sd.pool
    ci.lower.pool = difference.means - qnorm(1-alpha/2) * sd.pool
  }

  
  ### Prediction interval for SATT
  sd.adj.control = sqrt( k.n.m * var0 ) 
  ci.upper.control = difference.means + qnorm(1-alpha/2) * sd.adj.control
  ci.lower.control = difference.means - qnorm(1-alpha/2) * sd.adj.control
  #cat("Our CI:  [",ci.low, ci.upper,"]","\n")
  
  ### Prediction interval for SATC
  sd.adj.treated = sqrt( k.n.m * var1 ) 
  ci.upper.treated = difference.means + qnorm(1-alpha/2) * sd.adj.treated
  ci.lower.treated = difference.means - qnorm(1-alpha/2) * sd.adj.treated
  
  # CI using corr=1 
  sd.corr1 = sqrt(  var1/m + var0/(n-m) - ( sqrt(var1) - sqrt(var0) )^2/n  )
  ci.upper.corr1 = difference.means + qnorm(1-alpha/2) * sd.corr1
  ci.lower.corr1 = difference.means - qnorm(1-alpha/2) * sd.corr1
  
  # CI for SATE given rho
  if(!is.null(rho)){
    sd.rho = sqrt( k.n.m * ( p^2*var0 + (1-p)^2*var1 + 2*p*(1-p)*rho*sqrt(var0)*sqrt(var1) ) ) 
    ci.upper.rho = difference.means + qnorm(1-alpha/2) * sd.rho
    ci.lower.rho = difference.means - qnorm(1-alpha/2) * sd.rho
  }
  
  ### Shortest CI:
  length.satc = ci.upper.treated-ci.lower.treated
  length.satt = ci.upper.control-ci.lower.control
  length.sate = ci.upper.neyman-ci.lower.neyman
  #length.sharp = ci.upper.sharp - ci.upper.sharp
  
  length.vec = c(length.sate,length.satc,length.satt)
  parameter.vec = c("SATE","SATC","SATT")
  
  length.shortestCI = length.vec[which.min(length.vec)]
  sd.shortestCI = c(sd.neyman,sd.adj.treated, sd.adj.control)[which.min(length.vec)]
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
  
  ### results to return
  results = list(
    
    neymanCI = list(ci.lower.neyman = ci.lower.neyman,
                    ci.upper.neyman = ci.upper.neyman,
                    tstat.neyman = difference.means/sd.neyman
                    ),
    
    sattCI = list(ci.lower.satt = ci.lower.control,
                  ci.upper.satt = ci.upper.control,
                  tstat.satt = difference.means/sd.adj.control
                  ),
    
    satcCI = list(ci.lower.satc = ci.lower.treated,
                  ci.upper.satc = ci.upper.treated,
                  tstat.satc = difference.means/sd.adj.treated
                  ),
    
    shortestCI = list(ci.lower = ci.lower.shortestCI,
                      ci.upper = ci.upper.shortestCI,
                      tstat.shortest = difference.means/sd.shortestCI,
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
  
  if(!is.null(rho)){
    results$sate.rho = list(ci.upper.sate.rho = ci.upper.rho,
                                     ci.lower.sate.rho = ci.lower.rho)
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
      " Variance of control units: ",round(var0,dig=2),
      "\n",
      " Variance of treated units: ",round(var1,dig=2),
      "\n",
      "\n \n \n",
      "Additional confidence intervals: \n \n",
      " SATE (Corr=1 CI) - ", (1-alpha)*100," percent confidence interval: \n", 
      "      [",ci.lower.corr1,", ",ci.upper.corr1,"]", 
      "\n",
        sep="")
    
    if(variance.test){
      cat( "\n", " F Test to Compare Two Variances:  P-value = ",var.test(y[tr==1],y[tr==0])$p.value,"\n" )  
    }
    
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































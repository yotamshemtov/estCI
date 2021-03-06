% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aveCI.R
\name{aveCI}
\alias{aveCI}
\title{Confidence intervals for SATT and SATC}
\usage{
aveCI(outcome = NULL, treatment = NULL, alpha = 0.05, print = TRUE,
  Sharp.CI = TRUE, variance.test = TRUE, stats.only = FALSE,
  var1 = NULL, var0 = NULL, mu1 = NULL, mu0 = NULL, n = NULL,
  m = NULL, var.pool = NULL, rho = NULL)
}
\arguments{
\item{outcome}{Outcome of interest \cr \cr}

\item{treatment}{Treatmnet indicator that can take only two levels 0 or 1. \cr \cr}

\item{alpha}{The Type-I error rate (size) of the confidence interval.}
}
\value{
The function returns TO-ADD help manual description
}
\description{
The funciton calculates confidence/prediction intervals (CI and PI) for sample average treatment effects. It calculates a prediction intervals for the sample average treatment effect on the treated (SATT), the sample average treatment effect on the controls (SATC), and a confidence intervals for the average treatment effect (SATE).
}
\examples{

### Example 1:
# Estimate the average treatment effect on the treated
n=1000
y0 = rnorm(n,mean=10,sd=1)
y1 = rnorm(n,mean=13,sd=3)
tau = y1-y0
tr = rep(0,n)
tr.index = sample(c(1:n),size=n/2,replace=FALSE)
tr[tr.index]=1
y = tr*y1 + y0*(1-tr)
results = aveCI(y,tr, print=TRUE)

### Example 2:
# Real data example using data from Tunca and Egeli (1996) that was also used by Rosenbaum (2001)
data(tunca.and.egeli.1996)
head(tunca.and.egeli.1996)
results = aveCI(outcome=tunca.and.egeli.1996$y,treatment=tunca.and.egeli.1996$tr, print=TRUE)

### Example 3:
# Real data examle using data from Lalonde (1986)
data(lalonde1986)
head(lalonde1986)
y = lalonde1986$re78
tr= lalonde1986$treat
# A confidence interval at the 5 level - default.
results = aveCI(y,tr, print=TRUE, Sharp.CI = TRUE)
# A confidence interval at the 1 level.
results = aveCI(y,tr, print=TRUE, alpha=0.01, Sharp.CI = TRUE)


}


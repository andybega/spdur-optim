################################################################
# Optimizing splolog in rcpp
if(Sys.info()['user']=='janus829' | Sys.info()['user']=='s7m'){
  setwd("~/Research/WardProjects/spduration-optim/ilc-example")}

# Load libraries
loadPkg=function(toLoad){
  for(lib in toLoad){
  if(! lib %in% installed.packages()[,1])
    { install.packages(lib, repos='http://cran.rstudio.com/') }
  suppressMessages( library(lib, character.only=TRUE) ) }
}

packs=c("corpcor", 'microbenchmark', 'Rcpp', 'RcppArmadillo', 'ggplot2')
loadPkg(packs)

# Load core estimation functions
core_funcs = list.files("../core-functions", pattern="*.R")
for(func in core_funcs){source(paste0('../core-functions/', func))}

# Look only at likelihood calculation
# (Y, X, Z saved from debugging spdur)
load("XYZ-matrices.rda")
set.seed(6886)
theta <- rnorm(ncol(X) + ncol(Z) + 1)
################################################################

################################################################
# Make sure cpp versions return same result
sourceCpp("../core-functions/weib_lnlC.cpp")
weib_lnlCpp(theta, Y, X); weib_lnl(theta, Y, X)

sourceCpp("../core-functions/spweib_lnlC.cpp")
spweib_lnlCpp(theta, Y, X, Z); spweib_lnl(theta, Y, X, Z)
################################################################

################################################################
# Speed comparison
# First we pass cpp likelihood into existing functions from spweibull.R
#' Regular Weibull regression
weibullCpp <- function(Y, X, inits=NULL, max.iter, silent=TRUE) {
  if (is.null(inits)) {
    inits <- c(rep(0, ncol(X)), 0)
  }
  
  trace <- !silent
  est <- optim(inits, weib_lnlCpp, method="BFGS", 
               control=list(trace=trace, maxit=max.iter), 
               hessian=T, y=Y, X=X)
  
  # Solve other results
  if (est$convergence!=0 & !silent) stop('Model did not converge')
  coef <- est$par
  vcv <- solve(est$hessian)
  vcv <- make.positive.definite(vcv)
  logL <- -est$value
  
  # Put together results
  return(list(coefficients = coef, vcv = vcv, logL = logL))
}

# Split Population
spweibullCpp <- function(Y, X, Z, max.iter, silent=FALSE) {  
  # Estimate base model
  if (!exists("base.inits")) {
    base.inits <- c(rep(0, ncol(X)), 0)
  }
  if (!silent) cat('Fitting base weibull...\n')
  base <- weibull(Y=Y, X=X, inits=base.inits, max.iter=200, silent=TRUE)
  
  # Estimate full model
  x.inits <- base$coefficients[1:ncol(X)]
  a.init <- base$coefficients[ncol(X)+1]
  if (!silent) cat('Fitting split weibull...\n')
  trace <- !silent
  est <- optim(c(x.inits, rep(0, ncol(Z)), a.init), spweib_lnlCpp, method="BFGS", 
    control=list(trace=trace, maxit=max.iter), hessian=T, y=Y, X=X, Z=Z)
  
  # Solve other results
  if (est$convergence!=0 & !silent) stop('Model did not converge')
  coef <- est$par
  vcv <- solve(est$hessian)
  vcv <- make.positive.definite(vcv)
  logL <- -est$value
  
  # Put together results
  return(list(coefficients = coef, vcv = vcv, logL = logL, base=base))
}

if(!'weibullTimes.rda' %in% list.files()){
  times=microbenchmark(
    weibull(Y, X, max.iter=200, silent=TRUE),
    weibullCpp(Y, X, max.iter=200, silent=TRUE),
    spweibull(Y, X, Z, max.iter=200, silent=TRUE),
    spweibullCpp(Y, X, Z, max.iter=200, silent=TRUE), 
    times=100, unit="s")
  save(times, file='weibullTimes.rda')
}

if(!'times' %in% ls()){ load('weibullTimes.rda') }
print(times)
################################################################
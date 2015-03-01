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
sourceCpp("../core-functions/loglog_lnlC.cpp")
loglog_lnlCpp(theta, Y, X); loglog_lnl(theta, Y, X)

sourceCpp("../core-functions/sploglog_lnlC.cpp")
sploglog_lnlCpp(theta, Y, X, Z); sploglog_lnl(theta, Y, X, Z)
################################################################

################################################################
# Speed comparison
# First we pass cpp likelihood into existing functions from sploglog.R
#' Regular Log-logistic regression
loglogCpp <- function(Y, X, inits=NULL, max.iter, silent=TRUE) {
  if (is.null(inits)) {
    inits <- c(rep(0, ncol(X)), 0)
  }
  
  trace <- !silent
  est <- optim(inits, loglog_lnlCpp, method="BFGS", 
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
sploglogCpp <- function(Y, X, Z, max.iter, silent=FALSE) {
  # Estimate base model
  if (!exists("base.inits")) {
    base.inits <- c(rep(0, ncol(X)), 1)
  }
  if (!silent) cat('Fitting base loglog...\n')
  base <- loglog(Y=Y, X=X, inits=base.inits, max.iter=200, silent=TRUE)
    
  # Estimate full model
  x.inits <- base$coefficients[1:ncol(X)]
  a.init <- base$coefficients[ncol(X)+1]
  if (!silent) cat('Fitting split loglog...\n')
  trace <- !silent
  est <- optim(c(x.inits, rep(0, ncol(Z)), a.init), sploglog_lnlCpp, method="BFGS", 
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

if(!'loglogTimes.rda' %in% list.files()){
  times=microbenchmark(
    loglog(Y, X, max.iter=200, silent=TRUE),
    loglogCpp(Y, X, max.iter=200, silent=TRUE),
    sploglog(Y, X, Z, max.iter=200, silent=TRUE),
    sploglogCpp(Y, X, Z, max.iter=200, silent=TRUE), 
    times=100, unit="s")
  save(times, file='loglogTimes.rda')
}

if(!'times' %in% ls()){ load('loglogTimes.rda') }
print(times)
################################################################
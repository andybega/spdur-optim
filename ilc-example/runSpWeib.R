################################################################
# Optimizing spweibull in rcpp
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
coreFuncs = list.files("../core-functions", pattern="*.R")
for(func in coreFuncs){source(paste0('../core-functions/', func))}

cppFuncs = list.files("../core-functions", pattern="*.cpp")
for(func in cppFuncs){sourceCpp(paste0('../core-functions/', func))}

# Look only at likelihood calculation
# (Y, X, Z saved from debugging spdur)
load("XYZ-matrices.rda")
set.seed(6886)
theta <- rnorm(ncol(X) + ncol(Z) + 1)
################################################################

################################################################
# Make sure cpp versions return same result
weib_lnlCpp(theta, Y, X); weib_lnl(theta, Y, X)
spweib_lnlCpp(theta, Y, X, Z); spweib_lnl(theta, Y, X, Z)
################################################################

################################################################
# Speed comparison
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
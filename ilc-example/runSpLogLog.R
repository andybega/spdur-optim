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
loglog_lnlCpp(theta, Y, X); loglog_lnl(theta, Y, X)
sploglog_lnlCpp(theta, Y, X, Z); sploglog_lnl(theta, Y, X, Z)
################################################################

################################################################
# Speed comparison
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
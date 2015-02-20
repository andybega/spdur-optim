# Optimizing splolog in rcpp

if(Sys.info()['user']=='janus829' | Sys.info()['user']=='s7m'){
  setwd("~/Research/WardProjects/spduration-optim/ilc-example")}

loadPkg=function(toLoad){
  for(lib in toLoad){
  if(! lib %in% installed.packages()[,1])
    { install.packages(lib, repos='http://cran.rstudio.com/') }
  suppressMessages( library(lib, character.only=TRUE) ) }
}

# Load libraries
packs=c("corpcor", 'lineprof', 'microbenchmark')
loadPkg(packs)

# Load core estimation functions
core_funcs <- list.files("../core-functions", pattern="*.R")
for(func in core_funcs){source(paste0('../core-functions/', func))}

# Load example data
load("irc-data-mod.rda")

# Theme 1 with all data
model1 <- function() {
  spdur(duration ~ 1 + log10(i.matl.conf.DIStGOV.l1+1) + 
          log10(i.matl.coop.GOVtGOV.l1+1),
        atrisk ~ 1 + ldr.irregular + ldr.foreign + log10(mths.in.power+1),
        data=irc.data, silent=TRUE)
}

# Run time; use times>1 to eliminate overhead variance
# microbenchmark(model1(), times=1, unit="s")

# # Line by line profiling; might give different results each run
# m1_l1 <- lineprof({
#   spdur(duration ~ 1 + log10(i.matl.conf.DIStGOV.l1+1) + 
#           log10(i.matl.coop.GOVtGOV.l1+1),
#         atrisk ~ 1 + ldr.irregular + ldr.foreign + log10(mths.in.power+1),
#         data=irc.data, silent=TRUE)
# }, interval=0.1)

# shine(m1_l1)


# Look only at likelihood calculation
# (Y, X, Z saved from debugging spdur)
load("XYZ-matrices.rda")

# weib_lnl
# spweib_lnl

set.seed(6886)
theta <- rnorm(ncol(X) + ncol(Z) + 1)
# spweib_l1 <- lineprof(spweib_lnl(theta, Y, X, Z), torture=T)

# Testing out stuff
library(Rcpp)
library(RcppArmadillo)
library(rbenchmark)

sourceCpp("../core-functions/loglog_lnlC.cpp")

loglog_lnlCpp(theta, Y, X); loglog_lnl(theta, Y, X)

sourceCpp("../core-functions/sploglog_lnlC.cpp")

sploglog_lnlCpp(theta, Y, X, Z); sploglog_lnl(theta, Y, X, Z)